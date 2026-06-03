library(dHRUM)
library(data.table)
library(SoilHyP)   # provides SCEoptim()
library(hydroGOF)  # provides NSE()

# ── Model setup ───────────────────────────────────────────────────────────────
nHrus  <- 1
Areas  <- 4.7 * 1000 * 1000
IdsHrus <- paste0("ID", seq_along(Areas))
dhrus  <- initdHruModel(nHrus, Areas, IdsHrus)

setGWtypeToAlldHrus(dHRUM_ptr = dhrus,
                    gwTypes = rep("LIN_RES", length(Areas)), hruIds = IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,
                          soilTypes = rep("PDM", length(Areas)), hruIds = IdsHrus)
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,
                             surfaceStorTypes = rep("SurfaceAll", length(Areas)), hruIds = IdsHrus)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,
                           fastResponseTypes = rep("SerialCascadeLinRes", length(Areas)), hruIds = IdsHrus)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("Rutter_Gash",times= length(Areas)),hruIds=IdsHrus,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep("Pitman",times= length(Areas)))


filname2 <- "./Calibrations/Amalie/indata/BP_1960_01_01.txt"
setPTInputsToAlldHrusFromFile(dHRUM_ptr = dhrus, filname2)
calcPetToAllHrus(dHRUM_ptr = dhrus, 50.1, "HAMON")

# ── Initial parameters ────────────────────────────────────────────────────────
ParDF <- data.frame(
  B_SOIL = 1.6, C_MAX = 35,   B_EVAP = 2.5,  KS = 0.01,  KF = 0.03,
  ADIV = 0.8,   CDIV = 0.2,   SDIV = 0.1,    CAN_ST = 2, STEM_ST = 1,
  CSDIV = 0.8,  TETR = 0,     DDFA = 0.75,   TMEL = 0.0, RETCAP = 10,
  D_BYPASS = 0.8, THR = 10,   KS2 = 0.1,     ALPHA = 0.5, FOREST_FRACT = 0.3,
  FC = 10,      KF_NONLIN = 10, KF2 = 0.01,  C = 10,     INFR_MAX = 10,
  RF = 0.5,     WP = 0.3,     CMIN = 25,     L = 0.1,    B_EXP = 0.3,
  KFR = 0.03
)

setParamsToAlldHrus(dHRUM_ptr = dhrus,
                    ParsVec = as.numeric(ParDF[1, ]), ParsNames = names(ParDF))
myPars <- getCurdHRUpars(dHRUM_ptr = dhrus, 0)

# myPars is a list:
#   [[1]] Cur_names – parameter names
#   [[2]] Cur_vals  – current values
#   [[3]] Upper     – upper bounds
#   [[4]] Lower     – lower bounds
parNames <- myPars$Cur_names
parLower <- as.numeric(myPars[[4]])
parUpper <- as.numeric(myPars[[3]])
parInit  <- as.numeric(myPars[[2]])

# ── Observed data ─────────────────────────────────────────────────────────────
# FDC targets (13 long-term quantile points)
dny   <- c(30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 355, 364)
p_OBS <- dny / 365.25
A     <- 4.7 * 1000 * 1000                   # catchment area [m²]
QmBP  <- c(26, 18, 14, 12, 10, 8.0, 7.0,
           6.0, 4.5, 3.5, 2.5, 1.0, 0.5)    # [l/s per day]
RmBP  <- QmBP * (3600 * 24) / A             # [mm/day]

# Observed daily discharge time series [mm/day], same length as model output.
# Set to NULL to fall back to FDC-only calibration.
#   obs_Q <- read.table("path/to/observed_Q_mm.txt")[[1]]
obs_Q <- NULL

# ── Objective function factory ────────────────────────────────────────────────
# Returns a ready-to-use function compatible with SCEoptim (single numeric arg).
#
# mode         : "fdc"      – SSE on flow duration curve only
#                "nse"      – (1 – NSE) on daily time series only
#                "combined" – weighted sum of both (requires obs_Q)
# obs_Q        : numeric vector of observed daily runoff [mm/day].
#                Must be the same length as the model output. Required for
#                "nse" and "combined" modes.
# warmup       : number of leading time steps to skip when computing NSE
#                (spin-up period, default 365 days).
# w_nse        : weight given to the NSE component in "combined" mode
#                (w_fdc = 1 - w_nse). Ignored in single-metric modes.
# log_transform: if TRUE, apply log() to sim/obs before computing NSE
#                (emphasises low-flow fit).

make_objective <- function(mode        = c("fdc", "nse", "combined"),
                           obs_Q       = NULL,
                           warmup      = 365L,
                           w_nse       = 0.5,
                           log_transform = FALSE) {

  mode <- match.arg(mode)

  if (mode %in% c("nse", "combined") && is.null(obs_Q))
    stop("obs_Q must be supplied for mode = '", mode, "'")

  # Pre-compute reference variance of the FDC targets (for normalisation).
  fdc_ref <- sum(RmBP^2)

  function(myPar) {
    setParamsToAlldHrus(dHRUM_ptr = dhrus,
                        ParsVec   = as.numeric(myPar),
                        ParsNames = parNames)
    calcHBInAlldHrus(dHRUM_ptr = dhrus)
    gatherHBdata(dHRUM_ptr = dhrus)

    dta <- getOutput(dHRUM_ptr = dhrus)
    dF  <- data.frame(dta$outDta)
    names(dF) <- dta$VarsNams
    sim <- dF$TOTR

    # ── FDC component: normalised SSE (dimensionless) ─────────────────────
    fdc_cost <- if (mode %in% c("fdc", "combined")) {
      simRM <- as.numeric(quantile(sim, probs = (1 - p_OBS), na.rm = TRUE))
      sum((simRM - RmBP)^2) / fdc_ref
    } else 0

    # ── NSE component: 1 – NSE (dimensionless, 0 = perfect) ──────────────
    nse_cost <- if (mode %in% c("nse", "combined")) {
      idx   <- seq(warmup + 1L, length(sim))
      s     <- sim[idx]
      o     <- obs_Q[idx]
      valid <- !is.na(o) & !is.na(s) & is.finite(o) & is.finite(s)
      s <- s[valid]; o <- o[valid]

      if (log_transform) {
        s <- log(pmax(s, 1e-6))
        o <- log(pmax(o, 1e-6))
      }

      ss_res <- sum((s - o)^2)
      ss_tot <- sum((o - mean(o))^2)
      # cap at 10 to avoid -Inf penalties dominating the search
      min(ss_res / ss_tot, 10)
    } else 0

    # ── Combined cost ─────────────────────────────────────────────────────
    w_fdc <- 1 - w_nse
    as.double(w_fdc * fdc_cost + w_nse * nse_cost)
  }
}

# ── Build the objective ───────────────────────────────────────────────────────
# Change `mode` to "nse" or "combined" once obs_Q is loaded.
obj <- make_objective(
  mode          = "fdc",     # "fdc" | "nse" | "combined"
  obs_Q         = obs_Q,
  warmup        = 365L,
  w_nse         = 0.5,       # only used in "combined" mode
  log_transform = FALSE
)

# ── SCE-UA calibration ────────────────────────────────────────────────────────
sce_control <- list(
  ncomplex = 5,      # number of complexes (increase for harder problems)
  maxit    = 10000,  # max SCE iterations
  maxeval  = Inf,    # max function evaluations (Inf = unlimited)
  reltol   = 1e-5,
  tolsteps = 20,     # consecutive iterations within reltol to confirm convergence
  trace    = 1
)

set.seed(4721)
u <- SCEoptim(
  FUN     = obj,
  par     = parInit,
  lower   = parLower,
  upper   = parUpper,
  control = sce_control
)

cat("Best cost:", u$value, "\n")

# ── Apply best parameters and run final simulation ────────────────────────────
setParamsToAlldHrus(dHRUM_ptr = dhrus,
                    ParsVec   = as.numeric(u$par),
                    ParsNames = parNames)
calcHBInAlldHrus(dHRUM_ptr = dhrus)
gatherHBdata(dHRUM_ptr = dhrus)

dta <- getOutput(dHRUM_ptr = dhrus)
dF  <- data.frame(dta$outDta)
names(dF) <- dta$VarsNams

# ── Diagnostic plots ──────────────────────────────────────────────────────────
simBest <- as.numeric(quantile(dF$TOTR, probs = (1 - p_OBS), na.rm = TRUE))

plot(p_OBS, RmBP, pch = 19,
     xlab = "Exceedance probability", ylab = "Runoff [mm/day]",
     main = "Flow Duration Curve – BP catchment")
points(p_OBS, simBest, col = "red", pch = 19)
legend("topright", legend = c("Observed", "Simulated (SCE-UA)"),
       col = c("black", "red"), pch = 19)

if (!is.null(obs_Q)) {
  nt  <- length(dF$TOTR)
  nse_val <- NSE(sim = dF$TOTR[366:nt], obs = obs_Q[366:nt])
  cat("NSE (post-warmup):", round(nse_val, 3), "\n")
  plot(obs_Q, type = "l", col = "black",
       xlab = "Day", ylab = "Runoff [mm/day]", main = paste("NSE =", round(nse_val, 3)))
  lines(dF$TOTR, col = "red")
  legend("topright", legend = c("Observed", "Simulated"), col = c("black", "red"), lty = 1)
}

# ── Save results ──────────────────────────────────────────────────────────────
BP_SCE_df        <- as.data.frame(t(u$par))
names(BP_SCE_df) <- parNames[seq_along(u$par)]
BP_SCE_df$ID     <- "BP"

# save(BP_SCE_df, file = "./Calibrations/Amalie/outdata/par_dHRUM_lumped_SCE_BP.rda")
BP_SCE_df

