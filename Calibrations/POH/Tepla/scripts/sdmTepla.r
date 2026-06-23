library(dHRUM)
library(data.table)
library(fst)
library(lubridate)
library(SoilHyP) # SCE-UA

Qmd <-  as.data.table(read.fst("/home/hubert/MyData/POH/data/rmvody/mQ_4chpPlus.fst"))
dta <- as.data.table(readRDS("/home/hubert/MyData/POH/data/Tepla/IN/INP_Tepla.rds"))
ddd <- as.data.table(readRDS("/home/hubert/MyData/POH/data/Tepla/IN/INF_Tepla.rds"))
ddd[,ID := chp_14_s]

teplaCHP <- as.data.table(read.table("/home/hubert/MyData/POH/data/OUT/catchments/Tepla_chp14s.txt",header = 1))

chps_Brezova <- as.data.table(teplaCHP[Povodi %in% "Brezova",])
chps_Brezova[, ID := chp14s]

resBrezova <- merge(dta,chps_Brezova, by = "ID", all.y = TRUE )

resBrezovaII <- merge(ddd,resBrezova, by = "ID", all.y = TRUE )

length(resBrezovaII[, unique(ID)])
nrow(chps_Brezova)

Days <- c(30,60,90,120,150,180,210,240,270,300,330,355,364)

nHrus <- 1
ntreads <- 1
IdsHrus <- paste0(basinCHP)
Areas <- resBrezovaII$a_km2[i]
latBas <- resBrezovaII$LAT[i]

i=1
basinCHP <- as.character(resBrezovaII[i,ID])
dta1basin <- dta[ID %in% basinCHP,]



dhrus <- initdHruModel(nHrus,Areas,IdsHrus,ntreads)

setSnowMeltModeltypeToAlldHrus(dHRUM_ptr = dhrus,snowMeltModelTypes = rep("DDF",times = length(Areas)),hruIds = IdsHrus)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("van_Dijk",times= length(Areas)),hruIds=IdsHrus,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep("Pitman",times= length(Areas)))
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep("SurfaceAll",times= length(Areas)),hruIds=IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, soilTypes = rep("PDM",times= length(Areas)), hruIds = IdsHrus)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep("SerialCascadeLinRes",times= length(Areas)),hruIds=IdsHrus)
setGWtypeToAlldHrus(dhrus,gwTypes = rep("LIN_RES", times =nHrus), IdsHrus)
setPTLInputsToAlldHrus(dhrus, Prec = dta1basin$PREP, Temp = dta1basin$TAVG, Lai = dta1basin$LAI, inDate =  as.Date("1961/01/01"))
calcPetToAllHrus(dHRUM_ptr = dhrus,latBas,"HAMON")
pars <- as.data.table(getCurdHRUpars(dHRUM_ptr = dhrus, singleHruId = 0))
pars = pars[-c(6,8,10,12,13,14),]
parscor = pars

parNames <- parscor$Cur_names
parLower <- as.numeric(parscor[[4]])
parUpper <- as.numeric(parscor[[3]])
parInit  <- as.numeric(parscor[[2]])

RMobs = as.numeric(Qmd[chp_4 %in% basinCHP,][,2:15])

dny   <- c(30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 355, 364)
p_OBS <- dny / 365.25
# A     <- 4.7 * 1000 * 1000                   # catchment area [m²]
# QmBP  <- c(26, 18, 14, 12, 10, 8.0, 7.0,
#            6.0, 4.5, 3.5, 2.5, 1.0, 0.5)    # [l/s per day]
# RmBP  <- QmBP * (3600 * 24) / A             # [mm/day]

# Observed daily discharge time series [mm/day], same length as model output.
# Set to NULL to fall back to FDC-only calibration.
# obs_Q <- read.table("path/to/observed_Q_mm.txt")[[1]]
obs_Q <- NULL


make_objective <- function(mode        = c("fdc", "nse", "combined"),
                           obs_Q       = NULL,
                           warmup      = 365L,
                           w_nse       = 0.5,
                           log_transform = FALSE) {

  mode <- match.arg(mode)

  if (mode %in% c("nse", "combined") && is.null(obs_Q))
    stop("obs_Q must be supplied for mode = '", mode, "'")

  # Pre-compute reference variance of the FDC targets (for normalisation).
  fdc_ref <- sum(RMobs^2)

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
      simRM <- c(mean(sim), simRM)
      sum((simRM - RMobs)^2) / fdc_ref
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

setParamsToAlldHrus(dHRUM_ptr = dhrus,
                    ParsVec   = as.numeric(u$par),
                    ParsNames = parNames)
calcHBInAlldHrus(dHRUM_ptr = dhrus)
gatherHBdata(dHRUM_ptr = dhrus)

dtaDhrum <- getOutput(dHRUM_ptr = dhrus)
dF  <- data.frame(dtaDhrum$outDta)
names(dF) <- dtaDhrum$VarsNams

# ── Diagnostic plots ──────────────────────────────────────────────────────────
simBest <- as.numeric(quantile(dF$TOTR, probs = (1 - p_OBS), na.rm = TRUE))

plot(p_OBS, RMobs[2:14], pch = 19,
     xlab = "Exceedance probability", ylab = "Runoff [mm/day]",
     main = paste0("Flow Duration Curve – ", basinCHP ,"catchment"))
points(p_OBS, simBest, col = "red", pch = 19)
legend("topright", legend = c("Observed", "Simulated (SCE-UA)"),
       col = c("black", "red"), pch = 19)

plot(dF$TOTR, type="l",ylab="R [mm/day]", xlab ="T [day]")
lines(dF$BASF, col="red")

j=27
nHrus <- j
ntreads <- j
Areas <- resBrezovaII$a_km2[1:j]
latBas <- resBrezovaII$LAT[1:j]

basinCHPs <- unique(chps_Brezova[,ID])[1:j]
dta2basins <- dta[ID %in% basinCHPs,]

# dtaIpns <- data.table(HruId = dta2basins$ID, P =dta2basins$PREP, T =dta2basins$TAVG, LAI = dta2basins$LAI, DTM = dta2basins$DTM)
dtaIpns <- data.table(HruId = dta2basins$ID, P =dta2basins$PREP, T =dta2basins$TAVG, LAI = dta2basins$LAI, DTM = as.Date(dta2basins$DTM))
Areas <- ddd$a_km2[1:j]
dhrus <- initdHruModel(nHrus,Areas,basinCHPs,ntreads)
setSnowMeltModeltypeToAlldHrus(dHRUM_ptr = dhrus,snowMeltModelTypes = rep("DDF",times = length(Areas)),hruIds = basinCHPs)
setInterceptiontypeToAlldHrus(dHRUM_ptr = dhrus,intcptnTypes=rep("van_Dijk",times= length(Areas)),hruIds=basinCHPs,InstStLai = rep(TRUE,times= length(Areas)),smaxlaiTypes = rep("Pitman",times= length(Areas)))
setSurfaceStortypeToAlldHrus(dHRUM_ptr = dhrus,surfaceStorTypes=rep("SurfaceAll",times= length(Areas)),hruIds=basinCHPs)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus, soilTypes = rep("PDM",times= length(Areas)), hruIds = basinCHPs)
setFastResponsesToAlldHrus(dHRUM_ptr = dhrus,fastResponseTypes=rep("SerialCascadeLinRes",times= length(Areas)),hruIds=basinCHPs)
setGWtypeToAlldHrus(dhrus,gwTypes = rep("LIN_RES", times =nHrus), basinCHPs)
setPTLInputsToDistdHRUM(dhrus, dtaIpns)
calcPetToAllHrusDist(dHRUM_ptr = dhrus,latBas,rep("HAMON", times =nHrus),basinCHPs)



Dta<-dHRUMrunDist(dHRUM_ptr = dhrus)
dF  <- data.frame(Dta$outDta)
dF <- cbind(dF,Dta$Ids)
names(dF) <- c(Dta$VarsNams)
dF = as.data.table(dF)




