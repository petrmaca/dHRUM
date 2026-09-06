# tests/pet_regression.r
#
# Regression checks for the six PET methods in dHRUM.
#
# Each method is run on a single HRU with synthetic two-year daily inputs
# (spanning leap year 2000 and a year boundary), and the computed PET series
# is compared against an R transcription of the formula implemented in
# src/data_HB_1d.cpp. The script calls stop() on any mismatch.
#
# Usage: Rscript tests/pet_regression.r   (or source() from an R session;
# requires the dHRUM package to be installed from this source tree)
#
# Note: the Thornthwaite transcription deliberately reproduces the current
# C++ behavior, including the quirk that the final year's annual heat index
# sums only Jan--Nov (see REFACTORING_PLAN.md). Fixing that quirk requires
# updating this test.

library(dHRUM)

# ---------- synthetic input ----------------------------------------------

dates <- seq(as.Date("2000-01-01"), by = "day", length.out = 731)
jday <- as.integer(format(dates, "%j"))
yr   <- as.integer(format(dates, "%Y"))
mon  <- as.integer(format(dates, "%m"))
temp <- 8 + 12 * sin(2 * pi * (jday - 90) / 365.25) # seasonal; < 0 in winter
lat  <- 50

run_pet <- function(method) {
  d <- initdHruModel(1, 10, "ID1")
  setGWtypeToAlldHrus(dHRUM_ptr = d, gwTypes = "LIN_RES", hruIds = "ID1")
  setSoilStorTypeToAlldHrus(dHRUM_ptr = d, soilTypes = "PDM", hruIds = "ID1")
  setPTDateInputsToAlldHrus(d, Prec = rep(0, length(dates)), Temp = temp,
                            DateVec = dates)
  calcPetToAllHrus(dHRUM_ptr = d, lat, method)
  od <- dHRUMrunDist(d)
  od$outDta[, which(od$VarsNams == "PET")[1]]
}

# ---------- R transcriptions of src/data_HB_1d.cpp ------------------------

radlat <- lat / 180 * pi
leap <- (yr %% 4 == 0 & yr %% 100 != 0) | yr %% 400 == 0
ndy  <- ifelse(leap, 366, 365)
dec  <- 0.409 * sin(2 * jday * pi / ndy - 1.39)
dr   <- 1 + 0.033 * cos(jday * 2 * pi / ndy)
om   <- acos(-tan(radlat) * tan(dec))
Ra   <- (24 * 60) / pi * 0.0820 * dr *
        (om * sin(radlat) * sin(dec) + cos(radlat) * cos(dec) * sin(om))

expected <- list(
  OUDIN = ifelse(temp + 5 >= 0, 0.408 * Ra * (temp + 5) / 100, 0),
  HAMON = 0.1651 * 216.7 * (24 / pi * om / 12) *
          (6.108 * exp(17.27 * temp / (temp + 237.3)) / (temp + 273.3)),
  BLANEYCRIDDLE = (24 / pi * om / (365 * 12) * 0.85) * 100 * (0.46 * temp + 8.13),
  JENSENHAISE = pmax(1000 * Ra * temp / (40 * 2450), 0),
  MCGUINNESSBORDNE = pmax(1000 * Ra * (temp + 5) / (68 * 2450), 0)
)

# Faithful transcription of data_HB_1d::ThornthwaitePET, including its
# current-behavior quirks (see file header note).
thorne_expected <- function(temp, year, month, jday, lat) {
  n <- length(temp)
  radlat <- lat / 180 * pi
  leap <- (year %% 4 == 0 & year %% 100 != 0) | year %% 400 == 0
  ndy <- ifelse(leap, 366, 365)
  dec <- 0.409 * sin((2 * pi) * jday / ndy - 1.39)
  omega <- acos(-tan(radlat) * tan(dec))
  Nn <- 24 / pi * omega

  nmonths <- sum(c(1, diff(month) != 0)) + 1
  nyears  <- sum(c(1, diff(year) != 0)) + 1
  tam <- numeric(nmonths); Nmeanmonth <- numeric(nmonths)
  numDaysMonth <- numeric(nmonths); helpyear <- numeric(nmonths)

  helptam <- temp[1]; helpNn <- Nn[1]; nd <- 1; hi <- 1
  for (tst in 2:n) {
    if (month[tst] > month[tst - 1] || year[tst] > year[tst - 1]) {
      helptam <- 0; helpNn <- 0; hi <- hi + 1; nd <- 0
    }
    helptam <- helptam + temp[tst]; helpNn <- helpNn + Nn[tst]; nd <- nd + 1
    tam[hi] <- helptam / nd; Nmeanmonth[hi] <- helpNn / nd
    numDaysMonth[hi] <- nd; helpyear[hi] <- year[tst]
  }

  i_heat <- numeric(nmonths)
  for (it in seq_len(nmonths)) {
    if (tam[it] < 0) {
      tam[it] <- 0                     # note: mutates tam, as in C++
    } else if (tam[it] > 0) {
      i_heat[it] <- (tam[it] * 0.2)^1.514
    }
  }

  I_annual <- numeric(nyears); acoeff <- numeric(nyears)
  helpsumI <- i_heat[1]; helpIn <- 1
  for (it in 2:nmonths) {
    helpsumI <- helpsumI + i_heat[it]
    I_annual[helpIn] <- helpsumI - i_heat[it]
    Ih <- I_annual[helpIn]
    acoeff[helpIn] <- 6.751e-07 * Ih^3 - 7.71e-05 * Ih^2 + 0.01729 * Ih + 0.49239
    if (helpyear[it] != helpyear[it - 1]) { helpIn <- helpIn + 1; helpsumI <- i_heat[it] }
  }

  amonthly <- numeric(nmonths); Imonthly <- numeric(nmonths)
  Epetraw <- numeric(nmonths)
  amonthly[1] <- acoeff[1]; Imonthly[1] <- I_annual[1]
  Epetraw[1] <- 16 * (10 * tam[1] / Imonthly[1])^amonthly[1]
  hyInd <- 1
  for (it in 2:nmonths) {
    if (helpyear[it] != helpyear[it - 1]) hyInd <- hyInd + 1
    amonthly[it] <- acoeff[hyInd]; Imonthly[it] <- I_annual[hyInd]
    Epetraw[it] <- 16 * (10 * tam[it] / Imonthly[it])^amonthly[it]
  }
  Epet <- Epetraw * (Nmeanmonth / 12) * (numDaysMonth / 30)

  PEt <- numeric(n)
  helpit <- 1
  PEt[1] <- Epet[1] / numDaysMonth[1]
  for (tst in 2:n) {
    if (month[tst] != month[tst - 1]) helpit <- helpit + 1
    PEt[tst] <- Epet[helpit] / numDaysMonth[helpit]
  }
  PEt
}
expected$THORNTHWAITE <- thorne_expected(temp, yr, mon, jday, lat)

# ---------- compare --------------------------------------------------------

tol <- 1e-9
results <- data.frame(
  method = names(expected),
  max_abs_diff = NA_real_,
  day1_computed = NA_real_,
  day1_expected = NA_real_,
  status = "FAIL",
  row.names = NULL
)
for (i in seq_along(expected)) {
  m <- names(expected)[i]
  got <- run_pet(m)
  results$max_abs_diff[i] <- max(abs(got - expected[[m]]))
  results$day1_computed[i] <- got[1]
  results$day1_expected[i] <- expected[[m]][1]
  results$status[i] <- if (results$max_abs_diff[i] < tol) "PASS" else "FAIL"
}

print(format(results, digits = 8), right = FALSE)

if (any(results$status == "FAIL")) {
  stop("PET regression test FAILED for: ",
       paste(results$method[results$status == "FAIL"], collapse = ", "))
}
cat("\nPET regression test: all methods PASS (tol =", tol, ")\n")
