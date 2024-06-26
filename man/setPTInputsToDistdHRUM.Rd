% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{setPTInputsToDistdHRUM}
\alias{setPTInputsToDistdHRUM}
\title{Sets the distributed Precipitation, Temperature vectors to distributed dHRUM and init's the date using beg. of period.}
\usage{
setPTInputsToDistdHRUM(dHRUM_ptr, DataDF)
}
\arguments{
\item{dHRUM_ptr}{pointer to dHRUM instance}

\item{DataDF}{dataframe with DTM, Precipitation, Temperature, and HRU Ids}
}
\description{
Setting the different vector of Precipitation and temperature to all single HRU.
Setting the calender using the first date fo period using the first date of period
}
\examples{
nHrus <- 2
Areas <- runif(nHrus,min = 1,max  = 10)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
setGWtypeToAlldHrus(dHRUM_ptr = dhrus,gwTypes=rep("LIN_2SE",times= length(Areas)),hruIds=IdsHrus)
setSoilStorTypeToAlldHrus(dHRUM_ptr = dhrus,soilTypes=rep("PDM",times= length(Areas)),hruIds=IdsHrus)
}
