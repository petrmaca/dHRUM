% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{printToFile}
\alias{printToFile}
\title{Write Output of overall dHru simulation to file}
\usage{
printToFile(dHRUM_ptr, namOutFilet)
}
\arguments{
\item{dHRUM_ptr}{pointer to dHRU instance}

\item{namOutFilet}{file and path to write output to}
}
\description{
Write Output of overall dHru simulation to file
}
\examples{
nHrus <- 200
Areas <- runif(nHrus,min = 1,max  = 10)
IdsHrus <- paste0("ID",seq(1:length(Areas)))
dhrus <- initdHruModel(nHrus,Areas,IdsHrus)
filname2 = "../data/inBP_1960_01_01.txt"
setInputsToAlldHrus(dhrus, filname2)
ParDF = data.frame( B_SOIL = 1.6, C_MAX = 100, B_EVAP = 2,  KS = 0.1, KF = 0.2, ADIV = 0.3, CDIV = 0.03,
SDIV = 0.03, CAN_ST = 2, STEM_ST = 2, CSDIV = 0.3, TETR = 5, DDFA = 0.5, TMEL = 0, RETCAP = 10 )
setParamsToAlldHrus(dhrus,as.numeric(ParDF[1,]),names(ParDF),TRUE,0)
calcPetToHrus(dHRUM_ptr = dhrus,50.1,"Hamon")
calcHBInAlldHrus(dHRUM_ptr = dhrus)
gatherHBdata(dHRUM_ptr = dhrus)
printToFile(dhrus,file)
}
