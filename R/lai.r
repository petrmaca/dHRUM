#' Read Lai Data from Package for Ohre river basin
#'
#' This function locates and reads the sample LAI file file stored within the
#' package's internal data directory.
#'
#' @return A data frame containing the contents of 'LAI_POH_Category.rds' stored in .
#'
#' @export
#'
#' @examples
#' \dontrun{
#' laiDT <- get_Lai_DataCat()
#' head(laiDT)
#' }
get_Lai_DataCat <- function() {
  # This finds the actual path on the user's computer
  path <- system.file("extdata", "LAI_POH_Category.rds", package = "dHRUM")

  # Safety check: if the file isn't found, system.file returns ""
  if (path == "") {
    stop("File not found in the package!")
  }

  readRDS(path)
}

#' Read Lai Data from Package for Ohre river basin
#'
#' This function locates and reads the sample LAI file file stored within the
#' package's internal data directory.
#'
#' @return A data frame containing the contents of 'LAI_POH_Category.rds' stored in .
#'
#' @export
#'
#' @examples
#' \dontrun{
#' laiDT <- get_Lai_DataPOh()
#' head(laiDT)
#' }
get_Lai_DataPOh <- function() {
  # This finds the actual path on the user's computer
  path <- system.file("extdata", "LAI_POH.rds", package = "dHRUM")

  # Safety check: if the file isn't found, system.file returns ""
  if (path == "") {
    stop("File not found in the package!")
  }

  readRDS(path)
}
