#' Encode Charge State to One-Hot Tensor
#'
#' @param uniqueDt A data.table with unique peptides and a column called "Charge"
#'
#' @return
#' @export
#'
#' @examples
numeric2oneHot <- function(uniqueDt){
  cnames <- colnames(uniqueDt)
  colMissing <- !("Charge" %in% cnames)
  if(colMissing){
    stop("Column \'Charge\' is missing")
  }
  
  m <- matrix(nrow = length(uniqueDt$Charge), ncol = max(uniqueDt$Charge), data = 0L)
  
  #Fill Data
  locs <- cbind(c(1:length(uniqueDt$Charge)), uniqueDt$Charge)
  m[locs] <- 1
  
  m
}
