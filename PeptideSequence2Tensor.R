#' One-Hot Tensor Encoding of Sequence Data
#'
#' @param uniqueDt a data.table containing the sequence "Name", the charge
#'   "Charge", and an index specifying unique sequence-charge combinations
#'   "name_charge_index" (one sequence can have multiple charge states)
#' @param aa a character vector of amino acid codes
#'
#' @return an array
#' @export
#' 
 
sequence2tensor_oneHot <- function(uniqueDt, aa){
  #number of sequences to process
  nSequences <- max(uniqueDt$name_charge_index)
  nChar_sequences <- nchar(uniqueDt$Name)
  maxSeqLength <- max(nChar_sequences)
  
  #Output array (uniqSeq, maxSeqLen, aa)
  m <- array(data = 0L, 
             dim = c(nSequences, maxSeqLength, length(aa)), 
             dimnames = list(peptide = paste0(uniqueDt$Name, "_", uniqueDt$Charge), 
                             sequenceIndex = c(1:maxSeqLength), 
                             aa = aa))
  
  #Generate lists of equal length containing all the info
  ##List of sequence indexs
  sequence <- strsplit(uniqueDt$Name, split = "")
  aa_index <- lapply(X = sequence, FUN = match, table = aa)
  sequence_position_index <- lapply(X = aa_index, FUN = function(X){c(1L:length(X))})
  
  #Index locations within array
  aa_mat_index <- mapply(FUN = cbind, uniqueDt$name_charge_index, sequence_position_index, aa_index, SIMPLIFY = FALSE)
  aa_mat_index <- do.call(what = rbind, args = aa_mat_index)
  
  m[aa_mat_index] <- 1L
  m
}

#' Integer Sequence Tensor Encoding of Sequence Data
#'
#' @param uniqueDt a data.table containing the sequence "Name", the charge
#'   "Charge", and an index specifying unique sequence-charge combinations
#'   "name_charge_index" (one sequence can have multiple charge states)
#' @param aa a character vector of amino acid codes
#'
#' @return
#' @export
#'

sequence2tensor_intSequence <- function(uniqueDt, aa){
  #number of sequences to process
  nSequences <- max(uniqueDt$name_charge_index)
  maxSeqLength <- max(nchar(uniqueDt$Name))
  
  #Output array (uniqSeq, maxSeqLen, aa)
  m <- array(data = 0L, 
             dim = c(nSequences, maxSeqLength, 1), 
             dimnames = list(peptide = paste0(uniqueDt$Name, "_", uniqueDt$Charge), 
                             sequenceIndex = c(1:maxSeqLength)))
  
  #Generate lists of equal length containing all the info
  ##List of sequence indexs
  sequence <- strsplit(uniqueDt$Name, split = "")
  aa_index <- lapply(X = sequence, FUN = match, table = aa)
  
  for(i in c(1:length(aa_index))){
    m[i,c(1:length(aa_index[[i]])),1] <- aa_index[[i]]
  }
  
  m
}























