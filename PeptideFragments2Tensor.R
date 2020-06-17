#' Fragments to Tensor Encoding
#'
#' @param fragmentDt a data.table containing the "fragType_chargeCode",
#'   "ionType", "fragChargeState", "nFragments", "fragLength", and
#'   "intensity_norm"
#' @param uniqueDt unique charge/sequence data.table used to generate input encodings
#'
#' @return
#' @export
#'

fragments2tensor <- function(fragmentDt, uniqueDt){
  nSequences <- nrow(uniqueDt) #number of unique precursors
  maxSeqLength <- max(nchar(uniqueDt$Name)) #Max length of fragment
  
  #Array dimensions
  #seqLen <- maxSeqLength - 1
  seqLen <- maxSeqLength
  ionTypes <- max(fragmentDt$fragType_chargeCode)
  
  #Lookup table for fragment ion type and fragment charge to label array
  ionTypeLabel <- unique(fragmentDt[,list(fragType_chargeCode, ionType, fragChargeState)])
  setorder(ionTypeLabel, fragType_chargeCode)
  ionType_dimLabel <- vector(mode = "character", length = ionTypes)
  ionType_dimLabel[ionTypeLabel$fragType_chargeCode] <- paste0(ionTypeLabel$ionType, "_", ionTypeLabel$fragChargeState)
  
  m <- array(data = 0L, 
             dim = c(nSequences, seqLen, ionTypes), 
             dimnames = list(peptide = paste0(uniqueDt$Name, "_", uniqueDt$Charge), 
                             sequenceIndex = c(1:seqLen),
                             ionType = ionType_dimLabel))
  

  nFragments_perPepChargePair <- unique(fragmentDt[, list(name_charge_index, nFragments)])$nFragments
  
  sequenceIndex <- mapply(FUN = rep, x = c(1:nSequences), times = nFragments_perPepChargePair)
  sequenceIndex <- unlist(sequenceIndex)
  
  frag_mat_index <- cbind(sequenceIndex, fragmentDt$fragLength, fragmentDt$fragType_chargeCode)
  
  m[frag_mat_index] <- fragmentDt$intensity_norm
  m
}
