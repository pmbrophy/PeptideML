#' Import NIST MSP File
#'
#' @description Import a NIST formatted .msp file and convert columns to numeric
#'   data types
#'
#' @param mspPath file path to NIST .msp file
#'
#' @return a data.table
#' @export
#' 

import_NIST_MSP <- function(mspPath){
  library(mspReader)
  library(data.table)
  
  #Import
  mspData <- readMsp(path = mspPath, commentType = "NIST")
  
  #Merge data and info
  mspData <- merge(x = mspData$data, 
                   y = mspData$info[, .(index, Name, Fullname, Charge, Mods, Mz_exact, Parent, Pep)], 
                   by = "index")
  
  setorder(x = mspData, Fullname, Charge, Mods)
  
  #Convert to Numeric
  mspData[, mz := as.numeric(mz)]
  mspData[, Mz_exact := as.numeric(Mz_exact)]
  mspData[, Charge := as.numeric(Charge)]
  mspData[, intensity := as.numeric(intensity)]
  mspData[, Parent := as.numeric(Parent)]
  
  mspData
}

#' Cleanup Sequence/Name
#'
#' @param mspData a data.table
#'
#' @return same data.table with Name containing only a aa sequence
#' @export
#'

cleanupSequence <- function(mspData){
  cleanName <- strsplit(x = mspData$Name, "/")
  cleanName <- sapply(X = cleanName, "[[", 1)
  
  mspData[, Name := cleanName]
  
  remove(cleanName)
  
  mspData
}

#' Calculate sequence Length
#'
#' @param mspData a data.table returned by cleanupSequence()
#'
#' @return the same data.table also containing a column "seqLen" containing the length of the peptide sequence
#' @export
#'

calcSequenceLength <- function(mspData){
  mspData[, seqLen := nchar(Name)]
  
  mspData
}

#' Remove Modifications
#'
#' @param mspData a data.table returned by import_NIST_MSP
#'
#' @return a data.table with unmodified peptides 
#' @export
#'

removeMods <- function(mspData){
  mspData[, Mods := as.numeric(Mods)]
  mspData <- mspData[Mods == "0",]
  
  mspData
}

#' Process Annotations
#'
#' @description Process ion annotations and remove ions without ion-type
#'   assignment
#'
#' @param mspData a data.table
#'
#' @return return the same data.table only containing ions with annotations.
#'   Annotations are contained in two columns b/c some ions have multiple
#'   possible annotations
#' @export
#'

processAnnotations <- function(mspData){
  #1 Remove ions without annotation from data ("?")
  mspData <- mspData[annotation != "\"?\""]
  
  #2 Remove leading and closing quotes from annotation
  mspData[, annotation_nchar := nchar(mspData$annotation)]
  mspData[, annotation := substr(x = annotation, start = 2, stop = (annotation_nchar-1))]
  mspData[, annotation_nchar := NULL]
  
  #3 - String split by `,`
  #no more than 2 annotations exist
  splits <- strsplit(x = mspData$annotation, split = ",")
  mspData[,annotation := NULL]
  splits <- do.call(rbind, splits)
  
  #4 Split by `\`
  splits <- apply(FUN = strsplit, X = splits, split = "/", MARGIN = 2)
  
  #5 Save Data
  mspData[, annotation_1 := sapply(splits[[1]], "[[", 1)]
  mspData[, annotation_1_err := as.numeric(sapply(splits[[1]], "[[", 2))]
  mspData[, annotation_2 := sapply(splits[[2]], "[[", 1)]
  mspData[, annotation_2_err := as.numeric(sapply(splits[[2]], "[[", 2))]
  remove(splits)
  
  mspData
}


#' Remove ion types other than y and b
#'
#' @description Filter out all ions except pure y/b ion types. y/b ions that
#'   undergo additional neutral loss are removed. All ions start with either 1
#'   or 2 annotations. For ions with 2 annotations, only one of the annotations
#'   must be a pure y/b ion to be allowed through.
#'
#' @param mspData the data.table modified by processAnnotations()
#'
#' @return return the same data.table with only pure-y/b ions
#' @export
#' 

select_y_b_ions <- function(mspData){
  #1) Unique annotations
  targetIonTypes <- unique(mspData$annotation_1)
  
  #2) Remove neutral loss
  targetIonTypes <- targetIonTypes[!grepl(x = targetIonTypes, pattern = "[-]")]
  targetIonTypes <- targetIonTypes[!grepl(x = targetIonTypes, pattern = "[+]")]
  
  #3) Select only y/b
  targetIonTypes1 <- targetIonTypes[grep(x = targetIonTypes, pattern = "[yb]")]
  
  #Repeat steps 1-3 on annotation_2
  targetIonTypes <- unique(mspData$annotation_2) #1
  targetIonTypes <- targetIonTypes[!grepl(x = targetIonTypes, pattern = "[-]")] #2
  targetIonTypes <- targetIonTypes[!grepl(x = targetIonTypes, pattern = "[+]")] 
  targetIonTypes2 <- targetIonTypes[grep(x = targetIonTypes, pattern = "[yb]")] #3
  
  #Get unique ion types from annotation_1 and annotation_2
  targetIonTypes <- unique(targetIonTypes1, targetIonTypes2) 
  
  #4) Filter both annotations to those containing targetIonTypes
  mspData <- mspData[annotation_1 %in% targetIonTypes | annotation_2 %in% targetIonTypes,]
  mspData[, goodAnnotation := ifelse(test = annotation_1 %in% targetIonTypes, 
                                  yes = annotation_1, 
                                  no = annotation_2)]
  
  #5) Remove original annotation data
  remove(targetIonTypes, targetIonTypes1, targetIonTypes2)
  mspData[, annotation_1 := NULL]
  mspData[, annotation_1_err := NULL]
  mspData[, annotation_2 := NULL] 
  mspData[, annotation_2_err := NULL]
  
  mspData
}

#' Parse y/b annotations
#'
#' @description y/b ion types have the form "b10^2". Produce a columns for: 1)
#'   the ion type "goodAnnotation": y or b. 2) the fragment ion charge state
#'   "fragChargeState". 3) fragmentation location "fragIndex".
#'
#' @param mspData
#'
#' @return
#' @export
#'
#' @examples
parse_y_b_annotations <- function(mspData){
  #1) Get y/b ionType from goodAnnotation(s)
  mspData[, ionType := substr(goodAnnotation, start = 1, stop = 1)]
  
  #2) Get fragment charge state "b3^2" --> "b3", 2
  splits <- strsplit(x = mspData$goodAnnotation, split = "[[:punct:]]") 
  
  #fill in singly charged...
  splits <- lapply(X = splits, 
                   FUN = function(X){
                     if(length(X) == 1){
                       c(X, "1")
                     }else{
                       X
                     }
                   })
  
  mspData[, goodAnnotation := sapply(splits, "[[", 1)]
  mspData[, fragChargeState := as.numeric(sapply(splits, "[[", 2))]
  remove(splits)
  
  #3) Fragment location from b#/y#
  mspData[, fragIndex := as.numeric(substr(goodAnnotation, start = 2, stop = nchar(goodAnnotation)))]
}

#' Calculate fragment ion aa sequences
#'
#' @description calculate the fragment ion amino acid sequence and length of the
#'   fragment sequence. Only works for "y" and "b" ions without neutral loss or
#'   modifications
#'
#' @param mspData a data.table
#'
#' @return a data.table containing fragment ion sequence and length of sequence
#'
#' @export
#' 

calculate_fragIonSequences <- function(mspData){
  
  sequences <- mapply(FUN = calculate_fragIonSequence, 
                      ionType = mspData$ionType, 
                      fragmentLocation = mspData$fragIndex, 
                      precursorSequence = mspData$Name, 
                      sequenceLength = mspData$seqLen)
  
  mspData[,fragSequence := sequences]
  remove(sequences)
  mspData[,fragLength := nchar(fragSequence)]
  
  mspData
}

#' Calculate fragment ion aa sequence 
#'
#' @param ionTypes the type of fragment: "y" or "b"
#' @param fragmentLocations the location of fragmentation
#' @param precursorSequences the amino acid sequence of the precursor
#' @param sequenceLengths the length of the amino acid sequence
#'
#' @return a character vector of the fragment ion sequence
#' @export
#'

calculate_fragIonSequence <- function(ionType,  fragmentLocation, precursorSequence, sequenceLength){
  if(ionType == "b"){
    substr(x = precursorSequence, start = 1, stop = fragmentLocation)
    
  }else if(ionType == "y"){
    substr(x = precursorSequence, start = sequenceLength - fragmentLocation + 1, stop = sequenceLength)
    
  }else{
    stop("ion type not defined")
  }
}


#' Calculate the exact mass of a fragment ion from an amino acid sequence
#'
#' @param mspData the data.table
#'
#' @return return data.table with extra column "fragExactMass"
#' @export
#'

calculateSequence_ExactMass <- function(mspData){
  #get element object
  elems <- lapply(FUN = OrgMassSpecR::ConvertPeptide, X = mspData$fragSequence, IAA = FALSE)
  
  #exact mass
  fragExactMass <- sapply(FUN = OrgMassSpecR::MonoisotopicMass, X = elems, charge = 0)
  remove(elems)
  
  mspData[, fragExactMass := fragExactMass]
  remove(fragExactMass)
  
  mspData
}

#' Calculate mz of each fragment
#'
#' @param mspData a data.table
#'
#' @return the data.table plus the column "theorFragMz"
#' @export
#'

calculateSequence_mzs <- function(mspData){
  mz <- mapply(FUN = calculateSequence_mz,
               ionType = mspData$ionType, 
               exactMass = mspData$fragExactMass, 
               chargeState = mspData$fragChargeState)
  
  mspData[, theorFragMz := mz]
  
  mspData
}

#' Calculate mz of a single ion 
#'
#' @param ionType the type of the ion "y" or "b"
#' @param exactMass the calculated exact mass 
#' @param chargeState the charge state
#'
#' @return the m/z value
#' @export
#'

calculateSequence_mz <- function(ionType, exactMass, chargeState){
  if(ionType == "y"){
    mz <- (exactMass + (chargeState * 1.007825))/chargeState

  }else if(ionType == "b"){
    mz <- (exactMass - 18.01517 + (chargeState * 1.007825))/chargeState
    
  }else{
    stop("ion type not implemented")
  }
  
  mz
}
