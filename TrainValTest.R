#' Split integer sequence array and fragment array
#'
#' @param sequenceArray integer sequence data to be split
#' @param fragArray fragment array to be split
#' @param frac_test fraction of total data to be reserved for testing
#' @param frac_train fraction of data after removing test samples to be used for training. 
#'
#' @return a list containing train, validate, and test data
#' @export
#'

split_trainValTest <- function(sequenceArray, fragArray, frac_test = 0.1, frac_train = 0.7){
  #Index 
  sequenceArray_n <- nrow(sequenceArray)
  fragArray_n <- nrow(fragArray)
  array_i <- c(1:sequenceArray_n)
  
  if(sequenceArray_n != fragArray_n){
    stop("sequence array and fragment array differ in size")
  }
  
  ###Test Locations
  testLocs <- sample(x = array_i, size = floor(sequenceArray_n * frac_test), replace = FALSE)
  
  #Test arrays
  test_x <- sequenceArray[testLocs,,]
  test_y <- fragArray[testLocs,,]
  test_x_dims <- dim(test_x)
  test_x <- keras::array_reshape(x = test_x, dim = c(test_x_dims[1], test_x_dims[2], test_x_dims[3]))
  
  #Remove test samples from original data structure
  sequenceArray <- sequenceArray[-testLocs,,]
  fragArray <- fragArray[-testLocs,,]
  
  sequenceArray_dims <- dim(sequenceArray)
  fragArray_dim <- dim(fragArray)
  sequenceArray <- keras::array_reshape(x = sequenceArray, dim = c(sequenceArray_dims[1], sequenceArray_dims[2], test_x_dims[3]))
  
  ###Train/Validate
  #Split remaining into train/validate
  array_i <- c(1:sequenceArray_dims[1])
  trainLocs <- sample(x = array_i, size = floor(sequenceArray_dims[1] * frac_train), replace = FALSE)
  
  #Train
  train_x <- sequenceArray[trainLocs,,]
  train_y <- fragArray[trainLocs,,]
  train_x_dims <- dim(train_x)
  train_x <- keras::array_reshape(x = train_x, dim = c(train_x_dims[1], train_x_dims[2], sequenceArray_dims[3]))
  
  #Validate
  val_x <- sequenceArray[-trainLocs,,]
  val_y <- fragArray[-trainLocs,,]
  val_x_dims <- dim(val_x)
  val_x <- keras::array_reshape(x = val_x, dim = c(val_x_dims[1], val_x_dims[2], sequenceArray_dims[3]))
  
  list(train = list(x = train_x, y = train_y),
       validate = list(x = val_x, y = val_y),
       test = list(x = test_x, y = test_y))
}