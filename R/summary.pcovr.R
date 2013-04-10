summary.pcovr <-
function(object, ...){
  cat("SELECTED SETTINGS BY USER\n")
  cat("Model selection procedure: ",object$modsel,"\n")
  cat("Rotation criterion: ",object$rot,"\n")
  y <- ifelse(object$prep=="cent","Centered data", "Standardized data")#hmm
  cat("Preprocessing: ",y,"\n\n")
  cat("DATA CHARACTERISTICS\n")
  cat("Number of observations: ",nrow(object$Te),"\n")
  cat("Number of predictors: ",ncol(object$Px),"\n")
  cat("Number of criteria: ",ncol(object$Py),"\n\n")
  cat("MODEL SELECTION\n")
  cat("Number of components: ",object$R,"\n")
  cat("Weighting parameter value: ",object$alpha,"\n\n")
  cat("OUTPUT\n")
  cat("Loading matrix:\n")
  print(t(object$Px))
  cat("\nRegression weight matrix:\n")
  print(t(object$Py))
  cat("\nExplained variance in X:",object$Rx2,"\n")
  cat("\nExplained variance in Y:",object$Ry2,"\n")
  if (object$modsel=="sim"){
    cat("\nCross-validation fit:\n")
    print(object$Qy2)
  }
  if (object$modsel=="seq"){
    cat("\nWeighted variance accounted for:\n")
    print(object$VAFsum)
  }
  if (object$modsel=="seqAcv"){
    cat("\nWeighted variance accounted for:\n")
    print(object$VAFsum)
    cat("\nCross-validation fit:\n")
    print(object$Qy2)
  }
  if (object$modsel=="seqRcv"){
    cat("\nCross-validation fit:\n")
    print(object$Qy2)
  }
}
