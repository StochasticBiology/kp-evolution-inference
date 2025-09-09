#' Short featurenames
#' @param feature.names A character vector
#' @return A character vector
#' @export
shorten.featurenames <- function(feature.names) {
  short.names <- gsub("(.*)_acquired","\\1_a",feature.names)
  short.names <- gsub("(.*)_mutations","\\1_m",short.names)
  short.names <- gsub("(Bla_)?(.*_[a|m])","\\2",short.names)
  short.names
}

#' Extract features from model
#' @param model A list object returned by HyperTraPS
#' @return A character vector of featurenames
#' @export
extract.features <- function(model) {
  features <- model$bubbles[model$bubbles$Time == 0,]
  features <- features[order(features$OriginalIndex),]
  features$Name
}