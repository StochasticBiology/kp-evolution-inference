#' Preprocess experiment data
#' 
#' This function dichotomizes the data by any resistance gene group present
#' 
#' @param threshold number from 0 to 1 denoting the lower bound threshold for 
#'  gene presence in order to consider a feature.
#' @examples
#' prepocess.kleborate()
#' @export
preprocess.kleborate <- function(threshold = 0.01) {
  if (!file.exists("raw/Klebsiella pneumoniae__kleborate.csv")) {
    stop("Missing download: Klebsiella pneumoniae__kleborate.csv")
  }
  
  kleborate.df <- read.csv("raw/Klebsiella pneumoniae__kleborate.csv")
  # All resistance columns
  selection <- colnames(kleborate.df)[endsWith(colnames(kleborate.df),
                                                   "_acquired") | 
                                      endsWith(colnames(kleborate.df),
                                                 "_mutations") |
                                      colnames(kleborate.df) == "Bla_chr"]
  
  resistance.df <- kleborate.df[, selection]
  
  # Dichotomize by presense or absense
  resistance.df <- as.data.frame(apply(resistance.df,
                                       c(1,2),
                                       \(x) ifelse(x=="-", 0, 1)))
  
  threshold <- lapply(colnames(resistance.df), 
                      \(x) sum(resistance.df[,c(x)]) / length(resistance.df[,1]) >= threshold)
  
  resistance.df <- resistance.df[,colnames(resistance.df)[unlist(threshold)]]

  resistance.df <- cbind(id=kleborate.df$Genome.ID,
                         resistance.df)
  
  dir.create("clean/", showWarnings=FALSE)
  write.csv(resistance.df,
            "clean/kleborate-dichotomized.csv",
            row.names=FALSE)
}
