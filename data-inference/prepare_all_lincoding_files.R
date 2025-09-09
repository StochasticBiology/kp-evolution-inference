#' Create a LIN-file for a specific country
#' 
#' This function prepares the input data for LINcoding.py
#' 
#' @param country A country present in the dataset
#' @export
create.country.lin.file <- function(country="Germany") {
  if (!file.exists("raw/Klebsiella pneumoniae__cgmlst.csv")) {
    download.pathogen.watch()
  }
  if (!file.exists("clean/metadata.csv")) {
    preprocess.metadata()
  }
  
  cgmlst.df <- utils::read.csv("raw/Klebsiella pneumoniae__cgmlst.csv")
  meta.df <- utils::read.csv("clean/metadata.csv")
  
  samples <- meta.df$id[which(meta.df$Country == country)]
  cgmlst.subset <- cgmlst.df[which(cgmlst.df$Genome.ID %in% samples),]
  
  lin.codes <- stringr::str_split(cgmlst.subset$code, "_")
  lin.codes.numeric <- sapply(lin.codes, as.numeric)
  
  cbind(Genome.ID=cgmlst.subset$Genome.ID, data.frame(t(lin.codes.numeric)))
}

#' Prepare all LINcoding.py input files
#' 
#' This function creates a .tsv file for all countries listed in the metadata
#' file. 
#' 
#' @param overwrite Logical.
#' @export
prepare.all.lincoding.files <- function(overwrite=FALSE) {
  df <- read.csv("clean/metadata.csv")
  countries <- unique(df$Country)
  for (country in countries) {
    fn_country <- chartr(" ", "_", country)
    if (!file.exists(paste0("clean/", fn_country, ".tsv")) || overwrite) {
      print(paste("Doing", country))
      
      df <- create.country.lin.file(country)
      if (length(df[,1]) == 1) next
      
      write.table(df, paste0("clean/",country,".tsv"), 
                row.names=FALSE, 
                quote=FALSE,
                sep='\t', 
                na="")
    } else {
      print(paste("File already done and overwrite=FALSE", country))
    }
  }
}
