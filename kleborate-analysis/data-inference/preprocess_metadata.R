#' Preprocess the metadata
#' 
#' This cleans inconsistencies in the metadata. Primarily different variations
#' for names of different countries.
#' 
#' @examples
#' prepocess.metadata()
#' @export
preprocess.metadata <- function() {
  if (!file.exists("raw/Klebsiella pneumoniae__metadata.csv")) {
    stop("Missing download: raw/Klebsiella pneumoniae__metadata.csv")
  }
  
  metadata.df <- read.csv("raw/Klebsiella pneumoniae__metadata.csv")
  
  for (i in seq_along(metadata.df$Country)) {
    country <- metadata.df$Country[[i]]
    renamed <- switch(country,
                      "UK" = "United Kingdom",
                      "United Kingdom (Scotland)" = "United Kingdom",
                      "United Kingdom (England/Wales/N. Ireland)" = "United Kingdom",
                      "Viet Nam" = "Vietnam",
                      "UAE" = "United Arab Emirates",
                      "KSA" = "Saudi Arabia",
                      "West Bank" = "Palestine",
                      "The Gambia" = "Gambia",
                      "Other (International Space Station)" = "",
                      "missing" = "",
                      country
    ) 
    metadata.df$Country[[i]] <- gsub(" ", "_", renamed)
  }
  
  dir.create("clean/", showWarnings=FALSE)
  write.csv(metadata.df,
            "clean/metadata.csv",
            row.names=FALSE)
}
