#' Download pathogen watch data
#' 
#' This function downloads all datasets used in the study
#' 
#' @param overwrite Logical
#' @examples
#' download.pathogen.watch()
#' @export
download.pathogen.watch <- function(overwrite=FALSE) {
  if (.Platform$OS.type == "windows") {
  	download.method = "wininet"
  
  } else if (.Platform$OS.type == "unix") {
  	download.method = "wget"
  
  } else {
  	stop("Unknown Operating System")
  }
  
  pw.endpoint <- "https://pathogenwatch-public.ams3.cdn.digitaloceanspaces.com/"
  
  files <- c("Klebsiella pneumoniae__kleborate.csv",
  	         "Klebsiella pneumoniae__cgmlst.csv",
  	         "Klebsiella pneumoniae__metadata.csv")
  
  for (filename in files) {
  	if (!file.exists(paste0("raw/",filename)) || overwrite) {
  		dir.create("raw/", showWarnings=FALSE)
  		utils::download.file(url=paste0(pw.endpoint,filename,".gz"),
  			      destfile=paste0("raw/",filename,".gz"),
  			      method=download.method)
  		R.utils::gunzip(paste0("raw/",filename,".gz"), overwrite = overwrite)
  		
  	} else {
  	  print(paste("File", filename, "is already downloaded in folder raw."))
  	}
  }
}

