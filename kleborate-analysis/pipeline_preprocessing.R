# This script downloads and preprocesses all files necessary for the analysis

path <- "data-inference"  # replace with your folder path
files <- list.files(path, pattern = "\\.R$", full.names = TRUE)
sapply(files, source)

download.pathogen.watch()
preprocess.kleborate(threshold=0)
preprocess.metadata()

prepare.all.lincoding.files()

setwd("trees/")
system("chmod +x make_trees.sh; ./make_trees.sh")
