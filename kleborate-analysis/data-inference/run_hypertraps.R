#' Run the model
#' 
#' This function runs HyperTraPS.
#' 
#' @param country which country to run
#' @param seed the random seed
#' @param length steps in 10 ** length
#' @param hypertraps.path the location of the cpp source file
#' @examples
#' run.hypertraps() 
#' @export
run.hypertraps <- function(country = "Germany", 
                           my.seed = 1, 
                           run.length = 2,
                           walkers = 200) {
  
  if (!(require("Rcpp") &&
        require("ape"))) stop("Missing required dependencies Rcpp and ape")
  
  tree.path <- paste0("clean/",country,".nwk")
  if (!file.exists(tree.path)) {
    stop(paste("No newick-tree for", country, "in clean directory"))
  }
  
  if (!file.exists("clean/kleborate-dichotomized.csv")) {
      stop("Run preprocess_kleborate.R first!")
  }

  resistance.df <- read.csv("clean/kleborate-dichotomized.csv")
  featurenames <- setdiff(colnames(resistance.df), "id")
  
  ctree <- hypertrapsct::curate.tree(tree.path, 
                       "clean/kleborate-dichotomized.csv")
  
  model.fit <- hypertrapsct::HyperTraPS(ctree$dests, # set to one thousand walkers
                       initialstates=ctree$srcs,
                       length = run.length,
                       walkers = walkers,
                       penalty = 1,
                       seed = my.seed,
                       featurenames=featurenames)
  
  return (model.fit)
}
