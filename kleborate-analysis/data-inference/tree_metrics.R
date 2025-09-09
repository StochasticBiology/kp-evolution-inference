#' Extract number of transitions and feature statistics for a curated tree
#' 
#' @param ctree A list. Output from hypertraps-ct::curate.tree
#' @return A list of n_transition and mean acquisition of all features
#' @export
tree.metrics <- function(c.tree) {
  features <- colnames(c.tree$data)[-1]
  counts <- vector("list", length(features))
  
  for (i in 1:length(c.tree$transitions$from)) {
    from <- c.tree$transitions$from[[i]]
    to <- c.tree$transitions$to[[i]]
    
    n_present_ancestor = stringr::str_count(from, "1")
    for (j in 1:stringr::str_length(from)) {
      if (substring(from,j,j) != substring(to, j,j)) {
        counts[[j]] <- c(counts[[j]], n_present_ancestor)
      }
    }
  }
  
  mean.na <- function(n) {
    if (is.numeric(n)) {
      return (mean(n))
    } else {
      return (NA)
    }
  }
  
  return (list(transitions = length(c.tree$transitions$from), 
               mean.acquisition = sapply(counts, mean.na)))
}

#' Calculate tree metrics for all newicks in clean and return a data frame
#' 
#' @return A data frame containing tree metrics for all newick trees
#' @export
get.tree.metrics.df <- function() {
  source("../hypertraps-ct/hypertraps.R")
  kleborate.df <- read.csv("clean/kleborate-dichotomized.csv")
  metrics.df <- data.frame(matrix(ncol = length(colnames(kleborate.df)) + 1, 
                                  nrow = 0))
  
  for (file in list.files("clean", pattern = "*\\.nwk")) {
    print(paste("parsing", file))
    c.tree <- curate.tree(paste0("clean/",file), 
                          "clean/kleborate-dichotomized.csv")
    metrics <- tree.metrics(c.tree)
    data.row <- c(gsub("(*)\\.nwk","\\1",file), 
                  as.numeric(metrics$transitions), 
                  as.numeric(metrics$mean.acquisition))
    
    metrics.df <- rbind(metrics.df, data.row)
  }
  
  colnames(metrics.df) <- c("country", 
                            "n_transitions", 
                            colnames(kleborate.df)[-c(1)])
  
  metrics.df$country <- factor(metrics.df$country)
  metrics.df[is.na(metrics.df)] <- Inf
  
  return (metrics.df)
}