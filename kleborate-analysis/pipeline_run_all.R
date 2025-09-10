# Before running this all files and newicks must be prepared. I.e run 
# "Rscript pipeline_preprocessing.R" and "bash make_trees.sh" in exec-folder.

library(hypertrapsct)

args <- commandArgs(trailing=TRUE)

if (length(args) != 1)
  stop("USAGE: Rscript pipeline_run_all.R [cluster_id]")

path <- "data-inference"  # replace with your folder path
files <- list.files(path, pattern = "\\.R$", full.names = TRUE)
sapply(files, source)

length = 6
walkers = 1000
counter = 0
cluster_id = as.numeric(args[[1]])
n_clusters = 8

# Create parameters and metadata for all runs
countries <- gsub("\\.nwk", "", list.files("clean", pattern=".*\\.nwk"))
seeds <- c(1,2,3)
load("misc-data/uninformativecountries.Rdata")
runs <- expand.grid(uninformative.countries$country, seeds)
colnames(runs) <- c("country", "seed")
load("data/tree_metrics.Rdata")
runs <- dplyr::left_join(runs, trees.df, keep=FALSE)

# Run models from smallest to largest transition set, dividing them evenly among
# the cores.
dir.create("models/length6-0L")
for (i in order(runs$n_transitions, runs$seed)) {
  # Check if this thread should run model
  counter <- counter + 1
  if ((cluster_id + counter) %% n_clusters != 0) {
    next
  }

  country <- runs[i, "country"]
  seed <- runs[i, "seed"]

  # Do some logging
  sink(file = paste0("logfile-",cluster_id,".txt"),
       append = TRUE)
  print(paste(Sys.time(), country, seed, runs[i, "n_transitions"]))
  sink()
  
  # Run model
  model <- run.hypertraps(country, seed, length, walkers)
  
  # Rename model to country.seed in global env and save
  new_name = paste(country, seed, sep=".")
  list2env(setNames(list(model), new_name), .GlobalEnv)
  save(list=new_name,
       file=paste0("models/length6-0L/",
                   country,
                   "-length-",length,
                   "-seed-",seed,
                   "-walkers-",walkers,
                   ".Rdata"))
}

# Log end time
sink(file = paste0("logfile-",cluster_id,".txt"),
     append = TRUE)
print(paste(Sys.time(), "DONE", runs[i, "n_transitions"]))
sink()
