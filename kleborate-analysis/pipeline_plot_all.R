library(hypertrapsct)
library(ggpubr)

path <- "data-inference"  # replace with your folder path
files <- list.files(path, pattern = "\\.R$", full.names = TRUE)
sapply(files, source)

read.models <- function(path, pattern="*.Rdata") {
  models <- list()

  for (file in list.files(path, pattern = pattern)) {
    object.name <- load(file=paste0(path,file))
    models <- c(models, object.name)
  }
  
  countries <- unique(gsub("\\.[1-3]","",models))
  country.list <- vector("list", length(countries))
  names(country.list) <- countries
  
  unfinished <- character()
  for (country in countries) {
    if (all(sapply(c(".1", ".2", ".3"), \(x) exists(paste0(country, x))))) {
      country.list[[country]]$seed.1 <- get(paste0(country, ".1"))
      country.list[[country]]$seed.2 <- get(paste0(country, ".2"))
      country.list[[country]]$seed.3 <- get(paste0(country, ".3"))
    } else {
      print(paste("Not all seeds are done for", country))
      unfinished <- c(country, unfinished)
    }
  }
  return (country.list)
  countries <- setdiff(countries,unfinished)
  country.list <- country.list[countries]
  names(country.list) <- countries
  return (country.list)

}

models <- read.models("models/length6-0L/")
save(models, file="all_models.Rdata")

name <- load("all_models.Rdata")
country.list <- get(name)
countries <- names(country.list)

plot.path <- "plots/"
dir.create(plot.path)

name <- load("misc-data/tree_metrics.Rdata")
trees.df <- get(name)
transitions <- as.numeric(trees.df$n_transitions)
names(transitions) <- trees.df$country
transitions <- transitions[names(transitions) %in% countries]

cluster.cubes(country.list, transitions)
ggsave(paste0(plot.path, "mds-embedding-countries-",length(country.list),".png"), width = 15, height = 12)

plt <- plot.bubbles.by.feature(country.list, transitions)
plt
ggsave(paste0(plot.path, "features-sort.png"), width = 30, height = 15)

featurenames <- extract.features(country.list[[1]]$seed.1)
s.names <- shorten.featurenames(featurenames)

for (country in countries) {
  if (file.exists(paste0(plot.path, country, ".png"))) next
  models <- country.list[[country]]
  ggarrange(plotHypercube.lik.trace(models$seed.1),
            plotHypercube.lik.trace(models$seed.2),
            plotHypercube.lik.trace(models$seed.3),
            plotHypercube.bubbles(models$seed.1) + theme(legend.position = "none"),
            plotHypercube.bubbles(models$seed.2) + theme(legend.position = "none"),
            plotHypercube.bubbles(models$seed.3) + theme(legend.position = "none"),
            plotHypercube.influences(models$seed.1) + theme(legend.position = "none"),
            plotHypercube.influences(models$seed.2) + theme(legend.position = "none"),
            plotHypercube.influences(models$seed.3) + theme(legend.position = "none")
            )
  ggsave(paste0(plot.path, country, ".png"), width = 15, height = 12)
}

