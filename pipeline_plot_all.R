devtools::load_all()
source(path.hypertraps)
library(ggpubr)

read.models <- function(path, pattern="") {
  models <- list()

  #path <- "/media/olav/VERBATIM/length4/"
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

#save(models, file="all_models.Rdata")
setwd("..")
name <- load("all_models.Rdata")
country.list <- get(name)
countries <- names(country.list)
#models.bayesian <- read.models("~/hyperevol/remote/length6-1000/", "walkers-1000")
#models.root <- read.models("~/hyperevol/remote/length6-0L/")
#models <- models.bayesian[setdiff(names(models.bayesian), names(models.root))]
#models <- c(models, models.root)
#models.sa <- read.models("models/length4-new/", "-sa")
#just.sa <- setdiff(names(models.sa), names(models.bayesian))
#models.sa <- models.sa[just.sa]

#save(models.bayesian, models.sa, file = "all_102_models.Rdata")

#countries <- c(names(models.sa), names(models.bayesian))
#country.list <- c(models.sa, models.bayesian)

plot.path <- "plots/june/"
dir.create(plot.path)

name <- load("data/tree_metrics.Rdata")
trees.df <- get(name)
transitions <- as.numeric(trees.df$n_transitions)
names(transitions) <- trees.df$country
transitions <- transitions[names(transitions) %in% countries]


cluster.cubes(country.list, transitions)
ggsave(paste0(plot.path, "mds-embedding-countries-",length(country.list),".png"), width = 15, height = 12)

#cluster.cubes(country.list[transitions > 300], transitions[transitions > 300])
#ggsave(paste0("plots/length6/mds-embedding-countries-",length(country.list),"-length-6.png"), width = 15, height = 12)

plt <- plot.bubbles.by.feature(country.list, transitions)
plt
ggsave(paste0(plot.path, "features-sort.png"), width = 30, height = 15)

#plt.threshold.gt.1000 <- plot.bubbles.by.feature(country.list[transitions > 1000])
#plt.threshold.gt.1000 + ggtitle("Transition sets larger than 1000")

#plt.threshold.gt.300 <- plot.bubbles.by.feature(country.list[transitions > 300])
#plt.threshold.gt.300 + ggtitle("Transition sets larger than 300")

#plt.threshold.lt.100 <- plot.bubbles.by.feature(country.list[transitions < 100])
#plt.threshold.lt.100 + ggtitle("Transition sets smaller than 100")
#ggsave("plots/length6/features-threshold.png", width = 30, height = 5)

#plt.threshold.lt.50 <- plot.bubbles.by.feature(country.list[transitions < 50])
#plt.threshold.lt.50 + ggtitle("Transition sets smaller than 50")

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


#plotHypercube.bubbles(country.list[["Norway"]]$seed.1)
#ggplot(country.list[["Norway"]]$seed.1$bubbles, aes(x=Time, y=Name, size=Probability)) + geom_point() +
#  theme_light() 
#names(country.list[["Italy"]]$seed.1)
#names(country.list[["Norway"]]$seed.1)
#names(models.bayesian[[1]]$seed.1)

#models.old <- models.bayesian[names(models.root)]

#for (name in names(models.root)) {
#  old = models.old[name]
#  new = models.root[name]
  
#  bub.old <- old[[1]]$seed.1$bubbles
#  bub.new <- new[[1]]$seed.1$bubbles
#  bub.old$enforce.root = "false"
#  bub.new$enforce.root = "true"
#  bub = rbind(bub.old, bub.new)
#  ggplot(bub, aes(x=Time, y=Name, size=Probability, color=enforce.root)) + geom_point(alpha=0.5) + ggtitle(name)
#  ggsave(paste0(plot.path, "compare-",name, ".png"))
#}

