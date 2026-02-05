library(phytools)
library(hypertrapsct)
library(ggpubr)
library(parallel)

country = "Romania"

tree = read.tree(paste0("clean/", country, ".nwk", collapse=""))
resistance.df <- read.csv("clean/kleborate-dichotomized.csv")
df = resistance.df[resistance.df$id %in% tree$tip.label,]
ct = curate.tree(tree, df)

# wrapper function for HyperTraPS analysis, setting random seed and which data subset we're using
parallel.fn = function(seed) {
  if(seed < 4) {
    HyperTraPS(ct$dests,
               length = 4.5, kernel = 4, penalty=1, walkers = 500)
  } else {
    HyperTraPS(ct$dests,
               initialstates = ct$srcs,
               length = 4.5, kernel = 4, penalty=1, walkers = 500)
  }
}

# run these experiments in parallel. 
parallelised.runs <- mcmapply(parallel.fn, seed=1:6, 
                              SIMPLIFY = FALSE,
                              mc.cores = min(detectCores()-1, 6))

for(i in 1:6) {
  parallelised.runs[[i]]$featurenames = colnames(df)[2:ncol(df)]
}

sf = 2
fname = paste0(country, "-compare-both.png", collapse="")
png(fname, width=900*sf, height=500*sf, res=72*sf)
ggarrange(plotHypercube.curated.tree(ct, hjust = 1, scale.fn = NULL, font.size = 2.5) + theme(
  plot.margin = margin(t = 10, r = 10, b = 100, l = 10)) + coord_cartesian(clip = "off"),
  plotHypercube.bubbles.compare(parallelised.runs, sqrt.trans = TRUE),
  widths=c(1,2), labels=c("A", "B"))
dev.off()

fname = paste0(country, "-data.Rdata", collapse="")
save(parallelised.runs, file=fname)
