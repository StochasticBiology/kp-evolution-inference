library(phytools)
library(hypertrapsct)
library(ggpubr)
library(parallel)

run.inference = TRUE
ct = list()

for(country in c("Romania", "Senegal")) {
  
  tree = read.tree(paste0("clean/", country, ".nwk", collapse=""))
  resistance.df <- read.csv("clean/kleborate-dichotomized.csv")
  df = resistance.df[resistance.df$id %in% tree$tip.label,]
  ct[[country]] = curate.tree(tree, df)
  
  # wrapper function for HyperTraPS analysis, setting random seed and which data subset we're using
  parallel.fn = function(seed) {
    if(seed < 4) {
      HyperTraPS(ct[[country]]$dests,
                 length = 4.5, kernel = 4, penalty=1, walkers = 500)
    } else {
      HyperTraPS(ct[[country]]$dests,
                 initialstates = ct$srcs,
                 length = 4.5, kernel = 4, penalty=1, walkers = 500)
    }
  }
  
  if(run.inference == TRUE) {
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
    ggarrange(plotHypercube.curated.tree(ct[[country]], hjust = 1, scale.fn = NULL, font.size = 2.5) + theme(
      plot.margin = margin(t = 10, r = 10, b = 100, l = 10)) + coord_cartesian(clip = "off"),
      plotHypercube.bubbles.compare(parallelised.runs, sqrt.trans = TRUE),
      widths=c(1,2), labels=c("A", "B"))
    dev.off()
    
    sf = 4
    fname = paste0(country, "-compare-both-2.png", collapse="")
    png(fname, width=400*sf, height=200*sf, res=72*sf)
    ggarrange(
      plotHypercube.curated.tree(ct[[country]], hjust = 1, scale.fn = NULL, font.size = 0) + 
        #theme( plot.margin = margin(t = 10, r = 10, b = 100, l = 10)) + 
        coord_cartesian(clip = "off"),
      plotHypercube.bubbles.compare(parallelised.runs[c(1,4)], 
                                    sqrt.trans = TRUE, 
                                    expt.names = c("No phylo", "Phylo") ) +
        theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank()),
      widths=c(1,2))
    dev.off()
    
    fname = paste0(country, "-data.Rdata", collapse="")
    save(parallelised.runs, file=fname)
  }
  
}

#######

