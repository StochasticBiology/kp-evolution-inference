library(ape)
library(phangorn)
library(parallel)
library(hypertrapsct)

# how much of a problem is reversibility?
x.list = my.ct = plot.ct = list()
for(expt in 1:6) {
  set.seed(1)
  
  # we will construct a dataset synthetically with known dynamics, then allow reversible transitions and try to infer the true behaviour
  L = 10
  if(expt == 1) {
    accumulation.rate = 1
    lambda.rev = rep(0, L)
  } else if(expt == 2) {
    accumulation.rate = 2
    lambda.rev = c(rep(0.1, L/2), rep(0, L/2))
  } else if(expt == 3) {
    accumulation.rate = 2
    lambda.rev = 0.5*(1:L)/L
  } else if(expt == 4) {
    accumulation.rate = 2
    lambda.rev = 2*(1:L)/L
  } else if(expt == 5) {
    accumulation.rate = 2
    lambda.rev = 0.5*runif(10)
  } else if(expt == 6) {
    accumulation.rate = 2
    lambda.rev = runif(10)
  }
  
  acc.rule = "hard"
  
  # parameterisation for tree construction
  tree.size = 128
  birth.rate = 1
  death.rate = 0.5
  
  # create random phylogeny with tree.size nodes from birth-death process parameterised as above
  my.tree = rphylo(tree.size, birth=birth.rate, death=death.rate)
  my.tree$node.label = as.character(1:my.tree$Nnode)
  tree.labels = c(my.tree$tip.label, my.tree$node.label)
  
  # generate state for all nodes traversing tree breadth-first
  # and setting the state of the child nodes according to
  # accumulation rate and branch length
  my.root = getRoot(my.tree)
  to.do = c(my.root)
  # initialise state list
  x = list()
  x[[my.root]] = rep(0,L)
  # while we still have vertices to simulate
  while(length(to.do) > 0) {
    # initialise a new to-do list for next iteration
    new.to.do = c()
    # loop through each node in current to-do list
    for(i in to.do) {
      this.outgoing.edges = which(my.tree$edge[,1] == i)
      # loop over this node's children
      for(j in this.outgoing.edges) {
        this.child = my.tree$edge[j,2]
        this.branch.length = my.tree$edge.length[j]
        # construct state for this child based on its parent
        x[[this.child]] = x[[i]]
        if(acc.rule == "hard") {
          # find leftmost zero in current state, and change with some probability
          # recall dynamics here are 00000 -> 10000 -> 11000 -> 11100 -> 11110 -> 11111
          ## (see first paragraph of section "Synthetic case studies")
          ref = which(x[[this.child]] == 0)[1]
          if(runif(1) < accumulation.rate*this.branch.length) { x[[this.child]][ref] = 1 } 
        } else if(acc.rule == "soft") {
          for(locus in 1:L) {
            if(runif(1) < accumulation.rate*((L+1-locus)/L)*this.branch.length) { x[[this.child]][locus] = 1 }
          }
        }
        # in the reversible case, allow the leftmost feature ("first feature" in the ms.:
        # second paragraph of "Synthetic case studies") to revert with some probability
        
        for(locus in 1:L) {
          if(runif(1) < lambda.rev[locus]*this.branch.length) { x[[this.child]][locus] = 0 }
        }
        # add this child to to state list, and to next iteration's to-do
        new.to.do = c(new.to.do, this.child)
      }
    }
    # update to-do list
    to.do = new.to.do
  }
  
  x.list[[expt]] = x
  
  # construct the dataframe of binary labels for tree tips
  df = data.frame()
  for(i in 1:length(my.tree$tip.label)) {
    df = rbind(df, data.frame(
      cbind(data.frame(id=my.tree$tip.label[i]), matrix(x[[i]], ncol=L))
    ))
  }
  
  # HyperTraPS pipeline
  my.ct[[expt]] = curate.tree(my.tree, df)
  
  plot.ct[[expt]] = plotHypercube.curated.tree(my.ct[[expt]])
}

sf = 2
png("new-rev-control-data.png", width=1000*sf, height=800*sf, res=72*sf)
ggarrange(plotlist=plot.ct)
dev.off()

# wrapper function for HyperTraPS analysis
parallel.fn = function(ref, seed) {
  return(HyperTraPS(my.ct[[ref]]$dests, initialstates = my.ct[[ref]]$srcs,
                    length = 5, kernel = 5,
                    seed = seed))
}

if(run.inference == TRUE) {
# run these experiments in parallel. should take a few core minutes each
n.seed = 2
parallelised.runs <- mcmapply(parallel.fn, 
                              ref=rep(1:6, each=n.seed),
                              seed=rep(1:n.seed, 6),
                              SIMPLIFY = FALSE,
                              mc.cores = min(detectCores(), 6*n.seed))

save(parallelised.runs, file="new-rev-control.Rdata")
} else {
  load("new-rev-control.Rdata")
}

png("new-rev-control-reps.png", width=1000*sf, height=800*sf, res=72*sf)
ggarrange(plotHypercube.bubbles.compare(parallelised.runs[1:2]),
          plotHypercube.bubbles.compare(parallelised.runs[3:4]),
          plotHypercube.bubbles.compare(parallelised.runs[5:6]),
          plotHypercube.bubbles.compare(parallelised.runs[7:8]),
          plotHypercube.bubbles.compare(parallelised.runs[9:10]),
          plotHypercube.bubbles.compare(parallelised.runs[11:12]))
dev.off()

png("new-rev-control-all.png", width=1000*sf, height=800*sf, res=72*sf)
ggarrange(plotHypercube.curated.tree(my.ct[[1]]),
          plotHypercube.bubbles.compare(parallelised.runs[1:2]) +
            theme(legend.position="none"),
          plotHypercube.curated.tree(my.ct[[2]]),
          plotHypercube.bubbles.compare(parallelised.runs[3:4]) +
            theme(legend.position="none"),
          plotHypercube.curated.tree(my.ct[[3]]),
          plotHypercube.bubbles.compare(parallelised.runs[5:6]) +
            theme(legend.position="none"),
          plotHypercube.curated.tree(my.ct[[4]]),
          plotHypercube.bubbles.compare(parallelised.runs[7:8]) +
            theme(legend.position="none"),
          plotHypercube.curated.tree(my.ct[[5]]),
          plotHypercube.bubbles.compare(parallelised.runs[9:10]) +
            theme(legend.position="none"),
          plotHypercube.curated.tree(my.ct[[6]]),
          plotHypercube.bubbles.compare(parallelised.runs[11:12]) +
            theme(legend.position="none"), ncol=4, nrow=3,
          labels=c("i","","ii","","iii","","iv","","v","","vi",""))
dev.off()

png("new-rev-control-singles.png", width=1000*sf, height=800*sf, res=72*sf)
ggarrange(plotHypercube.curated.tree(my.ct[[1]]),
          plotHypercube.bubbles(parallelised.runs[1]) +
            theme(legend.position="none"),
          plotHypercube.curated.tree(my.ct[[2]]),
          plotHypercube.bubbles(parallelised.runs[3]) +
            theme(legend.position="none"),
          plotHypercube.curated.tree(my.ct[[3]]),
          plotHypercube.bubbles(parallelised.runs[5]) +
            theme(legend.position="none"),
          plotHypercube.curated.tree(my.ct[[4]]),
          plotHypercube.bubbles(parallelised.runs[7]) +
            theme(legend.position="none"),
          plotHypercube.curated.tree(my.ct[[5]]),
          plotHypercube.bubbles(parallelised.runs[9]) +
            theme(legend.position="none"),
          plotHypercube.curated.tree(my.ct[[6]]),
          plotHypercube.bubbles(parallelised.runs[11]) +
            theme(legend.position="none"), ncol=4, nrow=3,
          labels=c("i","","ii","","iii","","iv","","v","","vi",""))
dev.off()
