library(ape)
library(phangorn)
library(parallel)
source("hypertraps.R")

# how much of a problem is reversibility?
set.seed(1)

# we will construct a dataset synthetically with known dynamics, then allow reversible transitions and try to infer the true behaviour
L = 10
lambda.rev = 0.5
n.rev = L/2

# parameterisation for tree construction
tree.size = 128
birth.rate = 1
death.rate = 0.5

# accumulation rate for features (and loss rate, for reversible setup)
accumulation.rate = 1
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
      # find leftmost zero in current state, and change with some probability
      # recall dynamics here are 00000 -> 10000 -> 11000 -> 11100 -> 11110 -> 11111
      ## (see first paragraph of section "Synthetic case studies")
      ref = which(x[[this.child]] == 0)[1]
      if(runif(1) < accumulation.rate*this.branch.length) { x[[this.child]][ref] = 1 } 
      x.loss[[this.child]] = x[[this.child]]
      # in the reversible case, allow the leftmost feature ("first feature" in the ms.:
      # second paragraph of "Synthetic case studies") to revert with some probability
      if(FALSE) {
      for(locus in 1:L) {
        if(runif(1) < loss.rate*this.branch.length) { x[[this.child]][locus] = 0 }
      }
      }
      # add this child to to state list, and to next iteration's to-do
      new.to.do = c(new.to.do, this.child)
    }
  }
  # update to-do list
  to.do = new.to.do
}

# keep the original (irreversible) picture in x.true, while allowing x to have reversibility
x.true = x

# model reversibility by reverting the first several loci randomly
for(i in 1:length(my.tree$tip.label)) {
  n = rpois(1,lambda.rev)
  to.0 = sample(1:n.rev, n)
  x[[i]][to.0] = 0
}

# construct the dataframes of binary labels for tree tips
df.true = data.frame()
for(i in 1:length(my.tree$tip.label)) {
  df.true = rbind(df.true, data.frame(
    cbind(data.frame(id=my.tree$tip.label[i]), matrix(x.true[[i]], ncol=L))
  ))
}
# construct the dataframe of binary labels for tree tips
df = data.frame()
for(i in 1:length(my.tree$tip.label)) {
  df = rbind(df, data.frame(
    cbind(data.frame(id=my.tree$tip.label[i]), matrix(x[[i]], ncol=L))
  ))
}

# HyperTraPS pipeline
my.ct = curate.tree(my.tree, df)
my.ct.true = curate.tree(my.tree, df.true)

plot.ct = plotHypercube.curated.tree(my.ct)
plot.ct.true = plotHypercube.curated.tree(my.ct.true)

ggarrange(plot.ct, plot.ct.true)

# wrapper function for HyperTraPS analysis
parallel.fn = function(seed) {
  if(seed <= 3) {
  return(HyperTraPS(my.ct$dests, initialstates = my.ct$srcs,
                    length = 5, kernel = 5,
                    seed = seed))
  } else {
    return(HyperTraPS(my.ct.true$dests, initialstates = my.ct.true$srcs,
               length = 5, kernel = 5,
               seed = seed))
    }
}

# run these experiments in parallel. should take a few core minutes each
n.seed = 6
parallelised.runs <- mcmapply(parallel.fn, seed=1:n.seed,
                              SIMPLIFY = FALSE,
                              mc.cores = min(detectCores(), n.seed))

plotHypercube.bubbles.compare(parallelised.runs[1:3])
plotHypercube.bubbles.compare(parallelised.runs[4:6])

if(FALSE) {
h.inf = HyperTraPS(my.ct$dests, initialstates = my.ct$srcs,
                   length = 5, kernel = 5,
                   seed = 1, samplegap = 100)
h.inf.true = HyperTraPS(my.ct.true$dests, initialstates = my.ct.true$srcs,
                   length = 5, kernel = 5,
                   seed = 1, samplegap = 100)
}

h.inf = parallelised.runs[[1]]
h.inf.true = parallelised.runs[[4]]

plot.bubbles = plotHypercube.bubbles(h.inf)
plot.bubbles.true = plotHypercube.bubbles(h.inf.true)

ggarrange(plot.ct.true, plot.bubbles.true,
          plot.ct, plot.bubbles,
          nrow = 2, ncol = 2, widths = c(2,1))

#plotHypercube.sampledgraph2(h.inf, no.times = TRUE, node.labels = FALSE)

