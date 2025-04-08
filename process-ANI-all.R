library(phangorn)
library(ape)
source("hypertraps.R")
library(ggtree)
options(ignore.negative.edge=TRUE)

# read summary ANI scores from all-all comparison (old and new)
df = read.csv("138-summary.txt", skip=1, sep=" ", header=FALSE)
colnames(df) = c("isolate.1", "isolate.2", "ani")

# phrase as distances and initialise a distance matrix
df$distance = 1 - df$ani/100
ids = unique(df$isolate.1)
am = matrix(0, nrow=length(ids), ncol=length(ids))
rownames(am) = ids
colnames(am) = ids

# populate distance matrix from pairwise ANI scores
for(i in 1:nrow(df)) {
  x = which(ids == df$isolate.1[i])
  y = which(ids == df$isolate.2[i])
  am[x,y] = am[y,x] = df$distance[i]*10
}

# get a tree by clustering
treeNJ = upgma(am)
#treeNJ = NJ(am)

plot(treeNJ)

write.tree(treeNJ, "tree-all.phy")

# label old and new differently
tip.cols = rep("blue", length(treeNJ$tip.label))
tip.cols[grep("SAMN", treeNJ$tip.label)] = "red"
tip.cols[grep("2457", treeNJ$tip.label)] = "green"

ggtree(treeNJ) +  geom_tiplab(size=3,color = tip.cols)
ggtree(treeNJ) +  geom_tiplab(size=3,color = tip.cols) + layout_circular()

# pull the feature sets from the new dataset
tmpdf = read.csv("../From_Olav_fixed/dichotomized_1_3_8.csv")
tmpdf = cbind(tmpdf$strain, tmpdf)
colnames(tmpdf)[c(1,2)] = c("X", "id")
idset = sapply(strsplit(tmpdf$id, "[.]"), `[`, 1)
f.df = data.frame(id = idset, tmpdf[,3:ncol(tmpdf)])
new.set = curate.tree(treeNJ, f.df) # XXX why doesn't this throw an error -- as we don't have the old features yet?

# get the old tree based on LIN codes
old.tree = read.tree("../From_Olav_fixed/Tanzania.nwk")
# and the old dataset of features
old.refs = read.csv("../From_Olav_fixed/tanzania-resistance-profiles.csv")
for(i in 1:length(old.tree$tip.label)) {
  r = which(old.refs$id == old.tree$tip.label[i])
  old.tree$tip.label[i] = old.refs$displayname[r]
}
plot(old.tree)
ggtree(old.tree) + geom_tiplab2() + layout_circular()

# curate data for the old dataset (should match Olav's input data)
old.ct = curate.tree(old.tree, old.refs[,3:ncol(old.refs)])
plotHypercube.curated.tree(old.ct)

# do the inference (should match Olav's output)
old.fit = HyperTraPS(old.ct$dests, initialstates = old.ct$srcs, 
                     length = 5, kernel = 3, walkers = 400,
                     seed = 1)
old.fit$featurenames = colnames(f.df)[2:ncol(f.df)]
plotHypercube.sampledgraph2(old.fit, node.labels = FALSE, no.times = TRUE, thresh=0.05, truncate = 6)

# combine the old and new feature sets
f2.df = old.refs[,3:ncol(old.refs)]
colnames(f2.df)[1] = "id"

fboth.df = rbind(f.df, f2.df)

# now curate the all tree with the all dataset
all.ct = curate.tree(treeNJ, fboth.df)
plotHypercube.curated.tree(all.ct)

all.ct$data[grepl("SAMN", all.ct$data$label),2:ncol(all.ct$data)] = 
  3*all.ct$data[grepl("SAMN", all.ct$data$label),2:ncol(all.ct$data)]

all.ct$data[grepl("2457", all.ct$data$label),2:ncol(all.ct$data)] = 
  2*all.ct$data[grepl("2457", all.ct$data$label),2:ncol(all.ct$data)]

all.data.plot = plotHypercube.curated.tree(all.ct, hjust=1, font.size = 2) +  
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.05)))

all.data.plot
### need to look into this Rif_acquired behaviour

# some comparisons of ANI to LIN code trees
mean(am)
# very highly related by LIN
am["SAMN10390519","SAMN10390521"]
am["SAMN10390498","SAMN10390483"]
am["SAMN10390455","SAMN10390461"]
# moderate LIN
am["SAMN10390519","SAMN10390463"]
am["SAMN10390537","SAMN10390446"]
am["SAMN10390525","SAMN10390453"]
# distant LIN
am["SAMN10390519","SAMN10390542"]
am["SAMN10390499","SAMN10390526"]
am["SAMN10390509","SAMN10390515"]

# get sets of nodes that only have new data as descendants
n.tips = length(treeNJ$tip.label) 
safe.new = c()
for(i in (n.tips+1):(n.tips+treeNJ$Nnode)) {
  kids = Descendants(treeNJ, i, type =  "tips")[[1]]
  if(all(grepl("_", treeNJ$tip.label[kids]))) {
    safe.new = c(safe.new, i)
  }
}
# collect set of transitions that involve these nodes (independent new transitions)
safe.trans = c()
for(i in 1:nrow(all.ct$transitions)) {
  trans = all.ct$transitions[i,]
  if(trans$from.node %in% safe.new | trans$to.node %in% safe.new) {
    safe.trans = c(safe.trans, i)
  }
}

# XXX NEXT DO PREDICTIONS FOR DIFFERENT TIMES OF DATA

# construct predictions from (differently trained!) inference, and test with respect to new independent transitions
freqs = colSums(old.refs[,4:ncol(old.refs)])/nrow(old.refs)
preds.ranks = data.frame()
for(j in 1:length(safe.trans)) {
  changes = which(all.ct$srcs[safe.trans[j],] != all.ct$dests[safe.trans[j],])
  already = which(all.ct$srcs[safe.trans[j],] == 1)
  if(length(changes) > 0) {
    # if this is a non-trivial transition
    predict.fit = PosteriorAnalysis(old.fit, startstate = all.ct$srcs[safe.trans[j],])
    p.rates = predict.fit$predictrates
    n.poss = length(p.rates) - length(already)
    # rank predicted possible next steps 
    meaningful = p.rates
    meaningful[already] = NA
    meaningful.ranks = rank(meaningful, na.last = "keep")
    changes.ranks = (1+n.poss-meaningful.ranks)[changes]
    preds.ranks = rbind(preds.ranks, data.frame(ref=j, model="trained", ranks=changes.ranks))
    # do this without trained model just based on frequencies
    untrained = freqs
    untrained[already] = NA
    untrained.ranks = rank(untrained, na.last = "keep")
    changes.ranks.untrained = (1+n.poss-untrained.ranks)[changes]
    preds.ranks = rbind(preds.ranks, data.frame(ref=j, model="untrained", ranks=changes.ranks.untrained))
    # do this with null model
    null = rep(1, length(freqs)) #runif(length(freqs))
    null[already] = NA
    null.ranks = rank(null, na.last = "keep")
    changes.ranks.null = (1+n.poss-null.ranks)[changes]
    preds.ranks = rbind(preds.ranks, data.frame(ref=j, model="null", ranks=changes.ranks.null))
  }
}

#ggplot(preds.ranks, aes(x=ranks)) + geom_histogram() + facet_wrap(~ model, nrow=3)
# ^ well, not hugely impressive

predict.plot = ggplot(preds.ranks[preds.ranks$model != "untrained",], aes(x=ranks, color=model, fill=model)) + 
  geom_histogram(aes(y=..density..), binwidth=1, position="dodge")+
  geom_density(alpha=0.4) + 
  labs(x = "Predicted rank of true next steps", y = "Probability", fill="Model") + 
  scale_color_discrete(guide="none") +
  theme_minimal()
predict.plot

sf = 3
png("predictions-138.png", width=800*sf, height=400*sf, res=72*sf)
ggarrange(all.data.plot, predict.plot, labels=c("A", "B"))
dev.off()

ggplot(preds.ranks, aes(x=ranks, fill=model)) + 
  geom_density(alpha=0.4) #+ geom_histogram(aes(y=..density..)) #+ facet_wrap(~ model, nrow=3)


# do the inference and produce some summary fits
new.fit = HyperTraPS(all.ct$dests[safe.trans,], initialstates = all.ct$srcs[safe.trans,], 
                     length = 5, kernel = 3,
                     seed = 1)
new.fit$featurenames = colnames(f.df)[2:ncol(f.df)]
plotHypercube.sampledgraph2(new.fit, node.labels = FALSE, no.times = TRUE, thresh=0.05, truncate = 6)

plotHypercube.bubbles(new.fit, featurenames=new.fit$featurenames)



