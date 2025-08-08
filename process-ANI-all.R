library(phangorn)
library(ape)
library(hypertrapsct)
library(igraph)
library(stringr)
library(ggtree)
library(ggraph)
library(ggplot2)
library(ggpubr)
library(dplyr)

run.inference = FALSE
options(ignore.negative.edge=TRUE)

sf = 2

# read summary ANI scores from all-all comparison (old and new)
df = read.csv("data/138-summary.txt", skip=1, sep=" ", header=FALSE)
colnames(df) = c("isolate.1", "isolate.2", "ani")
df.mis = read.csv("data/138-bases.txt", skip=1, sep=" ", header=FALSE)
colnames(df.mis) = c("isolate.1", "isolate.2", "bases")

# "2457" = old
# "SAMN" = Kleborate
# "247"etc = new

df.mis = cbind(df.mis, data.frame(ani=df$ani))
df.mis$class = 0
df.mis$class[grep("SAMN", df.mis$isolate.1)] = 1+df.mis$class[grep("SAMN", df.mis$isolate.1)] 
df.mis$class[grep("SAMN", df.mis$isolate.2)] = 1+df.mis$class[grep("SAMN", df.mis$isolate.2)] 
df.mis$class[grep("2457", df.mis$isolate.1)] = 10+df.mis$class[grep("2457", df.mis$isolate.1)] 
df.mis$class[grep("2457", df.mis$isolate.2)] = 10+df.mis$class[grep("2457", df.mis$isolate.2)] 

df.mis$class[df.mis$class == 0] = "New-new"
df.mis$class[df.mis$class == 1] = "Kleb-new"
df.mis$class[df.mis$class == 2] = "Kleb-Kleb"
df.mis$class[df.mis$class == 10] = "New-old"
df.mis$class[df.mis$class == 11] = "Kleb-old"
df.mis$class[df.mis$class == 20] = "Old-old"

diag.1 = ggplot(df.mis, aes(x=bases, y=ani, color=factor(class))) + 
  geom_point(size=1) + theme_minimal() +
  labs(x = "% bases unaligned", y = "ANI among aligned bases")
diag.1

big.base = unique(df.mis$isolate.1[df.mis$bases>25])
unique(df.mis$isolate.1[df.mis$ani<95])

diag.2 = ggplot(df.mis[!(df.mis$isolate.1 %in% big.base),], aes(x=bases, y=ani, color=factor(class))) + 
  geom_point(size=1) + theme_minimal() +
  labs(x = "% bases unaligned", y = "ANI among aligned bases")
diag.2

png("diagnosis-plot.png", width=600*sf, height=600*sf, res=72*sf)
ggarrange(diag.1, diag.2)
dev.off()

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

vert.tree = ggtree(treeNJ) +  geom_tiplab(size=3,color = tip.cols)
png("tree-all.png", width=800*sf, height=2000*sf, res=72*sf)
print(vert.tree)
dev.off()

circ.tree = ggtree(treeNJ) +  geom_tiplab(aes(label=label,angle=angle), size=2,color = tip.cols) + 
  layout_circular() +
  theme(
    plot.margin = unit(c(3, 0, 3, 0), "cm")  # top, right, bottom, left
  )

png("tree-all-circ.png", width=800*sf, height=800*sf, res=72*sf)
print(circ.tree)
dev.off()

# pull the feature sets from the new dataset
tmpdf = read.csv("data/dichotomized_1_3_8.csv")
tmpdf = cbind(tmpdf$strain, tmpdf)
colnames(tmpdf)[c(1,2)] = c("X", "id")
idset = sapply(strsplit(tmpdf$id, "[.]"), `[`, 1)
f.df = data.frame(id = idset, tmpdf[,3:ncol(tmpdf)])
new.set = curate.tree(treeNJ, f.df) # XXX why doesn't this throw an error -- as we don't have the old features yet?

# get the old tree based on LIN codes
old.tree = read.tree("data/Tanzania.nwk")
# and the old dataset of features
old.refs = read.csv("data/tanzania-resistance-profiles.csv")
for(i in 1:length(old.tree$tip.label)) {
  r = which(old.refs$id == old.tree$tip.label[i])
  old.tree$tip.label[i] = old.refs$displayname[r]
}
plot(old.tree)
ggtree(old.tree) + geom_tiplab2() + layout_circular()

# curate data for the old dataset (should match Olav's input data)
old.ct = curate.tree(old.tree, old.refs[,3:ncol(old.refs)])
plotHypercube.curated.tree(old.ct)

plot.old = plotHypercube.curated.tree(old.ct)

# do the inference (should match Olav's output)
if(run.inference == TRUE) {
old.fit = HyperTraPS(old.ct$dests, initialstates = old.ct$srcs, 
                     length = 5, kernel = 3, walkers = 400,
                     seed = 1)
  save(old.fit, file="fitted-analysed-138-fit.Rdata")
} else {
  load("fitted-analysed-138-fit.RData")
}
old.fit$featurenames = colnames(f.df)[2:ncol(f.df)]
plotHypercube.sampledgraph2(old.fit, node.labels = FALSE, no.times = TRUE, thresh=0.05, truncate = 6)

plot.old.graph = plotHypercube.sampledgraph2(old.fit, node.labels = FALSE, 
                                             no.times = TRUE, thresh=0.025, truncate = 6)

plot.old.bubbles = plotHypercube.bubbles(old.fit) #+ scale_y_discrete(labels = colnames(f.df)[2:ncol(f.df)])
ggarrange(plot.old,
          ggarrange(plot.old.bubbles, plot.old.graph, nrow=2, labels=c("B", "C")),
          nrow=1, labels = c("A", ""))

# combine the old and new feature sets
f2.df = old.refs[,3:ncol(old.refs)]
colnames(f2.df)[1] = "id"

fboth.df = rbind(f.df, f2.df)

# now curate the all tree with the all dataset
all.ct = curate.tree(treeNJ, fboth.df)
plotHypercube.curated.tree(all.ct)

# "2457" = old
# "SAMN" = Kleborate
# "247"etc = new

old.sabrina = fboth.df[grep("^2457", fboth.df$id),]
new.sabrina = fboth.df[-c(grep("^2457", fboth.df$id),grep("SAMN", fboth.df$id)),]
kleb.df = fboth.df[grep("SAMN", fboth.df$id),]
old.sabrina.ct = curate.tree(treeNJ, old.sabrina)
new.sabrina.ct = curate.tree(treeNJ, new.sabrina)
kleb.df.ct = curate.tree(treeNJ, kleb.df)

ct.plots = ggarrange(plotHypercube.curated.tree(old.sabrina.ct, hjust = 1, font.size = 2) +
            coord_cartesian(clip = "off") + theme(
              plot.margin = unit(c(0, 0, 3, 0), "cm")  # top, right, bottom, left
            ),
          plotHypercube.curated.tree(new.sabrina.ct, hjust=1, font.size = 2)+
            coord_cartesian(clip = "off") + theme(
              plot.margin = unit(c(0, 0, 3, 0), "cm")  # top, right, bottom, left
            ),
          plotHypercube.curated.tree(kleb.df.ct, hjust=1, font.size = 2)+
            coord_cartesian(clip = "off") + theme(
              plot.margin = unit(c(0, 0, 3, 0), "cm")  # top, right, bottom, left
            ), 
          nrow=1, labels=c("A", "B", "C"))

get_proportions <- function(df, name) {
  data.frame(
    column = names(df),
    proportion = colMeans(df == 1),
    dataset = name
  )
}
df_summary <- bind_rows(
  get_proportions(old.sabrina, "2001"),
  get_proportions(new.sabrina, "2017"),
  get_proportions(kleb.df, "Pathogenwatch")
)


comp.plot = ggplot(df_summary[df_summary$column != "id",], aes(x=column, y=proportion, fill=dataset)) + 
  geom_col(position="dodge", width=0.65) +  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1)) +
   scale_fill_manual(values=c("#8888FF", "#880000", "#88888855")) + 
  labs(x="Character", y="Proportion\nwith character", fill = "Dataset")

png("new-data-comp.png", width=600*sf, height=600*sf, res=72*sf)
ggarrange( ct.plots, comp.plot, heights = c(2, 1), labels=c("", "D"), nrow=2)
dev.off()

zanzibar.df = read.csv("data/zanzibar-dichotomized.csv")

zanzibar.tree <- curate.tree("data/zanzibar-4-tree-1.phy", "data/zanzibar-dichotomized.csv")
png("new-data-zanzibar.png", width=600*sf, height=600*sf, res=72*sf)
plotHypercube.curated.tree(zanzibar.tree, hjust=1, font.size = 3) +
  coord_cartesian(clip = "off") + theme(
    plot.margin = unit(c(0, 0, 3, 0), "cm")  # top, right, bottom, left
  )
dev.off()


ct.plots.z = ggarrange(plotHypercube.curated.tree(zanzibar.tree, hjust=1, font.size = 2) +
                         coord_cartesian(clip = "off") + theme(
                           plot.margin = unit(c(0, 0, 3, 0), "cm")  # top, right, bottom, left
                         ),
                       plotHypercube.curated.tree(old.sabrina.ct, hjust = 1, font.size = 2) +
                       coord_cartesian(clip = "off") + theme(
                         plot.margin = unit(c(0, 0, 3, 0), "cm")  # top, right, bottom, left
                       ),
                     plotHypercube.curated.tree(new.sabrina.ct, hjust=1, font.size = 2)+
                       coord_cartesian(clip = "off") + theme(
                         plot.margin = unit(c(0, 0, 3, 0), "cm")  # top, right, bottom, left
                       ),
                     plotHypercube.curated.tree(kleb.df.ct, hjust=1, font.size = 2)+
                       coord_cartesian(clip = "off") + theme(
                         plot.margin = unit(c(0, 0, 3, 0), "cm")  # top, right, bottom, left
                       ), 
                     ncol=2, nrow=2, heights=c(1,2), labels=c("A", "B", "C", "D"))
ct.plots.z

df_summary.z <- bind_rows(
  get_proportions(old.sabrina, "2001"),
  get_proportions(new.sabrina, "2017"),
  get_proportions(kleb.df, "Pathogenwatch"),
  get_proportions(zanzibar.df, "Zanzibar")
)

comp.plot.z = ggplot(df_summary.z[!(df_summary.z$column %in% c("id","strain")),], aes(x=column, y=proportion, fill=dataset)) + 
  geom_col(position="dodge", width=0.7) +  theme_minimal() + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("#88CCFF", "#880000", "#CCCCCC", "#FF8800", "#005500")) + 
  labs(x="Character", y="Proportion\nwith character", fill = "Dataset")

png("new-data-comp-z.png", width=600*sf, height=800*sf, res=72*sf)
ggarrange( ct.plots.z, comp.plot.z, heights = c(3, 1), labels=c("", "E"), nrow=2)
dev.off()

# SAMN = Kleborate = highest multiplier (x3)
# 2457 = old = next multiplier (x2)
# remaining 247 = new = no multiplier (1)

all.ct$data[grepl("^247", all.ct$data$label),2:ncol(all.ct$data)] = 
  3*all.ct$data[grepl("^247", all.ct$data$label),2:ncol(all.ct$data)]

all.ct$data[grepl("^2457", all.ct$data$label),2:ncol(all.ct$data)] = 
  2*all.ct$data[grepl("^2457", all.ct$data$label),2:ncol(all.ct$data)]

all.data.plot = plotHypercube.curated.tree(all.ct, hjust=1, font.size = 2) +  
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.05)))

# so new = red, old = blue, Kleborate = grey
all.data.plot = plotHypercube.curated.tree(all.ct, factor.vals = TRUE, hjust = 1, font.size=2) +
  scale_fill_manual(values = c("white", "grey", "#4444FF", "#FF8888")) +  
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.05)))


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
png("predictions-138-new.png", width=500*sf, height=600*sf, res=72*sf)
ggarrange(ggarrange(all.data.plot, predict.plot, widths=c(1.,1), nrow=1, labels=c("A", "B")),
          comp.plot.z, labels=c("", "C"), nrow=2, heights=c(1.5,1))
dev.off()

png("new-data-summaries.png", width=1000*sf, height=600*sf, res=72*sf)
ggarrange(ct.plots.z, circ.tree, nrow=1, labels=c("", "E"), widths=c(1,1.3))
dev.off()

