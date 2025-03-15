library(phangorn)
library(ape)
source("hypertraps.R")
library(ggtree)
options(ignore.negative.edge=TRUE)

# read summary ANI scores from all-all comparison (old and new)
df = read.csv("summary-all.txt", skip=1, sep=" ", header=FALSE)
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
tip.cols = ifelse(grepl("SAMN", treeNJ$tip.label), "red", "blue")

ggtree(treeNJ) +  geom_tiplab(size=3,color = tip.cols)
#ggtree(treeNJ) +  geom_tiplab(size=3,color = tip.cols) + layout_circular()

# pull the feature sets from the new dataset
tmpdf = read.csv("From_Olav_fixed/kleborate_new_tanzania_samples_output_2.csv")
idset = sapply(strsplit(tmpdf$id, "[.]"), `[`, 1)
f.df = data.frame(id = idset, tmpdf[,3:ncol(tmpdf)])
new.set = curate.tree(treeNJ, f.df). # XXX why doesn't this throw an error -- as we don't have the old features yet?

# get the old tree based on LIN codes
old.tree = read.tree("From_Olav_fixed/Tanzania.nwk")
# and the old dataset of features
old.refs = read.csv("From_Olav_fixed/tanzania-resistance-profiles.csv")
for(i in 1:length(old.tree$tip.label)) {
  r = which(old.refs$id == old.tree$tip.label[i])
  old.tree$tip.label[i] = old.refs$displayname[r]
}
plot(old.tree)
ggtree(old.tree) + geom_tiplab2() + layout_circular()

# combine the old and new feature sets
f2.df = old.refs[,3:ncol(old.refs)]
colnames(f2.df)[1] = "id"

fboth.df = rbind(f.df, f2.df)

# now curate the all tree with the all dataset
all.ct = curate.tree(treeNJ, fboth.df)
plotHypercube.curated.tree(all.ct)

all.ct$data[grepl("SAMN", all.ct$data$label),2:ncol(all.ct$data)] = 
  2*all.ct$data[grepl("SAMN", all.ct$data$label),2:ncol(all.ct$data)]

plotHypercube.curated.tree(all.ct, hjust=1) +  scale_y_continuous(expand = expansion(mult = c(0.2, 0.05)))

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


############# pull set of subtrees that have no SAMN (old) ancestry
# Load your tree (replace "your_tree.nwk" with your actual file)
tree <- treeNJ

# Identify tips that begin with "SAMN"
samn_tips <- grep("^SAMN", tree$tip.label, value = TRUE)

# Function to check if a clade contains only non-"SAMN" tips
is_valid_clade <- function(node, tree, samn_tips) {
  tips_in_clade <- extract.clade(tree, node)$tip.label
  return(!any(tips_in_clade %in% samn_tips))  # TRUE if no "SAMN" labels
}

# Get all internal nodes
internal_nodes <- (length(tree$tip.label) + 1):max(tree$edge)

# Find valid subtrees
valid_clades <- Filter(function(node) is_valid_clade(node, tree, samn_tips), internal_nodes)

# Extract the subtrees
subtrees <- lapply(valid_clades, function(node) extract.clade(tree, node))

# Print number of extracted subtrees
length(subtrees)

# (Optional) Save subtrees as newick files
#for (i in seq_along(subtrees)) {
#  write.tree(subtrees[[i]], file = paste0("subtree_", i, ".nwk"))
#}

# (Optional) Plot one of the subtrees
plot(subtrees[[12]], main = "Example Subtree")

# get the set of transitions we acquire from curating these individual subtrees
for(i in 1:length(subtrees)) {
  tmp.ct = curate.tree(subtrees[[i]], f.df) 
  if(i == 1) {
    all.srcs = tmp.ct$srcs
    all.dests = tmp.ct$dests
    } else {
      all.srcs = rbind(all.srcs, tmp.ct$srcs)
      all.dests = rbind(all.dests, tmp.ct$dests)
    }
}

#new.fit = HyperTraPS(all.dests, initialstates = all.srcs, 
#                     length = 5, kernel = 3,
#                     seed = 1)

# summarise the changes we see in this new-independent-of-old datasets
new.changes = which(all.srcs != all.dests, arr.ind = TRUE)
sub.srcs = all.srcs[unique(new.changes[,1]),]
sub.dests = all.dests[unique(new.changes[,1]),]
sets = sub.dests
sets[sub.srcs == 0 & sub.dests == 1] = 2
plotHypercube.curated.tree(new.set, hjust = 1, font.size=2) +  scale_y_continuous(expand = expansion(mult = c(0.2, 0.05)))
#new.fit = HyperTraPS(new.set$dests, initialstates = new.set$srcs, 
#           length = 5, kernel = 3,
#           seed = 1)
#new.fit$featurenames = colnames(f.df)[2:ncol(f.df)]
#plotHypercube.sampledgraph2(new.fit, node.labels = FALSE, no.times = TRUE, thresh=0.05, truncate = 6)

