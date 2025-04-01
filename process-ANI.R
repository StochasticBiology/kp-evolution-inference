library(phangorn)
source("hypertraps.R")

experiment = "4"

# read summary ANI scores
df = read.csv(paste0(experiment, "-summary.txt", collapse=""), skip=1, sep=" ", header=FALSE)
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
treeUPGMA = upgma(am)
treeNJ = NJ(am)

plot(treeUPGMA)
plot(treeNJ)

write.tree(treeNJ, paste0(experiment, "-tree.phy", collapse=""))

if(experiment == "1") {
  tmpdf = read.csv("From_Olav_fixed/kleborate_new_tanzania_samples_output_2.csv") 
} else {
  tmpdf = read.csv("From_Olav_fixed/kleborate_new_tanzania_samples_output_1.csv") 
}
idset = sapply(strsplit(tmpdf$id, "[.]"), `[`, 1)
f.df = data.frame(id = idset, tmpdf[,3:ncol(tmpdf)])
new.set = curate.tree(treeNJ, f.df)

plotHypercube.curated.tree(new.set, hjust = 1, font.size=2) +  scale_y_continuous(expand = expansion(mult = c(0.2, 0.05)))
new.fit = HyperTraPS(new.set$dests, initialstates = new.set$srcs, 
           length = 5, kernel = 3,
           seed = 1)
new.fit$featurenames = colnames(f.df)[2:ncol(f.df)]
plotHypercube.sampledgraph2(new.fit, node.labels = FALSE, no.times = TRUE, thresh=0.1)

