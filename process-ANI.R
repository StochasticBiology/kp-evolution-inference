library(phangorn)

# read summary ANI scores
df = read.csv("1-summary.txt", skip=1, sep=" ", header=FALSE)
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
  am[x,y] = am[y,x] = df$distance[i]
}

# get a tree by clustering
treeUPGMA = upgma(am)
