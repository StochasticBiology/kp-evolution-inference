library(ggplot2)

if (file.exists("data/tree_metrics.Rdata")) {
  object <- load("data/tree_metrics.Rdata")
  if (object != "trees.df") stop("Invalid data/tree_metrics.Rdata")
} else {
  print("could not find tree_metrics.Rdata, generating it from scratch")
  devtools::load_all()
  trees.df <- get.tree.metrics.df()
  save(trees.df, file="data/tree_metrics.Rdata")
}


dir.create("plots/trees_metrics", showWarnings = FALSE, recursive = TRUE)
ggplot(trees.df, 
       aes(log(n_transitions), 
           reorder(country, n_transitions))) + 
  geom_bar(stat="identity") + 
  geom_text(aes(label=n_transitions), size=2, nudge_x=0.11) + 
  theme(text=element_text(size=6)) + 
  ggtitle("Ranked number of transitions in each newick tree", 
          subtitle="Bars in log scale")
ggsave("plots/trees_metrics/ranked_transitions.png", width=10, height=7)

ggplot(trees.df, aes(AGly_acquired, reorder(country, -AGly_acquired))) + 
  geom_bar(aes(fill=log(n_transitions)), stat="identity") + 
  theme(text=element_text(size=6)) + 
  geom_text(aes(label=round(AGly_acquired, digits=2)), size=2, nudge_x=0.11) + 
  ggtitle("Countries ordered by the mean number of features acquired before aminoglycosides",
          subtitle="Not acquired are designated by having mean Inf")
ggsave("plots/trees_metrics/agly.png", width=10, height=7)

ggplot(trees.df, aes(Bla_ESBL_acquired, reorder(country, -Bla_ESBL_acquired))) + 
  geom_bar(aes(fill=log(n_transitions)), stat="identity") + 
  theme(text=element_text(size=6)) + 
  geom_text(aes(label=round(Bla_ESBL_acquired, digits=2)), size=2, nudge_x=0.11) + 
  ggtitle("Countries ordered by the mean number of features acquired before ESBL",
          subtitle="Not acquired are designated by having mean Inf")
ggsave("plots/trees_metrics/bla_esbl.png", width=10, height=7)

ggplot(trees.df, aes(Flq_acquired, reorder(country, -Flq_acquired))) + 
  geom_bar(aes(fill=log(n_transitions)), stat="identity") + 
  theme(text=element_text(size=6)) + 
  geom_text(aes(label=round(Flq_acquired, digits=2)), size=2, nudge_x=0.11) + 
  ggtitle("Countries ordered by the mean number of features acquired before fluoroquinolones",
          subtitle="Not acquired are designated by having mean Inf")
ggsave("plots/trees_metrics/flq.png", width=10, height=7)

ggplot(trees.df, aes(Tmt_acquired, reorder(country, -Tmt_acquired)))+ 
  geom_bar(aes(fill=log(n_transitions)), stat="identity") + 
  theme(text=element_text(size=6)) + 
  geom_text(aes(label=round(Tmt_acquired, digits=2)), size=2, nudge_x=0.11) + 
  ggtitle("Countries ordered by the mean number of features acquired before trimethoprim",
          subtitle="Not acquired are designated by having mean Inf")
ggsave("plots/trees_metrics/tmt.png", width=10, height=7)