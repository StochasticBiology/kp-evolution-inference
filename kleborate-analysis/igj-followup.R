library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(ggpubr)
library(ggupset)
library(countrycode)
library(hypertrapsct)
library(tidyr)
library(dplyr)
library(viridisLite)
library(phytools)

# read in results of inference
if (!exists("country.list")) {
  name <- load("all_models.Rdata")
  country.list <- get(name)
}
# scaling factor for graphics
sf = 2

# broad contents
# TASK 1 -- preprocessing and curating
# TASK 2 -- PCA plots of summarised "bubble plot" dynamics (Fig 2B)
# TASK 3 -- bubble plots corresponding to PCA extremes (Fig 2C)
# TASK 4 -- individual features and PCA behaviour (Fig 3)
# TASK 5 -- geographical regions and PCA stats (Fig 4A-C)
# TASK 6 -- drug use data and acquisition (Fig 4D)
# TASK 7 -- summary global data (Fig 1A-B, Fig 2A)
# TASK 8 -- case study plot (Fig 1C)

############### TASK 1 -- preprocessing and curating

# read in prepared collection of bubble plots, and process
df = read.csv("igj-all-bubbles.csv", sep=";", stringsAsFactors = FALSE)
df = df[df$country != "USA2" & df$country != "USA3",]
region <- read.csv("misc-data/data-Ib27t.csv") # gbd superregion
region$oldval = region$val
region$val = gsub("Central Europe, Eastern Europe, and Central Asia", "CEEECA", region$val)
region$val = gsub("High-income", "HI", region$val)
region$val = gsub("Latin America and Caribbean", "LAC", region$val)
region$val = gsub("North Africa and Middle East", "NAME", region$val)
region$val = gsub("South Asia", "SA", region$val)
region$val = gsub("Southeast Asia, East Asia, and Oceania", "SEEAO", region$val)
region$val = gsub("Sub-Saharan Africa", "SSAf", region$val)

# process feature names, abbreviating
tmp = df[1:22,]
feature.names = tmp$Name[order(tmp$OriginalIndex)]
feature.names = gsub("_acquired", "-a", feature.names)
feature.names = gsub("_mutations", "-m", feature.names)

# list of countries
countries <- unique(df$country)

# pull the bubble plots into a wide format
wide_df <- df %>%
  mutate(ProbTime = Probability * Time) %>%
  group_by(country, OriginalIndex) %>%
  summarise(SumProbTime = sum(ProbTime, na.rm = TRUE), 
            .groups = "drop") %>%
  pivot_wider(names_from = OriginalIndex, values_from = SumProbTime)

# build dataframe of country name conversions
gbd.map = data.frame(country=NULL, val=NULL)
gbd.map <- rbind(gbd.map, c("Burkina Faso","Burkina_Faso"))
gbd.map <- rbind(gbd.map, c("Czech Republic","Czech_Republic"))
gbd.map <- rbind(gbd.map, c("Arab Republic of Egypt","Egypt"))
gbd.map <- rbind(gbd.map, c("The Gambia","Gambia"))
gbd.map <- rbind(gbd.map, c("Hong Kong","Hong_Kong"))
gbd.map <- rbind(gbd.map, c("Islamic Republic of Iran","Iran"))
gbd.map <- rbind(gbd.map, c("Lao People's Democratic Republic","Laos"))
gbd.map <- rbind(gbd.map, c("New Zealand","New_Zealand"))
gbd.map <- rbind(gbd.map, c("West Bank and Gaza","Palestine"))
gbd.map <- rbind(gbd.map, c("Russian Federation","Russia"))
gbd.map <- rbind(gbd.map, c("Saudi Arabia","Saudi_Arabia"))
gbd.map <- rbind(gbd.map, c("Slovak Republic","Slovakia"))
gbd.map <- rbind(gbd.map, c("South Africa","South_Africa"))
gbd.map <- rbind(gbd.map, c("Republic of Korea","South_Korea"))
gbd.map <- rbind(gbd.map, c("United Arab Emirates","United_Arab_Emirates"))
gbd.map <- rbind(gbd.map, c("R. B. de Venezuela","Venezuela"))
gbd.map <- rbind(gbd.map, c("United Kingdom","United_Kingdom"))
gbd.map <- rbind(gbd.map, c("United States of America", "USA"))

# assign GBD regions
region = rbind(region, data.frame(country="Hong Kong", val="SEEAO", oldval="Southeast Asia, East Asia, and Oceania"))
region = rbind(region, data.frame(country="Guadeloupe", val="LAC", oldval="Latin America and Caribbean"))
region = rbind(region, data.frame(country="Zanzibar", val="SSAf", oldval="Sub-Saharan Africa"))
for(i in 1:nrow(region)) {
  ref = which(gbd.map[,1] == region$country[i])
  if(length(ref) > 0) {
    region$country[i] = gbd.map[ref,2]
  }
}

# pull regions into dataframe
wide_df$region = ""
for(i in 1:nrow(wide_df)) {
  ref = which(region$country == wide_df$country[i])
  if(length(ref) > 0) {
    wide_df$region[i] = region$val[ref]
  }
}

# pull country codes into dataframe
wide_df$ccode = countrycode(wide_df$country, origin="country.name", destination = "iso3c")
wide_df$ccode[wide_df$country == "Zanzibar"] = "TZA*"

probtime_matrix <- as.matrix(wide_df[,-c(1, (ncol(wide_df))-0:1)])  # Remove 'country' column
rownames(probtime_matrix) <- wide_df$country

############### TASK 2 -- PCA plots of summarised "bubble plot" dynamics

# do the PCA
pca_result <- prcomp(probtime_matrix, scale. = TRUE)
pve <- (pca_result$sdev)^2 / sum(pca_result$sdev^2) * 100

# set up a dataframe storing PCA results
pca.df = data.frame(country=rownames(probtime_matrix),
                    region="",
                    pca1=pca_result$x[,1], 
                    pca2=pca_result$x[,2], 
                    pca3=pca_result$x[,3], 
                    pca4=pca_result$x[,4], 
                    pca5=pca_result$x[,5])
pca.df$ccode = countrycode(pca.df$country, origin="country.name", destination = "iso3c")
for(i in 1:nrow(pca.df)) {
  ref = which(region$country == pca.df$country[i])
  if(length(ref) > 0) {
    pca.df$region[i] = region$val[ref]
  }
}

# start producing the graphics
# useful colour scheme for the ellipses
pca.cols = viridis(7, option = "magma")
pca.cols[7] = "#FF0000"
pca.cols[1] = "#888888"

# PCA1 vs PCA2 plot, plus ellipses 
g.pca1.alt =  ggplot(pca.df, aes(x=pca1, y=pca2, color=region, label=country, fill=region)) + 
  geom_point(size = 1, alpha = 0.9, shape = 21, stroke = 0.5) +
  geom_point(data=pca.df[pca.df$country=="Zanzibar",], aes(x=pca1, y=pca2, label=country, fill=region), size=2, pch = 21, color="black") +
  stat_ellipse(level = 0.95, alpha = 0.5, linewidth = 1.2) +
  geom_text_repel(size=2.5, max.overlaps = 20, alpha=0.5, segment.alpha = 0.3)   + 
  labs(x="PCA1", y="PCA2", fill = "GBD Region", color = "GBD Region") +
  scale_color_manual(values = pca.cols) + scale_fill_manual(values = pca.cols) +
  theme_minimal()
g.pca1.alt

# PCA2 vs PCA3 plot, plus ellipses 
g.pca2.alt =  ggplot(pca.df, aes(x=pca3, y=pca2, color=region, label=country, fill=region)) + 
  geom_point(size = 1, alpha = 0.9, shape = 21, stroke = 0.5) +
  geom_point(data=pca.df[pca.df$country=="Zanzibar",], aes(x=pca1, y=pca2, label=country, fill=region), size=2, pch = 21, color="black") +
  stat_ellipse(level = 0.95, alpha = 0.5, linewidth = 1.2) +
  geom_text_repel(size=2.5, max.overlaps = 20, alpha=0.5, segment.alpha = 0.3)   + 
  labs(x="PCA3", y="PCA2", fill = "GBD Region", color = "GBD Region") +
  scale_color_manual(values = pca.cols) + scale_fill_manual(values = pca.cols) +
  theme_minimal()
g.pca2.alt

png("pca-alts.png", width=400*sf, height=600*sf, res=72*sf)
ggarrange(g.pca1.alt + theme(legend.position="none"), g.pca2.alt + theme(legend.position="bottom"), heights=c(1,1.3), nrow=2)
dev.off()

all.df = cbind(wide_df, pca.df[,2:ncol(pca.df)])

############### TASK 3 -- bubble plots corresponding to PCA extremes

# function for bubble plot for a particular country
b.plot = function(country, y.ticks="left") {
  sub = df[df$country == country,]
  sub$Name = gsub("_mutations", "-m", sub$Name)
  sub$Name = gsub("_acquired", "-a", sub$Name)
  ggplot(sub[sub$Probability > 1/22,], aes(x=Time, y=Name, size=Probability, color=Time)) + 
    geom_point(shape=16) + theme_minimal() + theme(legend.position = "none") +
    scale_colour_gradientn(
      colours = c("darkblue", "blue", "red", "darkred"),
      values = scales::rescale(c(0, 7, 14, 21))
    ) + labs(y=NULL, x="Ordinal Time") +
    scale_y_discrete(position=y.ticks)
}

# function to highlight particular rows in a bubble plot
highlight_layer = function(p, yvals, alpha=0.35) {
  highlight_data = data.frame(
    ymin = yvals-0.5,  # start of the highlighted y-ranges
    ymax = yvals+0.5   # end of the highlighted y-ranges
  )
  hlayer = geom_rect(
    data = highlight_data,
    aes(ymin = ymin, ymax = ymax),
    xmin = -Inf, xmax = Inf,
    fill = "yellow",
    alpha = alpha,
    inherit.aes = FALSE
  )
  p$layers <- c(list(hlayer), p$layers)
  return(p)
}

# create bubble plots for the particular countries of interest
bplot.1 = b.plot("Venezuela")+xlab(NULL)
bplot.2 = b.plot("Nigeria", y.ticks="right")+xlab(NULL)
bplot.3 = b.plot("South_Korea")+xlab(NULL)
bplot.4 = b.plot("Gambia", y.ticks="right")+xlab(NULL)

# combine into a plot object
b.plots = ggarrange(bplot.1 + xlab("Venezuela\n(low PCA1)"), ggplot()+theme_void(), 
                    bplot.2 + xlab("Nigeria\n(high PCA1)"), 
                    highlight_layer(bplot.3, c(2, 12, 15, 20), alpha=0.5) + xlab("South Korea\n(low PCA2)"), ggplot()+theme_void(),
                    highlight_layer(bplot.4, c(2, 12, 15, 20), alpha=0.5) + xlab("Gambia\n(high PCA2)"), ncol=3, nrow = 2, widths=c(1,0.2,1))

# produce plots
png("bplots.png", width=500*sf, height=500*sf, res=72*sf)
b.plots
dev.off()

############### TASK 4 -- individual features and PCA behaviour

# create dataframe with variables and PCA values
long_df <- all.df %>%
  pivot_longer(
    cols = -c("country", "ccode", "region", "pca1", "pca2", "pca3", "pca4", "pca5"),                     # all columns except 'country'
    names_to = "Variable",         # name for new column holding former column names
    values_to = "Value"           # name for new column holding the values
  )
long.df = as.data.frame(long_df)
long.df$Variable = as.numeric(long.df$Variable)
long.df$name = feature.names[long.df$Variable+1]

# subset particular features
set.1 = c(16, 9, 5, 21, 19)
set.2 = c(18, 3, 7)

# produce plots for these subsets: plot relationships between PCAs and expected orderings
pca.corrs.s1 = ggplot(long.df[!(long.df$Variable %in% c(set.1,set.2)),]) +
  geom_point(aes(x=Value,y=pca1), size=0.5, color="#BBBBBB") + 
  geom_smooth(aes(x=Value,y=pca1), method="lm", color="#444444") +
  geom_point(alpha=0.3, aes(x=Value,y=pca2-10), size=0.5, color="#FF888822") +
  geom_smooth(alpha=0.3, aes(x=Value,y=pca2-10), method="lm", color="#AA444422") +
  geom_point(alpha=0.3, aes(x=Value,y=pca3-20), size=0.5, color="#8888FF22") +
  geom_smooth(aes(x=Value,y=pca3-20), method="lm", color="#6666AA22") +
  geom_vline(xintercept=11) +
  facet_wrap(~name, nrow=2) + theme_minimal() + labs(x="Expected acquisition ordering", y="PCA projection (shifted)")
pca.corrs.s2 = ggplot(long.df[(long.df$Variable %in% c(set.1)),]) +
  geom_point(aes(x=Value,y=pca1), size=0.5, color="#BBBBBB22") + 
  geom_smooth(aes(x=Value,y=pca1), method="lm", color="#44444422") +
  geom_point(aes(x=Value,y=pca2-10), size=0.5, color="#FF8888") +
  geom_smooth(aes(x=Value,y=pca2-10), method="lm", color="#AA4444") +
  geom_point(aes(x=Value,y=pca3-20), size=0.5, color="#8888FF22") +
  geom_smooth(aes(x=Value,y=pca3-20), method="lm", color="#6666AA22") +
  geom_vline(xintercept=11) +
  facet_wrap(~name, nrow=1) + theme_minimal() + labs(x="Expected acquisition ordering", y="PCA projection (shifted)")
pca.corrs.s3 = ggplot(long.df[(long.df$Variable %in% c(set.2)),]) +
  geom_point(aes(x=Value,y=pca1), size=0.5, color="#BBBBBB22") + 
  geom_smooth(aes(x=Value,y=pca1), method="lm", color="#44444422") +
  geom_point(aes(x=Value,y=pca2-10), size=0.5, color="#FF888822") +
  geom_smooth(aes(x=Value,y=pca2-10), method="lm", color="#AA444422") +
  geom_point(aes(x=Value,y=pca3-20), size=0.5, color="#8888FF") +
  geom_smooth(aes(x=Value,y=pca3-20), method="lm", color="#6666AA") +
  geom_vline(xintercept=11) +
  facet_wrap(~name, nrow=1) + theme_minimal() + labs(x="Expected acquisition ordering", y="PCA projection (shifted)")

# pull together into a big figure
pca.corrs3.alt = ggarrange(pca.corrs.s1, 
                           ggarrange(pca.corrs.s2, pca.corrs.s3, nrow=1, labels=c("B", "C")), 
                           nrow=2, labels=c("A", ""))


long.df$Variable = as.numeric(long.df$Variable)

# scaling function for alternative, beeswarm PCA correlate plot
my.scale = function(x) {
  return(x)
  return((x-mean(x))/sd(x))
}

# set up dataframe for beeswarm plot
plot.long.df = long.df
plot.long.df$set = 1
plot.long.df$col[!(plot.long.df$Variable %in% c(set.1, set.2))] = my.scale(plot.long.df$pca1[!(plot.long.df$Variable %in% c(set.1, set.2))])
plot.long.df$set[plot.long.df$Variable %in% set.1] = 2
plot.long.df$col[plot.long.df$Variable %in% set.1] = my.scale(plot.long.df$pca2[plot.long.df$Variable %in% set.1])
plot.long.df$set[plot.long.df$Variable %in% set.2] = 3
plot.long.df$col[plot.long.df$Variable %in% set.2] = my.scale(plot.long.df$pca3[plot.long.df$Variable %in% set.2])

# beeswarm plot
pca.corrs.alt = ggplot(plot.long.df, aes(x=name, y=Value-22*(set-1), color=col)) + 
  geom_beeswarm(size=0.5,cex=0.5,priority="density") + scale_color_viridis_c(option="magma") + 
  geom_hline(yintercept = 11) + geom_hline(yintercept = 11-22) + geom_hline(yintercept = 11-44) +
  theme_minimal() +
  labs(x = "", y = "Expected acquistion timing (most covarying PCA)", color="PCA\nprojection") + theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_y_continuous(breaks = c(6, 11, 16, 6-22, 11-22, 16-22, 6-44, 11-44, 16-44), 
                     labels=c("earlier", "uniform", "later", 
                              "earlier", "uniform", "later",
                              "earlier", "uniform", "later")) +
  annotate("text", x = 1, y = 20, label = "PCA1") +
  annotate("text", x = 1, y = 20-22, label = "PCA2") +
  annotate("text", x = 1, y = 20-44, label = "PCA3") 

png("pca-corrs-alt.png", width=600*sf, height=400*sf, res=72*sf)
pca.corrs.alt
dev.off()

png("pca-corrs3-alt.png", width=600*sf, height=400*sf, res=72*sf)
pca.corrs3.alt
dev.off()

############### TASK 5 -- geographical regions and PCA stats

g.region.2 = ggplot(long.df[long.df$Variable == 1,], aes(x=region, y=pca2, fill=region, label=ccode)) +
  geom_boxplot(width=0.5) + geom_point() + geom_text_repel(size=2) + 
  scale_fill_viridis_d(option="magma") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none") 
g.region.3 = ggplot(long.df[long.df$Variable == 1,], aes(x=region, y=pca3, fill=region, label=ccode)) +
  geom_boxplot(width=0.5) + geom_point() + geom_text_repel(size=2) + 
  scale_fill_viridis_d(option="magma") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none") 
g.region.ind = ggplot(long.df[long.df$Variable %in% c(set.1,set.2),], 
                      aes(x=factor(name, levels=feature.names[c(set.1,set.2)+1]), y=Value, fill= region, label=ccode)) + 
  geom_boxplot() + 
  scale_fill_viridis_d(option="magma") +
  geom_hline(yintercept = 11) +
  geom_vline(xintercept = 5.5) +
  theme_minimal() + labs(x="Character", y="Expected acquisition ordering", fill="Region")

g.region.2z = g.region.2 + 
  geom_point(data = long.df[long.df$Variable == 1 & long.df$country == "Zanzibar",],
             aes(x=region, y=pca2, fill=region, label=ccode), pch=21, size=4)
g.region.3z = g.region.3 + 
  geom_point(data = long.df[long.df$Variable == 1 & long.df$country == "Zanzibar",],
             aes(x=region, y=pca3, fill=region, label=ccode), pch=21, size=4)
g.region.indz = g.region.ind + 
  geom_point(data = long.df[long.df$Variable %in% c(set.1, set.2) & long.df$country == "Zanzibar",],
             aes(x=factor(name, levels=feature.names[c(set.1,set.2)+1]), y=Value, fill=region, label=ccode), 
             pch=21, size=3, position=position_nudge(x = 0.33))

############### TASK 6 -- drug use data and acquisition

# load GLASS data
load("misc-data/glass_amc.Rdata")

# create drug use dataframe
d.use <- df.wide %>%
  group_by(CountryTerritoryArea) %>%
  summarise(across(-c(Year), median, na.rm = TRUE), .groups = "drop")

# pull drug codes
dcs = read.csv("misc-data/drug-codes.csv")
plot.dc = data.frame()
p.df = data.frame()
ref = 1

# relabel to produce readable annotations
row_strings <- paste0(dcs[,1], " (", dcs[,2], ")")
final_string <- paste(row_strings, collapse = "; ")
cat(gsub("-a", "", final_string))

# transformation function
myt = function(x) {
  return(sqrt(x))
}

# go through list of drugs, calculating linear model statistics and producing plot
for(i in 1:nrow(dcs)) {
  if(dcs$Relevant_ATC_Codes[i] %in% colnames(d.use)) {
    d.use.tmp = d.use[c("CountryTerritoryArea", dcs$Relevant_ATC_Codes[i])]
    colnames(d.use.tmp)[1] = "country"
    h.out.tmp = long.df[long.df$name == dcs$Resistance_Code[i], c("country", "region", "ccode", "name", "Value")]
    plot.dc.tmp = merge(d.use.tmp, h.out.tmp, by="country")
    colnames(plot.dc.tmp)[2] = "drug.use"
    my.lm = lm(plot.dc.tmp$Value ~ myt(plot.dc.tmp$drug.use))
    pval = round(summary(my.lm)$coefficients[2,4], digits=3)
    plot.dc.tmp$drug.label = paste0(dcs$Relevant_ATC_Codes[i], "/", dcs$Resistance_Code[i], ", p = ", pval, collapse="")
    plot.dc = rbind(plot.dc, plot.dc.tmp)
    p.df = rbind(p.df, data.frame(p = paste0("p=",pval,collapse=""), drug.label =paste0(dcs$Relevant_ATC_Codes[i], "/", dcs$Resistance_Code[i], ", p = ", pval, collapse=""), 
                                  x = mean(myt(plot.dc.tmp$drug.use)), y = 20))
  }
}

# pull all drug plots together
g.drugs = ggplot(plot.dc, aes(x=myt(drug.use), y=Value)) + 
  geom_smooth(method="lm", fill="#AAAAAA55", color="#AAAAAA99") +
  geom_point(size=1,aes(color=region, label=ccode)) + 
  geom_text_repel(size=2, alpha=0.7, aes(color=region, label=ccode)) + 
  theme_minimal() +
  facet_wrap(~drug.label, scales="free_x") +
  labs(x = "sqrt (median drug use)", y = "Expected\nacquisition ordering", color = "Region")

# reformat geographical covariate plot
g.covariates = ggarrange(
  ggarrange(
    g.region.2z + labs(x="Region", y="PCA2"), 
    g.region.3z + labs(x="Region", y="PCA3"), 
    nrow=1, labels=c("A", "B")), 
  g.region.indz + labs(y="Expected\nacquisition ordering"), 
  g.drugs + scale_color_viridis_d(option="magma") + theme(legend.position="none"), nrow=3, heights=c(0.7,0.7,1), labels=c("", "C", "D"))

# pull geographical and drug covariate plots together
png("covariates.png", width=500*sf, height=600*sf, res=72*sf)
g.covariates
dev.off()

############### TASK 7 -- summary global data 

if (!file.exists("clean/kleborate-dichotomized.csv")) {
  stop("Run preprocess_kleborate.R first!")
}

# read in resistance data
resistance.df <- read.csv("clean/kleborate-dichotomized.csv")
featurenames <- setdiff(colnames(resistance.df), "id")

# now working with global data, pull resistance feature info
resistance.df.l <- pivot_longer(resistance.df, setdiff(colnames(resistance.df), "id"))
resistance.df.l <- resistance.df.l[resistance.df.l$value==1,]

# create dataframe for UpSet plot
to.upset.plot = resistance.df.l %>% 
  group_by(id) %>% 
  dplyr::summarize(Genes = list(name)) %>%
  mutate(n_features = lengths(Genes))

# make UpSet plot
upset.plot = ggplot(to.upset.plot, aes(x = Genes, fill = n_features)) + 
  geom_bar() + 
  scale_x_upset(n_intersections=30) +
  scale_fill_gradientn(
    colours = c("blue", "orange", "red"),
    name = "Number of\nfeatures"
  ) +
  labs(x="Resistance character profiles", y="Number of genomes") +
  theme_minimal() + theme(legend.position = "none")

# create data frame for world map plot
c.df = data.frame()
cnames = names(country.list)
for(country in cnames) {
  tree.path <- paste0("clean/",country,".tsv")
  if (file.exists(tree.path)) {
    tmp = read.csv(tree.path, sep="\t")
    c.df = rbind(c.df, data.frame(country=country, count=nrow(tmp)))
  }
}

c.df$country = gsub("_", " ", c.df$country)
c.df$country[c.df$country=="United Kingdom"] <- "UK"

# get map data
world <- map_data("world")

# merge counts with map polygons
world_data <- world %>%
  left_join(c.df, by = c("region" = "country"))

# plot world map
plot.world = ggplot(world_data, aes(x = long, y = lat, group = group, fill = log10(count))) +
  geom_polygon(color = "gray70", size = 0.1) +
  scale_fill_viridis_c(option = "magma", na.value = "white", direction=-1) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(x=NULL, y=NULL, fill = "log10\n# genomes")

# initialise global data frame with an example (just for structure of data frame)
rdf = df[df$country=="Australia", c(1, 3, 5)]

# number of countries
ncountries = length(unique(df$country))

# loop over elements of the bubble plot, pulling average value over all countries for each
for(i in 1:nrow(rdf)) {
  refs = which(df$Time == rdf$Time[i] & df$OriginalIndex == rdf$OriginalIndex[i])
  rdf$Probability[i] = sum(df$Probability[refs])/ncountries
}

# produce global bubble plot
global.plot.alt = ggplot(data=rdf[rdf$Probability > 1/22,], aes(x=Time+1, y=feature.names[OriginalIndex+1], size=Probability, color=Time)) + 
  geom_point(shape=16) + theme_minimal() + theme(legend.position = "none") +
  scale_colour_gradientn(
    colours = c("darkblue", "blue", "red", "darkred"),
    values = scales::rescale(c(0, 7, 14, 21))
  ) + labs(y=NULL, x="Ordinal Time") 

# highlight bimodal features
global.plot.alt.hl = highlight_layer(global.plot.alt, c(2, 12, 18))

sf = 2
png("global-plot-alt.png", width=250*sf, height=250*sf, res=72*sf)
print(global.plot.alt.hl+xlab("Evolutionary ordering"))
dev.off()

# read in data on dates of samples and process
y.df = read.csv("misc-data/date-data.csv")
y2.df = y.df[y.df$isolation_year>2008,1:2]
y2.df = rbind(y2.df, data.frame(isolation_year = "Pre-2008", count=sum(y.df$count[y.df$isolation_year>0 & y.df$isolation_year<2008])))
y2.df = rbind(y2.df, data.frame(isolation_year = "NA", count=y.df$count[y.df$isolation_year==0 | y.df$isolation_year=="No value"]))

# produce plot
date.plot = ggplot(y2.df, aes(x=isolation_year, y=count)) + geom_col() +
  labs(x="Isolation year", y="Genomes") + theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) 

# combine world map, dates, and UpSet plots
g.fig.alt.1 = ggarrange(ggarrange(plot.world, 
                                  date.plot+ylab("\n\nGenomes"), labels=c("Ai", "ii"), nrow=2, heights=c(1.5,1)), 
                        upset.plot+ylab("\n\nNumber of\ngenomes")+
                          theme_combmatrix(combmatrix.panel.point.size = 1,
                                           combmatrix.panel.line.size  = 0.5) , nrow=1, labels=c("", "B"))

png("fig-alt-1.png", width=700*sf, height=300*sf, res=72*sf)
g.fig.alt.1
dev.off()

############### TASK 8 -- case study plot

# read in phylogeny for case study data
country = "Romania"
tree.path <- paste0("clean/",country,".nwk")
if (!file.exists(tree.path)) {
  stop(paste("No newick-tree for", country, "in clean directory"))
}

# produce curated tree for case study
ctree <- curate.tree(tree.path, 
                     "clean/kleborate-dichotomized.csv")

# curate and relabel data for simplicity
ctree.tmp = ctree
colnames(ctree.tmp$data) = gsub("_acquired", "-a", colnames(ctree.tmp$data))
colnames(ctree.tmp$data) = gsub("_mutations", "-m", colnames(ctree.tmp$data))
res.tmp = country.list$Romania$seed.3
res.tmp$bubbles$Name = gsub("_acquired", "-a", res.tmp$bubbles$Name)
res.tmp$bubbles$Name = gsub("_mutations", "-m", res.tmp$bubbles$Name)
res.tmp$featurenames = gsub("_acquired", "-a", res.tmp$featurenames)
res.tmp$featurenames = gsub("_mutations", "-m", res.tmp$featurenames)

g.2.alt = ggarrange(plotHypercube.curated.tree(ctree.tmp, font.size = 2.5, hjust=1) +
                      coord_cartesian(clip = "off") + theme(
                        plot.margin = unit(c(1, 1, 4, 1), "lines")  # top, right, bottom, left
                      ),
                    plotHypercube.sampledgraph2(res.tmp, truncate = 6, node.labels=FALSE, edge.label.size = 3,
                                                edge.check.overlap = FALSE, edge.label.angle = "none",
                                                no.times = TRUE),
                    b.plot("Romania") + theme(legend.position="right") + guides(color = "none") + xlab("Evolutionary ordering"),
                    # plotHypercube.bubbles(res.tmp, p.color = "#8888FF55") + labs(size="Probability"),
                    
                    nrow=1)
png("fig-g2alt.png", width=1100*sf*0.9, height=300*sf*0.7, res=72*sf)
g.2.alt
dev.off()

############### TASK 9 -- mega bubble plot

# get counts of genomes from each country
c.counts = unique(data.frame(name=world_data$region, count=world_data$count))
c.counts = c.counts[!is.na(c.counts$count),]
# use these to order country names
c.ordered = c.counts$name[order(c.counts$count)]
c.ordered = gsub("UK", "United_Kingdom", c.ordered)
c.ordered = gsub(" ", "_", c.ordered)
setdiff(df$country, c.ordered)

# edit names of features
df$new.name = gsub("_mutations", "-m", gsub("_acquired", "-a", df$Name))
# get mean ordering for each feature across dataset
df$mean = df$Time*df$Probability
f.names = unique(df$new.name)
f.means = c()
for(this.name in f.names) {
  f.means = c(f.means, sum(df$mean[df$new.name==this.name]))
}
# use this to order features
f.ordered = f.names[order(f.means)]

# make the megaplot
all.plot = ggplot(df[!(df$country %in% setdiff(df$country, c.ordered)),], 
       aes(x=Time, y=factor(country, levels=c.ordered), size=Probability, color = Time)) + 
  geom_point(shape =16) +
  facet_wrap(~factor(new.name, levels=f.ordered), nrow=1) +
  scale_color_gradientn(colours = c("darkblue", "blue", "red", "darkred")) +
  labs(x="Evolutionary ordering", y="Country") + guides(color = "none") + theme_minimal()

# output
sf = 3
png("all-plot.png", width=1700*sf, height=900*sf, res=72*sf)
print(all.plot)
dev.off()
