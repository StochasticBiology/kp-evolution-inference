library(ggplot2)
library(ggbeeswarm)
library(countrycode)
library(hypertrapsct)
library(tidyr)
library(dplyr)

df = read.csv("igj-all-bubbles.csv", sep=";", stringsAsFactors = FALSE)
region <- read.csv("raw/data-Ib27t.csv") # gbd superregion
region$oldval = region$val
region$val = gsub("Central Europe, Eastern Europe, and Central Asia", "CEEECA", region$val)
region$val = gsub("High-income", "HI", region$val)
region$val = gsub("Latin America and Caribbean", "LAC", region$val)
region$val = gsub("North Africa and Middle East", "NAME", region$val)
region$val = gsub("South Asia", "SA", region$val)
region$val = gsub("Southeast Asia, East Asia, and Oceania", "SEEAO", region$val)
region$val = gsub("Sub-Saharan Africa", "SSAf", region$val)

ggplot(df, aes(x=Time, y=Probability, group=country, color=factor(region))) +
  geom_line() + facet_wrap(~OriginalIndex) + scale_y_log10() + geom_hline(aes(yintercept=1/20)) +
  theme(legend.position = "none") +
  scale_color_viridis_d(option="magma")

tmp = df[1:22,]
feature.names = tmp$Name[order(tmp$OriginalIndex)]
feature.names = gsub("_acquired", "-a", feature.names)
feature.names = gsub("_mutations", "-m", feature.names)

countries <- unique(df$country)

wide_df <- df %>%
  mutate(ProbTime = Probability * Time) %>%
  group_by(country, OriginalIndex) %>%
  summarise(SumProbTime = sum(ProbTime, na.rm = TRUE), 
            .groups = "drop") %>%
  pivot_wider(names_from = OriginalIndex, values_from = SumProbTime)

wide_df2 <- df %>%
  mutate(Prob2Time = Probability * Time * Time, ProbTime = Probability * Time) %>%
  group_by(country, OriginalIndex) %>%
  summarise(SumProbTime = sum(ProbTime, na.rm = TRUE), 
            SD = sqrt(sum(Prob2Time, na.rm = TRUE) - sum(ProbTime, na.rm = TRUE)**2), 
            .groups = "drop") %>%
  pivot_wider(names_from = OriginalIndex, values_from = c(SumProbTime, SD))

gbd.map = data.frame()
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

region = rbind(region, data.frame(country="Hong Kong", val="Southeast Asia, East Asia, and Oceania"))
region = rbind(region, data.frame(country="Guadeloupe", val="Latin America and Caribbean"))

for(i in 1:nrow(region)) {
  ref = which(gbd.map[,1] == region$country[i])
  if(length(ref) > 0) {
    region$country[i] = gbd.map[ref,2]
  }
}

wide_df$region = ""
for(i in 1:nrow(wide_df)) {
  ref = which(region$country == wide_df$country[i])
  if(length(ref) > 0) {
    wide_df$region[i] = region$val[ref]
  }
}

wide_df$ccode = countrycode(wide_df$country, origin="country.name", destination = "iso3c")

# Convert to matrix (optional)
probtime_matrix <- as.matrix(wide_df[,-c(1, ncol(wide_df))])  # Remove 'country' column
rownames(probtime_matrix) <- wide_df$country

pca_result <- prcomp(probtime_matrix, scale. = TRUE)
pve <- (pca_result$sdev)^2 / sum(pca_result$sdev^2) * 100

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

g.pca = ggplot(pca.df, aes(x=pca1, y=pca2, color=region, label=country, fill=region)) + 
  geom_point() + geom_text_repel(size=2.5, max.overlaps = 20, alpha=0.5)   + labs(x="PCA1", y="PCA2")

#g.pca = ggplot(pca.df, aes(x=pca1, y=pca2, color=region, label=ccode, fill=region)) + 
#  geom_point() + geom_text_repel(max.overlaps = 40) 

g.pca.ellipse = g.pca + stat_ellipse(
  geom = "path",
  linewidth=2,
  alpha = 0.25      # transparency
  #color = NA        # optional: remove ellipse borders
)  +  theme_minimal()+ theme(legend.position = "bottom")

g.pca2 = ggplot(pca.df, aes(x=pca3, y=pca2, color=region, label=country, fill=region)) + 
  geom_point() + geom_text_repel(size=2.5, max.overlaps = 20, alpha=0.5) + theme_minimal()  + labs(x="PCA3", y="PCA2")

g.pca2.ellipse = g.pca2 + stat_ellipse(
  geom = "path",
  linewidth=2,
  alpha = 0.25      # transparency
 # color = NA        # optional: remove ellipse borders
) + theme(legend.position = "bottom")

pca.df$country == wide_df$country

all.df = cbind(wide_df, pca.df[,2:ncol(pca.df)])


######

bubble.plot = function(cnames) {
  ggplot(df[df$country %in% cnames,], aes(x=Time, y=Name, size=Probability, color=country)) + 
    geom_point(alpha=0.5) +
  theme_minimal()
}

bubble.plot(c("Venezuela", "Nigeria"))
bubble.plot(c("Senegal", "Argentina"))

g.pca.ex.1 = plotHypercube.bubbles.compare(list(country.list$Venezuela$seed.1, country.list$Nigeria$seed.1 ), 
                                           expt.names = c("Venezuela", "Nigeria"), sqrt.trans=TRUE,
                                           fill.name = "Country") +
  scale_y_continuous(breaks = 0:21, labels=feature.names) + 
  scale_fill_manual(values=c("#FF000088", "#0000FF88"))

plotHypercube.bubbles.compare(list(country.list$Cameroon$seed.1, country.list$Kenya$seed.1 ), sqrt.trans=TRUE)
g.pca.ex.2 = plotHypercube.bubbles.compare(list(country.list$Senegal$seed.1, country.list$Argentina$seed.1 ), 
                                           expt.names = c("Senegal", "Argentina"), sqrt.trans=TRUE, 
                                           fill.name = "Country") +
  scale_y_continuous(breaks = 0:21, labels=feature.names) + 
  scale_fill_manual(values=c("#FF000088", "#0000FF88"))

plotHypercube.bubbles.compare(list(country.list$Austria$seed.1, country.list$Myanmar$seed.1 ), sqrt.trans=TRUE)


my.ann = function(y) {
  return( annotate("rect",
                   xmin = -1, xmax = 22,     # x-range of the box
                   ymin = y-0.5, ymax = y+0.5,     # y-range of the box
                   color = "black",
                   fill=NA,
                   alpha = 0.5             # transparency
  ))
}

g.pca.part.1 = ggarrange(
  ggarrange(g.pca.ellipse, g.pca2.ellipse, nrow=2, labels=c("A", "B")),
  ggarrange(g.pca.ex.1, g.pca.ex.2  + my.ann(21) + my.ann(16) + my.ann(9) +my.ann(19) , nrow=2, labels=c("C", "D")),
  nrow=1
)

g.pca.part.1

ggplot(df[df$country=="Venezuela",], aes(x=Time, y=Name, size=Probability)) + geom_point()
ggplot(df[df$country=="Venezuela",], aes(x=Time, y=Name, size=Probability)) + geom_point()

long_df <- all.df %>%
  pivot_longer(
    cols = -c("country", "ccode", "region", "pca1", "pca2", "pca3", "pca4", "pca5"),                     # all columns except 'country'
    names_to = "Variable",         # name for new column holding former column names
    values_to = "Value"           # name for new column holding the values
  )
long.df = as.data.frame(long_df)
long.df$Variable = as.numeric(long.df$Variable)
long.df$name = feature.names[long.df$Variable+1]
ggarrange(
ggplot(long.df, aes(x=name, y=Value, color=pca1)) + 
  geom_beeswarm() + scale_color_viridis_b(option="magma") +
  theme(axis.text.x = element_text(angle=45,hjust=1)),
ggplot(long.df, aes(x=name, y=Value, color=pca2)) + 
  geom_beeswarm() + scale_color_viridis_b(option="magma") +
  theme(axis.text.x = element_text(angle=45,hjust=1)),
ggplot(long.df, aes(x=name, y=Value, color=pca3)) + 
  geom_beeswarm() + scale_color_viridis_b(option="magma") +
  theme(axis.text.x = element_text(angle=45,hjust=1)),
nrow=3
)

ggplot(long.df, aes(x = Value, y=pca1)) + geom_point() + facet_wrap(~name)
ggplot(long.df, aes(x = Value, y=pca2)) + geom_point() + facet_wrap(~name)
ggplot(long.df, aes(x = Value, y=pca3)) + geom_point() + facet_wrap(~name)
ggplot(long.df, aes(x = Value, y=pca4)) + geom_point() + facet_wrap(~name)
ggplot(long.df, aes(x = abs(Value-11), y=pca1)) + geom_point() + facet_wrap(~Variable)

pca.corrs = ggplot(long.df) +
  geom_point(aes(x=Value,y=pca1), size=0.5, color="grey") + 
  geom_smooth(aes(x=Value,y=pca1), method="lm", color="#444444") +
  geom_point(aes(x=Value,y=pca2-10), size=0.5, color="#FF8888") +
  geom_smooth(aes(x=Value,y=pca2-10), method="lm", color="#AA4444") +
  geom_point(aes(x=Value,y=pca3-20), size=0.5, color="#8888FF") +
  geom_smooth(aes(x=Value,y=pca3-20), method="lm", color="#6666AA") +
  geom_vline(xintercept=11) +
  facet_wrap(~name) + theme_minimal() + labs(x="Expected acquisition ordering", y="PCA projection (shifted)")
  
long.sub = long.df[long.df$Variable %in% c(16, 18, 20, 3, 5, 6, 7, 9),]

set.1 = c(16, 9, 5, 21, 19)
set.2 = c(18, 3, 7)

long.df$Variable = as.numeric(long.df$Variable)
g.pca.part.2 = ggarrange(
 ggplot(long.df[!(long.df$Variable %in% c(set.1, set.2)), ], aes(x=name, y=Value, color=pca1)) + 
   #ggplot(long.df, aes(x=name, y=Value, color=pca1)) + 
    geom_beeswarm(size=1) + scale_color_viridis_b(option="magma") + 
    geom_hline(yintercept = 11) +
    scale_x_discrete(limits=feature.names) + 
    theme_minimal() +
    labs(x = "", y = "PCA1") +
    theme(axis.text.x = element_text(angle=45,hjust=1), legend.position="none"),
  ggplot(long.df[long.df$Variable %in% set.1,], aes(x=name, y=Value, color=pca2)) + 
    geom_beeswarm(size=1) + scale_color_viridis_b(option="magma") + 
    geom_hline(yintercept = 11) +
    theme_minimal() +
    labs(x = "", y = "PCA2") +
    scale_x_discrete(limits=feature.names) + theme(axis.text.x = element_text(angle=45,hjust=1), legend.position="none"),
  ggplot(long.df[long.df$Variable %in% set.2,], aes(x=name, y=Value, color=pca3)) + 
    geom_beeswarm(size=1) + scale_color_viridis_b(option="magma") + 
    geom_hline(yintercept = 11) +
    theme_minimal() +
    labs(x = "", y = "PCA3") +
    scale_x_discrete(limits=feature.names) + theme(axis.text.x = element_text(angle=45,hjust=1), legend.position="none"),
  nrow=3, heights=c(1.5,1,1)
)

my.scale = function(x) {
  return(x)
  return((x-mean(x))/sd(x))
}
plot.long.df = long.df
plot.long.df$set = 1
plot.long.df$col[!(plot.long.df$Variable %in% c(set.1, set.2))] = my.scale(plot.long.df$pca1[!(plot.long.df$Variable %in% c(set.1, set.2))])
plot.long.df$set[plot.long.df$Variable %in% set.1] = 2
plot.long.df$col[plot.long.df$Variable %in% set.1] = my.scale(plot.long.df$pca2[plot.long.df$Variable %in% set.1])
plot.long.df$set[plot.long.df$Variable %in% set.2] = 3
plot.long.df$col[plot.long.df$Variable %in% set.2] = my.scale(plot.long.df$pca3[plot.long.df$Variable %in% set.2])

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



sf = 2
png("pca-igj.png", width=800*sf, height=600*sf, res=72*sf)
g.pca.part.1
dev.off()

png("pca-ids.png", width=500*sf, height=500*sf, res=72*sf)
g.pca.part.2
dev.off()

png("pca-corrs.png", width=600*sf, height=400*sf, res=72*sf)
pca.corrs
dev.off()

png("pca-corrs-alt.png", width=600*sf, height=400*sf, res=72*sf)
pca.corrs.alt
dev.off()

ggplot(long.sub, aes(x = Value, y=pca1)) + geom_point() + facet_wrap(~Variable)
ggplot(long.sub, aes(x = Value, y=pca2)) + geom_point() + facet_wrap(~Variable) # 16, 9, 5
ggplot(long.sub, aes(x = Value, y=pca3)) + geom_point() + facet_wrap(~Variable) # 18, 3, 7

ggarrange(
ggplot(long.df[long.df$Variable == 1,], aes(x=region, y=pca1, color=region, label=country)) +
  geom_boxplot() + geom_text_repel(),
ggplot(long.df[long.df$Variable == 1,], aes(x=region, y=pca2, color=region, label=country)) +
  geom_boxplot() + geom_text_repel(),
ggplot(long.df[long.df$Variable == 1,], aes(x=region, y=pca3, color=region, label=country)) +
  geom_boxplot() + geom_text_repel(),
nrow=3
)

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
  theme_minimal() + labs(x="Feature", y="Expected acquisition ordering", fill="Region")

png("pca-orders.png", width=500*sf, height=400*sf, res=72*sf)
ggarrange(ggarrange(g.region.2, g.region.3, nrow=1, labels=c("A", "B")), g.region.ind, nrow=2, labels=c("", "C"))
#ggarrange(g.region.2, g.region.3, g.region.ind, nrow=3)
dev.off()


#########

load("data/glass_amc.Rdata")

d.use <- df.wide %>%
  group_by(CountryTerritoryArea) %>%
  summarise(across(-c(Year), median, na.rm = TRUE), .groups = "drop")

dcs = read.csv("drug-codes.csv")
plot.dc = data.frame()
p.df = data.frame()
ref = 1

myt = function(x) {
  return(sqrt(x))
}

for(i in 1:nrow(dcs)) {
  if(dcs$Relevant_ATC_Codes[i] %in% colnames(d.use)) {
  d.use.tmp = d.use[c("CountryTerritoryArea", dcs$Relevant_ATC_Codes[i])]
  colnames(d.use.tmp)[1] = "country"
  h.out.tmp = long.df[long.df$name == dcs$Resistance_Code[i], c("country", "region", "ccode", "name", "Value")]
  plot.dc.tmp = merge(d.use.tmp, h.out.tmp, by="country")
  colnames(plot.dc.tmp)[2] = "drug.use"
  plot.dc.tmp$drug.label = paste0(dcs$Relevant_ATC_Codes[i], "/", dcs$Resistance_Code[i], collapse="")
  my.lm = lm(plot.dc.tmp$Value ~ myt(plot.dc.tmp$drug.use))
  pval = round(summary(my.lm)$coefficients[2,4], digits=3)
  plot.dc = rbind(plot.dc, plot.dc.tmp)
  p.df = rbind(p.df, data.frame(p = paste0("p=",pval,collapse=""), drug.label =paste0(dcs$Relevant_ATC_Codes[i], "/", dcs$Resistance_Code[i], collapse=""), 
                                x = mean(myt(plot.dc.tmp$drug.use)), y = 20))
  }
}

  ggplot(plot.dc, aes(x=myt(drug.use), y=Value)) + 
    geom_smooth(method="lm", fill="#AAAAFF55", color="#AAAAFF99") +
    geom_point(aes(color=region, label=ccode)) + 
    geom_text_repel(size=2.5, alpha=0.7, aes(color=region, label=ccode)) + 
    geom_text(data = p.df, aes(x = x, y = y, label = p),
              inherit.aes = FALSE) +
    theme_minimal() +
    facet_wrap(~drug.label, scales="free_x") +
    labs(x = "sqrt (median drug use)", y = "Expected acquisition ordering", color = "Region")
  
  ggplot(plot.dc[plot.dc$drug.label=="J01DH/Bla_Carb-a",], aes(x=myt(drug.use), y=Value)) + 
    geom_smooth(method="lm", fill="#AAAAFF55", color="#AAAAFF99") +
    geom_point(aes(color=region)) + 
    geom_text_repel(size=3, alpha=0.7, aes(color=region, label=country)) + 
    geom_text(data = p.df[p.df$drug.label=="J01DH/Bla_Carb-a",], aes(x = x, y = y, label = p),
              inherit.aes = FALSE) +
    theme_minimal() +
    facet_wrap(~drug.label, scales="free_x") +
    labs(x = "sqrt (median drug use)", y = "Expected acquisition ordering", color = "Region")
  
  library(lme4)
  
  my.lmm = lmer(Value ~ myt(drug.use) + (myt(drug.use) | drug.label), data = plot.dc)
  my.null = lmer(Value ~ 1 + (1 | drug.label), data = plot.dc)
  
  anova(my.lmm, my.null)
  
 ##########
  
  
  plot.dc.tmp = plot.dc[[ref]] 
#  dc.list[[ref]] = ggplot(plot.dc[[ref]], aes(x=plot.dc[[ref]][[2]], y=1/Value, label=ccode)) +
#                        geom_point() + geom_smooth(method="lm") +
#    labs(x=colnames(plot.dc[[ref]])[2], y=paste0("1/",plot.dc[[ref]][1,5],collapse=""))
  my.lm = lm(plot.dc.tmp$Value ~ sqrt(plot.dc.tmp[[2]]))
  pval = round(summary(my.lm)$coefficients[2,4], digits=3)
  dc.list[[ref]] = ggplot(plot.dc.tmp, aes(x=sqrt(plot.dc.tmp[[2]]), y=(Value), label=ccode)) +
                           geom_point() + geom_smooth(method="lm") +
        labs(x=paste0("sqrt ", colnames(plot.dc.tmp)[2], collapse=""), y=paste0("",plot.dc.tmp[1,5],collapse="")) +
    annotate("text", x=0.15, y = max(plot.dc.tmp$Value)-1, label=paste0(i," p=",pval,collapse=""))
  ref = ref+1
  }
}
ggarrange(plotlist=dc.list)

