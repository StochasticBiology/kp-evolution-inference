path <- "data-inference"  # replace with your folder path
files <- list.files(path, pattern = "\\.R$", full.names = TRUE)
sapply(files, source)

library(tidyverse)
library(ggrepel)

if (!exists("country.list")) {
  name <- load("all_models.Rdata")
  country.list <- get(name)
}

igj.get.pca.df = function() {
  if (!exists("country.list")) {
    name <- load("all_models.Rdata")
    country.list <- get(name)
  }
  countries <- names(country.list)
  
  # construct (a) a dataframe storing all bubble plots in long form
  # (b) a list of bubble plots for one seed per country in matrix form
  mlist = list()
  for(i in 1:100) {
#  for(i in 1:length(country.list)) {
    b.1 = country.list[[i]]$seed.1$bubbles
    b.1$Probability = b.1$Probability - 1/20
  #  b.1$Probability[b.1$Probability < 0.1] = runif(1,min=0, max=0.001)
    matrix_result <- reshape2::acast(b.1, Time ~ OriginalIndex, value.var = "Probability")
    mlist[[i]] = matrix_result
  }
  
  flattened_matrices <- lapply(mlist, as.vector)
  
  # Combine the vectors into a matrix where each row is a flattened matrix
  data_matrix <- do.call(rbind, flattened_matrices)
  # Perform PCA
  pca_result <- prcomp(data_matrix, scale. = TRUE)
  
  data.frame(country=countries[1:100], 
             pca1=pca_result$x[,1], 
             pca2=pca_result$x[,2], 
             pca3=pca_result$x[,3], 
             pca4=pca_result$x[,4], 
             pca5=pca_result$x[,5])
}

region <- read.csv("misc-data/data-Ib27t.csv") # gbd superregion
seed.1 <- list()
summary.stats <- data.frame()
all.bubbles = data.frame()
country.names = names(country.list)
i = 0

for (country in country.list) {
  i = i+1
  bubbles <- country$seed.1$bubbles
  seed.1 <- c(seed.1, list(bubbles))
  bubbles %>%
    group_by(Name) %>%
    summarise(expected.order = sum(Time * Probability)) %>%
    pivot_wider(names_from=Name, values_from=expected.order) -> obs
  summary.stats <- rbind(summary.stats, obs)
  bubbles$country = country.names[[i]]
  bubbles$region = ""
  ref = which(region$country == country.names[[i]])
  if(length(ref) > 0) {
    bubbles$region = region$val[ref]
  }
  all.bubbles = rbind(all.bubbles, bubbles)
}
bubbles <- country.list[["USA"]]$seed.2$bubbles
seed.1 <- c(seed.1, list(bubbles))
bubbles %>%
  group_by(Name) %>%
  summarise(expected.order = sum(Time * Probability)) %>%
  pivot_wider(names_from=Name, values_from=expected.order) -> obs
summary.stats <- rbind(summary.stats, obs)
bubbles$country = "USA2"
bubbles$region = ""
ref = which(region$country == "USA")
if(length(ref) > 0) {
  bubbles$region = region$val[ref]
}
all.bubbles = rbind(all.bubbles, bubbles)
bubbles <- country.list[["USA"]]$seed.3$bubbles
seed.1 <- c(seed.1, list(bubbles))
bubbles %>%
  group_by(Name) %>%
  summarise(expected.order = sum(Time * Probability)) %>%
  pivot_wider(names_from=Name, values_from=expected.order) -> obs
summary.stats <- rbind(summary.stats, obs)
bubbles$country = "USA3"
bubbles$region = ""
ref = which(region$country == "USA")
if(length(ref) > 0) {
  bubbles$region = region$val[ref]
}
all.bubbles = rbind(all.bubbles, bubbles)

load("new-data/zanzibar.Rdata")
bubbles <- model.zanzibar.1$bubbles
seed.1 <- c(seed.1, list(bubbles))
bubbles %>%
  group_by(Name) %>%
  summarise(expected.order = sum(Time * Probability)) %>%
  pivot_wider(names_from=Name, values_from=expected.order) -> obs
summary.stats <- rbind(summary.stats, obs)
bubbles$country = "Zanzibar"
bubbles$region = ""
ref = which(region$country == "Zanzibar")
if(length(ref) > 0) {
  bubbles$region = region$val[ref]
}
all.bubbles = rbind(all.bubbles, bubbles)

write.table(all.bubbles, "igj-all-bubbles.csv", sep=";", row.names=FALSE, quote=FALSE)
ggplot(all.bubbles, aes(x=Time, y=Probability, group=country, color=factor(region))) +
  geom_line() + facet_wrap(~OriginalIndex) + scale_y_log10() + geom_hline(aes(yintercept=1/20)) +
  theme(legend.position = "none") +
  scale_color_viridis_d(option="magma")

ggplot(all.bubbles[all.bubbles$Probability > 0.1,], aes(x = Time, y = OriginalIndex, color=region)) + 
  geom_jitter() +
  theme(legend.position="none")

length(unique(all.bubbles$country[all.bubbles$Probability > 0.1]))
length(unique(all.bubbles$country[all.bubbles$Probability > 0.0]))

names(seed.1) <- names(country.list)

pca.df <- igj.get.pca.df()

summary.stats = cbind(pca.df, summary.stats[1:100,])

load("misc-data/glass_amc.Rdata") # loads a data frame called df.wide
region <- read.csv("misc-data/data-Ib27t.csv") # gbd superregion

gbd.map = region %>%
  inner_join(summary.stats) %>%
  select(country)

summary.stats$region = ""
for(i in 1:nrow(summary.stats)) {
  ref = which(region$country == summary.stats$country[i])
  if(length(ref) > 0) {
    summary.stats$region[i] = region$val[ref]
  }
}

ggplot(summary.stats, aes(pca1,pca2,label=country,color=region)) + geom_point() + geom_text_repel()
  
############

gbd.map$country.gbd <- gbd.map$country
gbd.map$country.kp <- gbd.map$country
gbd.map <- gbd.map[,c("country.gbd", "country.kp")]

gbd.map <- rbind(gbd.map, c("Burkina Faso","Burkina_Faso"))
gbd.map <- rbind(gbd.map, c("Czech Republic","Czech_Republic"))
gbd.map <- rbind(gbd.map, c("Arab Republic of Egypt","Egypt"))
gbd.map <- rbind(gbd.map, c("The Gambia","Gambia"))
gbd.map <- rbind(gbd.map, c(NA,"Guadeloupe"))
gbd.map <- rbind(gbd.map, c(NA,"Hong_Kong"))
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

setdiff(summary.stats$country, region$country)
region.corrected <- full_join(region, gbd.map, by=join_by(country==country.gbd))
region.corrected$val[which(region.corrected$country.kp == "Guadeloupe")] = "Latin America and Caribbean"
region.corrected$val[which(region.corrected$country.kp == "Hong_Kong")] = "High-income"
setdiff(region.corrected$country, gbd.map$country.gbd)

summary.stats %>%
  left_join(region.corrected, by = join_by(country == country.kp)) -> summary.region

region.corrected$country
region.corrected$country.glass <- ifelse(region.corrected$country %in% df.wide$CountryTerritoryArea,
                                         region.corrected$country, NA)

region.corrected$country.glass[which(region.corrected$country == "United Kingdom")] = "United Kingdom of Great Britain and Northern Ireland"
region.corrected$country.glass[which(region.corrected$country == "Tanzania")] = "United Republic of Tanzania"
region.corrected$country.glass[which(region.corrected$country == "Côte d'Ivoire")] = "C<U+00F4>te d'Ivoire"
region.corrected$country.glass[which(region.corrected$country == "Czech Republic")] = "Czechia"
region.corrected$country.glass[which(region.corrected$country == "Arab Republic of Egypt")] = "Egypt"
region.corrected$country.glass[which(region.corrected$country == "Islamic Republic of Iran")] = "Iran (Islamic Republic of)"

region.selection <- region.corrected %>%
  select(val, country.kp, country.glass) %>%
  na.omit()

df.wide %>%
  inner_join(region.selection, by=join_by(CountryTerritoryArea == country.glass)) %>%
  left_join(summary.stats, by=join_by(country.kp == country))-> glass.gbd

summary.stats %>%
  left_join(region.corrected, by=join_by(country==country.kp)) %>%
  ggplot(aes(pca1,pca2,label=country,color=val)) + 
  geom_point() + geom_text_repel() + theme_light() -> pca.1.2
