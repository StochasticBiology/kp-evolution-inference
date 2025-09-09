library(ggplot2)

get.seeds <- function(fitted.cube.country) {
  list(seed.1 = fitted.cube.country[1:12],
       seed.2 = fitted.cube.country[13:24],
       seed.3 = fitted.cube.country[25:36])
}

full.length.featurename <- function(short.featurename) {
  switch(short.featurename,
         "AGly" = "Aminoglycosides",
         "Sul" = "Sulfoamides",
         "Tmt" = "Trimethoprim",
         "Bla_Carb" = "Bla_Carb",
         "Phe" = "Phenicols",
         "Bla" = "Bla",
         "Flq" = "Fluoroquinolones",
         "Bla_ESBL" = "Bla_ESBL",
         "Tet" = "Tetracycline",
         "MLS" = "Macrolides",
         "Rif" = "Rifampicin",
         "Gly" = "Glycopeptides",
         "Col" = "Colistin",
         "Fcyn" = "Fosfomycin",
         "Bla_inhR" = "Bla_inhR",
         "Bla_ESBL_inhR" = "Bla_ESBL_inhR",
         "Tgc" = "Tigecycline",
         short.featurename
         )
}

#' Plot bubbles for a list of countries grouped by feature
#' 
#' @param country.list A list of models with three seeds labeled seed.1, seed.2
#'  and seed.3
#' @return a ggplot
plot.bubbles.by.feature <- function(country.list, transition_set = NA) {
    features = unique(country.list[[1]]$seed.1$bubbles$Name)

    ft.df <- data.frame(Time=numeric(),
                        ReorderedIndex=numeric(),
                        OriginalIndex=numeric(),
                        Name=character(),
                        Probability=numeric(),
                        Seed=numeric(),
                        Country=character())

    for (ft in features) {
      for (country in names(country.list)) {
        seeds <- country.list[[country]]
        
        bubbles.1 <- seeds$seed.1$bubbles[which(seeds$seed.1$bubbles$Name == ft),]
        bubbles.2 <- seeds$seed.2$bubbles[which(seeds$seed.2$bubbles$Name == ft),]
        bubbles.3 <- seeds$seed.3$bubbles[which(seeds$seed.3$bubbles$Name == ft),]
        
        bubbles.1$Seed <- rep(1, length(bubbles.1$Time))
        bubbles.2$Seed <- rep(2, length(bubbles.1$Time))
        bubbles.3$Seed <- rep(3, length(bubbles.1$Time))
        
        bubbles = rbind(bubbles.1, bubbles.2, bubbles.3)
        
        bubbles$Country = rep(country, length(bubbles$Time))
        ft.df <- rbind(ft.df, bubbles)
      }  
    }

    ft.df$Seed <- factor(ft.df$Seed)
    
    if (anyNA(transition_set)) {
      ft.df$Country <- factor(ft.df$Country, levels=sort(unique(ft.df$Country), decreasing=TRUE))
    } else {
      ft.df$Country <- factor(ft.df$Country, levels=names(sort(transition_set)))
    }
    weighted_time = ft.df$Probability * (ft.df$Time + 1)
    weight <- \(x) mean(weighted_time[which(ft.df$Name == x)]) * 17
    id <- order(sapply(features,weight))
    sorted.features <- features[id]
    sorted.levels <- sapply(features, full.length.featurename)[id]

    ft.df$Feature.name <- factor(sapply(ft.df$Name, full.length.featurename),
                                 levels=sorted.levels)

    library(dplyr)
    plt <- ggplot(ft.df) + 
      geom_point(aes(Time, Country, 
                     size=Probability, 
                     color=Seed,
                     alpha=0.3)) +
      xlab("Ordered time") +
      ggtitle("") +
      facet_wrap(ft.df$Feature.name, nrow = 1)
    ggsave("plots/feature-plot.png", plot=plt, width=10, height=8)

    for (feature in unique(ft.df$Feature.name)) {
      next
      df <- ft.df[which(ft.df$Feature.name == feature),]
      ggplot(df) + 
        geom_point(aes(Time, Country, 
                       size=Probability, 
                       color=Seed,
                       alpha=0.3)) +
        xlab("Ordered time") +
        ggtitle(feature) + 
        theme_classic()
      ggsave(paste0("plots/bubbles_by_feature/",feature,".svg"), width=4, height=4)
    }
    plt
}
