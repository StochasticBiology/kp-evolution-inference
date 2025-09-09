#' Plot a word cloud of a feature
#' 
#' @return a list of plots
plot.feature.wordcloud <- function() {
    library(ggplot2)
    library(wordcloud)
    #library(wordcloud2)
    source("~/doc/grafisk-profil-uib.R")
    library(ggpubr)
    df <- read.csv("../klebsiella-datasets/raw/Klebsiella pneumoniae__kleborate.csv")
    cols <- colnames(df)
    cols <- cols[which(endsWith(cols, "_acquired"))]
    df <- df[,cols]
    colnames(df) <- gsub("_acquired", "", colnames(df))


    make.df <- function(col) {
      s <- paste(unlist(df[,col]), collapse=';')
      
      s <- gsub("(-;|;-)","", s)
      s <- unlist(strsplit(s, ";"))
      
      # Remove " +13V" like syntax (unknown significance)
      s <- gsub(" \\+[0-9][0-9][0-9]?[A-Z]$", "", s)
      
      # Remove asterix, hats and questionmarks all denoting distance to known alleles
      s <- gsub("[*^?]*$", "", s)
      
      # Remove .v1, .v2 etc
      s <- gsub("\\.v[1-9]$", "", s)
      
      words <- unique(s)
      freqs <- sapply(words, \(x) sum(x == s))
      
      data.frame(word = words, freq = freqs)
    }

    plot.wordcloud <- function(col) {
      feature.df <- make.df(col)
      svg(paste0("plots/wordclouds/",col,"-wordcloud.svg"))
      wordcloud(feature.df$word, feature.df$freq,
                colors = rep_len(uib.blue[1:6], length.out=nrow(feature.df)),
                ordered.colors=TRUE)
      dev.off()
    }

    plot.lollipop <- function(col) {
      feature.df <- make.df(col)
      feature.df$word <- factor(feature.df$word)
      feature.df$word <- reorder(feature.df$word, feature.df$freq, sum)
      ggplot(feature.df) + 
        geom_point(aes(freq, word)) +
        geom_segment(aes(x=0, xend=freq, y=word, yend=word))
      ggsave(paste0("plots/frequencies/",col,"-dots.png"), 
             width = max(length(feature.df$word) / 8, 5), 
             height = max(length(feature.df$word) / 3,4))
    }

    sapply(colnames(df), plot.wordcloud)
    sapply(colnames(df), plot.lollipop)
}
