#' Cluster cubes by their inferred acquisition order
#' 
#' @param country.list A nested list of model runs with one entry for each 
#'  country containing a sublist labeled seed.1, seed.2 and seed.3
#' @return ggplot with multidimentional scaling for the cubes
cluster.cubes <- function(model.list, n_transition=c()) {
    average.bubbles <- function(model.output) {
      sorted <- \(bub) bub[order(bub$OriginalIndex, bub$Name),]$Probability
      unlist(lapply(t(as.matrix(sorted(model.output$seed.1$bubbles), 
                                sorted(model.output$seed.2$bubbles), 
                                sorted(model.output$seed.3$bubbles))), 
                    mean))
    }
    
    L1.norm.difference.matrix <- function(cube.list) {
      cube.params <- lapply(cube.list, average.bubbles)
      diffs <- numeric()
      for (c1 in cube.params) {
        for (c2 in cube.params) {
          diffs[[length(diffs) + 1]] <- sum(abs(c1 - c2))
        }
      }
      
      dim(diffs) <- c(length(cube.list), length(cube.list))
      diffs
    }

    plot.MDS <- function(cube.list) {
      mat <- L1.norm.difference.matrix(cube.list)
      colnames(mat) <- names(cube.list)       
      rownames(mat) <- names(cube.list)
      mds <- cmdscale(mat)
      colnames(mds) <- c("MDS1", "MDS2")
      
      FAIL.ggplot <- require(ggplot2)
      FAIL.ggrepel <- require(ggrepel)
      if (!(FAIL.ggrepel) || !(FAIL.ggplot)) {
        print("Missing either ggrepel or ggplot2 for MDS ploting")
        return (-1)
      }
      if (length(n_transition) == 0) {
        ggplot(mds, aes(x = MDS1, y = MDS2)) + 
          geom_point() +
          geom_text_repel(aes(label=rownames(mds)),
                          box.padding = 0.1,
                          force_pull = 3,
                          force = 0.5) +
          theme_classic()
      } else {
        ggplot(mds, aes(x = MDS1, y = MDS2, color = log(n_transition))) + 
          geom_point() +
          geom_text_repel(aes(label=rownames(mds)),
                          box.padding = 0.1,
                          force_pull = 3,
                          force = 0.5) +
          theme_classic()
      }
    }
    
    plot.MDS(model.list) + 
      ggtitle("Classical (Metric) Multidimensional Scaling of Bubbles")
}
