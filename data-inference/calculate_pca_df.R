get.pca.df <- function() {
  if (!exists("country.list")) {
    name <- load("all_models.Rdata")
    country.list <- get(name)
  }
  countries <- names(country.list)
  
  # construct (a) a dataframe storing all bubble plots in long form
  # (b) a list of bubble plots for one seed per country in matrix form
  mlist = list()
  for(i in 1:length(country.list)) {
    b.1 = country.list[[i]]$seed.1$bubbles
    matrix_result <- reshape2::acast(b.1, Time ~ OriginalIndex, value.var = "Probability")
    mlist[[i]] = matrix_result
  }
  
  flattened_matrices <- lapply(mlist, as.vector)
  
  # Combine the vectors into a matrix where each row is a flattened matrix
  data_matrix <- do.call(rbind, flattened_matrices)
  # Perform PCA
  pca_result <- prcomp(data_matrix, scale. = TRUE)
  
  data.frame(country=countries, 
             pca1=pca_result$x[,1], 
             pca2=pca_result$x[,2], 
             pca3=pca_result$x[,3], 
             pca4=pca_result$x[,4], 
             pca5=pca_result$x[,5])
}
#ggplot2::ggplot(pca.df, ggplot2::aes(pca1, pca2, label=country)) + ggplot2::geom_text()
