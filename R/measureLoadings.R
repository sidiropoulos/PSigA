#' @title Measure variable loadings
#' @description Measures the Euclidean Norm of the PCA loadings (rotation matrix) of the first
#' \code{n} principal components, where the cumulative variance on the \code{n}-th principal component
#' is above a cutoff.
#'
#' @author Nikos Sidiropoulos
#'
#' @param pca \code{\link{prcomp}} object.
#' @param varCutoff cumulative variance cutoff.
#'
#' @return norms a vector with the Euclidean Norm of all the loadings.
#'
#' @export

measureLoadings <- function(pca, varCutoff = 0.75) {
    pc <- which(summary(pca)$importance[3,] > varCutoff)[1]
    norms <- c()
    for (i in 1:nrow(pca$rotation)){
        norms <- c(norms, norm(as.matrix(pca$rotation[i, 1:pc])))
    }
    norms
}
