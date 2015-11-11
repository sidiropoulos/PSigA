#' @title Measure variable loadings
#' @description Measures the Euclidean Norm of the PCA loadings (rotation
#' matrix) of the first \code{n} principal components, where the cumulative
#' variance on the \code{n}-th principal component is above a cutoff.
#'
#' @param pca \code{\link{prcomp}} object.
#' @param varCutoff cumulative variance cutoff.
#' @param pcs when set the function will compute the Euclidean norm for
#' selected principal components instead of using the first \code{n}
#' (overrides \code{varCutoff}). Vector of length 2.
#'
#' @return norms a vector with the Euclidean Norm of all the loadings.
#'
#' @export

measureLoadings <- function(pca, varCutoff = 0.75, pcs = NULL) {
    pc <- which(summary(pca)$importance[3,] > varCutoff)[1]
    norms <- c()
    for (i in 1:nrow(pca$rotation)){

        if (!is.null(pcs))
            norms <- c(norms, norm(as.matrix(pca$rotation[i, pcs[1]:pcs[2]])))
        else
            norms <- c(norms, norm(as.matrix(pca$rotation[i, 1:pc])))
    }

    norms
}
