#' @title Measure variable loadings
#' @description Measures the Euclidean Norm of the PCA loadings (rotation
#' matrix) of the first \code{n} principal components, where the cumulative
#' variance on the \code{n}-th principal component is above a cutoff.
#'
#' @param pca \code{\link{prcomp}} object.
#' @param varCutoff cumulative variance cutoff.
#' @param pc when set the function will compute the Euclidean norm using the
#' first \code{pc} principal components instead of using the first \code{n}
#' (overrides \code{varCutoff}).
#'
#' @return norms a vector with the Euclidean Norm of all the loadings.
#'
#' @examples
#'
#' data(iris)
#'
#' pca <- prcomp(iris[, -5], scale = FALSE)
#' measureLoadings(pca)
#'
#' #use only the first 2 principal components
#' measureLoadings(pca, pc = 2)
#'
#' @export

measureLoadings <- function(pca, varCutoff = 0.75, pc = NULL) {

    norms <- c()
    for (i in 1:nrow(pca$rotation)){

        if (is.null(pc))
            pc <- which(summary(pca)$importance[3,] > varCutoff)[1]

        norms <- c(norms, norm(as.matrix(pca$rotation[i, 1:pc])))
    }

    names(norms) <- rownames(pca$rotation)
    norms
}
