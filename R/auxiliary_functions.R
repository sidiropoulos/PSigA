#' @title PCA based on a given gene-signature
#'
#' @description \code{signaturePCA} performs principal component analysis on
#' the given data matrix and gene-signature and returns the results as an object
#' of class \code{prcomp}.
#'
#' @param signature character vector with the signature's gene identifiers
#' @param data Gene expression matrix where rownames correspond to unique gene
#' identifiers (HGNC format \link{http://www.genenames.org}) and columns
#' correspond to samples.
#' @param center a logical value indicating whether the variables should be
#' shifted to be zero centered. See \code{\link{prcomp}} for more details.
#' @param scale a logical value indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place. See
#' \code{prcomp} for more details.
#' @param ... arguments passed to \code{prcomp}.
#'
#' @export signaturePCA
signaturePCA <- function(signature, data, center = TRUE, scale = FALSE, ...){

    genes <- signature[signature %in% rownames(data)]
    prcomp(t(data[genes,]), center = center, scale = scale, ...)
}
