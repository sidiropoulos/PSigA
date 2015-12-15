#' @title PCA based on a given gene-signature
#'
#' @description \code{signaturePCA} performs principal component analysis on
#' the given data matrix and gene-signature and returns the results as an object
#' of class \code{prcomp}.
#'
#' @param signature character vector with the signature's gene identifiers.
#' @param data Gene expression matrix where rownames correspond to unique gene
#' identifiers in \code{signature} format and columns correspond to samples.
#' @param center a logical value indicating whether the variables should be
#' shifted to be zero centered. See \code{\link{prcomp}} for more details.
#' @param scale a logical value indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place. See
#' \code{prcomp} for more details.
#' @param ... arguments passed to \code{prcomp}.
#'
#' @return A list with class \code{prcomp}.
#'
#' @examples
#'
#' require(Biobase)
#' require(breastCancerVDX)
#' data(vdx)
#' data(MSigDB)
#'
#' #get the first 5000 probes of the vdx array
#' VDX <- readSamples(data = exprs(vdx)[1:5000,],
#'                    genes = fData(vdx)$Gene.symbol[1:5000])
#'
#' sigPCA <- signaturePCA(MSigDB[["DOANE_BREAST_CANCER_ESR1_UP"]], VDX)
#'
#' #plot using sigPlot
#' sigPlot(sigPCA)
#'
#' @export signaturePCA
signaturePCA <- function(signature, data, center = TRUE, scale = FALSE, ...){

    genes <- signature[signature %in% rownames(data)]
    prcomp(t(data[genes,]), center = center, scale = scale, ...)
}
