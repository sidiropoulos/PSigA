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
#' dummyData <- do.call(rbind, lapply(seq(0.1, 1, by = 0.1),
#'                      rnorm, n = 100, m = 6))
#' rownames(dummyData) <- paste(rep("gene", nrow(dummyData)),
#'                              seq(1, nrow(dummyData)), sep = "")
#' dummySig <- c("gene1", "gene8", "gene9", "gene10")
#'
#' dummyPca <- signaturePCA(dummyData, dummySig)
#'
#' @export signaturePCA
signaturePCA <- function(data, signature, center = TRUE, scale = FALSE, ...){

    inSet <- signature %in% rownames(data)
    if (sum(inSet) < 3){
        stop(strwrap("Less than 3 genes in the signature match the rows in the
                     data. Make sure your data and signature are in the same
                     gene format. Consider using parseData()", prefix = " ",
                     width = getOption("width")))
    }

    genes <- signature[signature %in% rownames(data)]
    prcomp(t(data[genes,]), center = center, scale = scale)
}
