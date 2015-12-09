#' @title Peak distance in a 2D signature PCA
#'
#' @description Locates and measures the Euclidean distance between the two
#' highest peaks of the 2D density function, estimated from the samples'
#' distribution in the first two principal components in a signature PCA
#' (see \code{\link{signaturePCA}}).
#'
#' @param signature character vector with the signature's gene identifiers in
#' HGNC format.
#' @param data Gene expression matrix where rownames correspond to unique gene
#' identifiers (HGNC format \link{http://www.genenames.org}) and columns
#' correspond to samples.
#' @param threshold density cutoff. Density values lower than the
#' \code{threshold} will not be considered peaks. Usefull when outliers are
#' present in the PCA.
#' @param n Number of grid points in each direction. Can be scalar or a
#' length-2 integer vector. See \link{kde2d}.
#'
#' @return A vector of length 2. The first value corresponds to the score and
#' the second to the number of genes in the \code{signature} that were found in
#' the \code{data}.
#'
#' @import MASS
#' @export
peakDistance2d <- function(signature, data, threshold, n = 200){

    #
    inSet <- signature %in% rownames(data)

    if (sum(inSet) < 3)
        return(c(-1,0))

    sigpca <- prcomp(t(data[signature[inSet], ]), center = TRUE,
                     scale = FALSE)

    d <- kde2d(sigpca$x[,1], sigpca$x[,2], n = n)
    peaks <- .find2Dpeaks(d, threshold)

    if (length(peaks) <= 1)
        return(0)
    if (nrow(peaks) < 2)
        return(0)

    maxPeaks <- peaks[sort(d$z[peaks], decreasing = TRUE,
                           index.return = TRUE)$ix[1:2], ]

    dist <- norm(matrix(c(d$x[maxPeaks[,1]],d$y[maxPeaks[,2]]), nrow = 2))

    c(dist, sum(inSet))
}

.find2Dpeaks <- function(d, threshold = 0.015){

    peaks <- c()
    for (i in 2:(nrow(d$z)-1)){
        for (j in 2:(ncol(d$z)-1)){

            if (d$z[i,j] < threshold)
                next

            if (sum(d$z[i,j] >= matrix(d$z[(i-1):(i+1), (j-1):(j+1)])) == 9)
                peaks <- rbind(peaks, c(i,j))
        }
    }

    peaks
}
