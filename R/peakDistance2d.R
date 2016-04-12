#' @title Peak distance in a 2D signature-PCA
#'
#' @description Measures the Euclidean distance between the two
#' highest peaks of the 2D density function. The density is estimated from the
#' samples' distribution in the first two principal components of a PCA, using
#' the genes of a given signature.
#'
#' @param signature character vector with the signature's gene identifiers
#' @param data Gene expression matrix where rownames correspond to unique gene
#' identifiers in \code{signature} fortmat and columns correspond to samples.
#' @param threshold density cutoff. Density values lower than the
#' \code{threshold} will not be considered peaks. Useful when outliers are
#' present in the PCA.
#' @param n Number of grid points in each direction. Can be scalar or a
#' length-2 integer vector. See \link{kde2d}.
#' @param magnitude When TRUE the score is multiplied by the cluster density.
#' Default: FALSE.
#' @param scale a logical value indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place.
#' See prcomp for more details. Default: FALSE.
#' @param filtered logical value indicating if genes in the supplied
#' \code{signature} list that are not present in the \code{data} have been
#' filtered out. Default: FALSE.
#'
#' @return A vector of length 2. The first value corresponds to the score and
#' the second to the number of genes in the \code{signature} that were found in
#' the \code{data}.
#'
#' @examples
#'
#' dummyData <- do.call(rbind, lapply(seq(0.1, 0.9, by = 0.1),
#'                      rnorm, n = 100, m = 6))
#'
#' #add a row with bimodal gene expression
#' dummyData <- rbind(dummyData, c(rnorm(70, 6, 0.1), rnorm(30, 9, 0.1)))
#'
#' rownames(dummyData) <- paste(rep("gene", nrow(dummyData)),
#'                              seq(1, nrow(dummyData)), sep = "")
#' rownames(dummyData)
#'
#' dummySig <- c("gene1", "gene8", "gene9", "gene10", "gene20", "gene30")
#'
#' peakDistance2d(dummySig, dummyData)
#'
#' #values correspond to the peak distance, the number of genes from the
#' #signature found in the data and the total number of genes in the signature
#' #respectively
#'
#' #removing the bimodal gene from the signature results to a lower score
#' peakDistance2d(dummySig[-4], dummyData)
#'
#' @import MASS
#' @export
peakDistance2d <- function(signature, data, threshold = 0.005, n = 200,
                           magnitude = FALSE, scale = FALSE, filtered = FALSE){

    if (!filtered)
        signature <- .filterSignatures(signature, rownames(data))

    #if less than 3 genes are found return -1
    if (length(signature) < 3) {
        stop(strwrap("Too few genes (< 3) from the selected signature are
                     present in the dataset."))
    }

    sigpca <- prcomp(t(data[signature, ]), center = TRUE,
                     scale = scale)

    d <- kde2d(sigpca$x[,1], sigpca$x[,2], n = n)
    peaks <- .find2Dpeaks(d, threshold)

    #if less than 2 peaks are found, return 0
    if (length(peaks) <= 1)
        return(c(0, length(signature)))

    if (nrow(peaks) < 2)
        return(c(0, length(signature)))

    maxPeaks <- peaks[sort(d$z[peaks], decreasing = TRUE,
                           index.return = TRUE)$ix[1:2], ]

    dist <- norm(matrix(c(d$x[maxPeaks[,1]],d$y[maxPeaks[,2]]), nrow = 2))

    if (magnitude)
        c(dist*d$z[maxPeaks][1]*d$z[maxPeaks][2], length(signature))
    else
        c(dist, length(signature))
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
