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
