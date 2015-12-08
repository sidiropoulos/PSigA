peakDistance <- function(signature, data, pc = 2, adjust = 3/4,
                           show.all = FALSE, scale = FALSE, denCutoffLow = 0.005)
{
    inSet <- signature %in% rownames(data)

    if (sum(inSet) < 3)
        return("NA")

    sigpca <- prcomp(t(data[signature[inSet], ]), center = TRUE,
                     scale = scale)

    for (i in 1:pc) {
        densities <- apply( sigpca$x[ , 1:pc ], 2, density, bw = "SJ",
                            adjust = adjust)
    }

    score <- 0

    for ( i in 1:pc )
    {
        curDensity <- densities[[ i ]]
        curDensity$y[curDensity$y < denCutoffLow ] <- 0

        tp <- which(diff(sign(diff(d$y)))==-2) + 1
        if (length(tp) < 2)
            next

        maxPeaks <- sort(curDensity$y[tp], decreasing = TRUE,
                       index.return = TRUE)$ix[1:2]
        range <- range(curDensity$x[tp][maxPeaks])
        # range <- range(curDensity$x[tp$tppos])

        score <- score + abs(range[2] - range[1])

        #consider scaling like below
        #scaledArea <- area * propVar[ i ]
    }

    #normalize for signature length
    score <- score / log(sum(inSet))

    if (show.all)
        c(score, sum(inSet))
    else
        score
}
