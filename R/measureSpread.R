#' @title Measures the area under the density curve on PC space
#'
#' @description \code{measureSpread} integrates the density curves of the first
#' \emph{n} principle components of a \code{\link{signaturePCA}}, where \emph{n}
#' is component where the cumulative variance is above a user specified cutoff.
#' Much more to be written...
#'
#' @author Nikos Sidiropoulos
#'
#' @param sigpca signature PCA. Object of class \code{\link{prcomp}}. Can be
#' produced with \code{\link{signaturePCA}}.
#' @param varCutoff cumulative variance cutoff. The function will measure the
#' area under the components that describe more or equal portion of the
#' variance defined by this parameter.
#' @param pc when set the function will compute the area in the first \code{pcs}
#' components instead of using the cumulative variance percentage.
#' (overrides \code{varCutoff}).
#' @param denCutoffLow lower density cutoff.
#' @param denCutoffHigh upper density cutoff.
#'
#' @return density area. If \code{show.all = TRUE} it returns a vector with the
#' number of genes used in the signature, the number of principal
#' components used and the density area.
#'
#' @import caTools
#' @export measureSpread

measureSpread <- function(signature, data, varCutoff = 0.75,
                          denCutoffLow = 0.005, denCutoffHigh = 0.05,
                          pc = NULL, show.all = FALSE)
{
    inSet <- signature %in% rownames(data)
    sigpca <- prcomp(t(data[signature[inSet], ]), center = TRUE,
                     scale = FALSE)

    #Get the proportion of the variance in the PC's
    propVar <- sigpca$sdev**2 / sum( sigpca$sdev**2 )

    if (!is.null(pc)) {
        idx <- pc
    } else {

        #Find for which component the cumulative proportion of the variance is
        #above the cutoff
        idx <- which(cumsum(sigpca$sdev**2) / sum(sigpca$sdev**2) > varCutoff)[1]
    }

    if (idx == 1){
        densities <- density(sigpca$x[, 1], bw = "SJ")
        densities <- list(densities)
    }else
        densities <- apply( sigpca$x[ , 1:idx ], 2, density, bw = "SJ" )

    areaSum <- 0

    for ( i in 1:idx )
    {
        curDensity <- densities[[ i ]]

        edges <- which( curDensity$y > denCutoffLow )
        edges <- c( min( edges ), max( edges ) )

        curDensity$y[curDensity$y > denCutoffHigh ] = denCutoffHigh
        curDensity$y[curDensity$y < denCutoffLow ] = 0

        area <- trapz( curDensity$x, curDensity$y )
        scaledArea <- area * propVar[ i ]

        areaSum <- areaSum + scaledArea
    }
    if (show.all)
        c(sum(inSet), idx, areaSum)
    else
        areaSum
}
