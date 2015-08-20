#' @title Measures the area under the density curve on PC space
#'
#' @description \code{measureSpread} integrates the density curves of the first \emph{n} principle components of a
#' \code{\link{signaturePCA}}, where \emph{n} is component where the cumulative variance is above a user specified cutoff.
#' Much more to be written...
#'
#' @param sigpca signature PCA. Object of class \code{\link{prcomp}}. Can be produced with \code{\link{signaturePCA}}.
#' @param varCutoff cumulative variance cutoff. The function will measure the area under the components that describe more
#' or equal portion of the variance defined by this parameter.
#' @param denCutoffLow lower density cutoff.
#' @param denCutoffHigh upper density cutoff.
#' @param sigName ?
#'
#' @return a data frame with columns: Signature Name, Number of principal components used, Density area.
#'
#' @import caTools
#' @export measureSpread

measureSpread <- function( sigpca, varCutoff = 0.75, denCutoffLow = 0.005, denCutoffHigh = 0.05,
                           sigName = "NA" )
{

    #Get the proportion of the variance in the PCA's
    propVar <- sigpca$sdev**2 / sum( sigpca$sdev**2 )

    #Find for which component the cumulative proportion of the variance is above the cutoff
    idx <- which( cumsum( sigpca$sdev**2 )  / sum( sigpca$sdev**2 ) > varCutoff )[ 1 ]

    densities <- apply( sigpca$x[ , 1:idx ], 2, density, bw = "SJ" )

    areaProd <- 1
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
    c( sigName, idx, areaSum )
}
