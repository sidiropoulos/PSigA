#' @title Signature Biplot for principal compontents
#'
#' @description Produces a bibplot of the selected signature's PCA. Based on
#' \code{\link{ggbiplot}} and \code{\link{ggplot2}}.
#'
#' @param sigPCA a \code{\link{prcomp}} object, preferrably produced using
#' \code{\link{signaturePCA}}
#' @param groups optional factor variable indicating the groups that the
#' observations belong to. If provided the points will be colored accordingly
#' @param labels optional vector of labels for the observations
#' @param pcs which PCs to plot. default: 1,2
#' @param obs.size size of the points for the observations
#' @param var.size size of the text for the variables
#' @param var.scaled logical value. When set to TRUE the variable text size is
#' proportional to the Euclidean Norm of the variable's loading
#' (see \code{\link{measureLoadings}}).
#' @param main plot title
#' @param palette colorbrewer palette scheme to be used. Applicable only when
#' \code{groups} are provided.
#' @param ... methods passed to \code{\link{ggbiplot}}.
#'
#' @return A ggplot2 plot
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
#' sigBiplot(sigPCA)
#' sigBiplot(sigPCA, groups = factor(pData(vdx)$grade),
#'           main = "DOANE_BREAST_CANCER_ESR1_UP")
#'
#' @import ggbiplot
#' @export

sigBiplot <- function(sigPCA, groups = NULL, labels = NULL, pcs = c(1,2),
                      main = "", obs.size = 2, var.size = 3,
                      var.scaled = FALSE, palette = "Paired", ...){

    if (var.scaled) {
        norms <- measureLoadings(sigPCA)
        var.size <- norms*var.size
    }

    if (!is.null(labels)){
        p <- ggbiplot(sigPCA, choices = pcs, groups = groups, labels = labels,
                      labels.size = obs.size, varname.size = var.size, ... )
        p <- p + theme(legend.position="none")
        p <- p + scale_color_brewer(palette = palette)
    }else if (!is.null(groups)){
        p <- ggbiplot(sigPCA, choices = pcs, groups = groups,
                      varname.size = var.size, ... )
        p <- p + geom_point(aes(colour=groups), size = obs.size)
        p <- p + scale_color_brewer(palette = palette)
    }else {
        p <- ggbiplot(sigPCA, choices = pcs, varname.size = var.size, ...)
        p <- p + geom_point(size = obs.size)
    }
    p + ggtitle(main)
}

#' @title Signature PCA plot
#'
#' @description Plot of the selected signature's PCA using
#' \code{\link{ggplot2}}.
#'
#' @param sigPCA a \code{\link{prcomp}} object, preferrably produced using
#' \code{\link{signaturePCA}}
#' @param groups optional factor variable indicating the groups that the
#' observations belong to.
#' @param text when TRUE it plots textual annotations instead of points
#'  based on \code{groups} parameter.
#' @param pcs which PCs to plot. default: 1,2
#' @param main plot title
#' @param palette colorbrewer palette scheme to be used.
#' @param ... methods passed to \code{\link{ggplot}}
#'
#' @return A ggplot2 plot
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
#' sigPlot(sigPCA)
#' sigPlot(sigPCA, groups = factor(pData(vdx)$grade),
#'           main = "DOANE_BREAST_CANCER_ESR1_UP")
#'
#' @import ggplot2
#' @export
#'
sigPlot <- function(sigPCA, groups = NULL, text = FALSE, pcs = c(1,2),
                    main = "", palette = "Paired", ...)  {

    # Groups flag
    gFlag <- TRUE

    if (is.null(groups)){
        groups <- rep("sample", nrow(sigPCA$x))
        gFlag <- FALSE
    }

    data <- data.frame(groups, sigPCA$x[, pcs])

    colnames(data) <- c("groups", "x", "y")
    xlab <- paste("PC", pcs[1], sep = "")
    ylab <- paste("PC", pcs[2], sep = "")

    #Prevent R CMD check "no global definition.." NOTE
    x <- y <- NULL

    if (gFlag) {
        p <- ggplot(data, aes(x = x, y = y, label = groups,
                              colour = factor(groups)))
        p <- p + scale_color_brewer(palette = palette) + labs(colour='Groups')
    } else {
        p <- ggplot(data, aes(x = x, y = y))
    }

    p <- p + xlab(xlab) + ylab(ylab) + ggtitle(main)

    if (text){
        p <- p + geom_text(show_grid = FALSE, size = 4)
        p <- p + theme(legend.position="none")
    }else {
        p <- p + geom_point()
    }
    p
}
