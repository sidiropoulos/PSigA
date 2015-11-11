#' @title Produces a bibplot of the selected signature's PCA
#'
#' @description To do
#'
#' @param sigPCA a \code{\link{prcomp}} object, preferrably produced using
#' \code{\link{signaturePCA}}
#' @param pcs which PCs to plot. default: 1,2
#' @param groups optional factor variable indicating the groups that the
#' observations belong to. If provided the points will be colored accordingly
#' @param labels optional vector of labels for the observations
#' @param obs.size size of the points for the observations
#' @param var.size size of the text for the variables
#' @param var.scaled logical value. When set to TRUE the variable text size is
#' proportional to the Euclidean Norm of the variable's loading
#' (see \code{\link{measureLoadings}}).
#' @param main plot title
#' @param palette colorbrewer palette scheme to be used. Applicable only when
#' \code{groups} are provided.
#' @param ... methods passed to \code{\link{ggbiplot}}
#'
#' @return signature PCA biplot
#'
#' @import ggbiplot
#' @export

sigBiplot <- function(sigPCA, pcs = c(1,2), groups = NULL, labels = NULL,
                      main = "Signature", obs.size = 2, var.size = 3,
                      var.scaled = FALSE, palette = "Paired", ...) {

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
#' @description To do..
#'
#' @param sigPCA a \code{\link{prcomp}} object, preferrably produced using
#' \code{\link{signaturePCA}}
#' @param pcs which PCs to plot. default: 1,2
#' @param groups optional factor variable indicating the groups that the
#' observations belong to.
#' @param text when TRUE it plots textual annotations instead of points
#'  based on \code{groups} parameter.
#' @param main plot title
#' @param palette colorbrewer palette scheme to be used.
#' @param ... methods passed to \code{\link{ggplot}}
#'
#' @import ggplot2
#' @export
#'
sigPlot <- function(sigPCA, pcs = c(1,2), groups = NULL, text = FALSE,
                    main = "Signature", palette = "Paired", ...)  {

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
