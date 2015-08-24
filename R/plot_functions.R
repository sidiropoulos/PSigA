#' @title Produces a bibplot of the selected signature's PCA
#'
#' @description To do
#'
#' @author Nikos Sidiropoulos
#'
#' @param sigPCA a \code{\link{prcomp}} object, preferrably produced using \code{\link{signaturePCA}}
#' @param pcs which PCs to plot. default: 1,2
#' @param groups optional factor variable indicating the groups that the observations belong to.
#' If provided the points will be colored according to groups
#' @param db microarray platform
#' @param genenames when TRUE, gene symbols are used instead of probe ids as variable names
#' @param main plot title
#' @param palette colorbrewer palette scheme to be used. Applicable only when groups are provided.
#' @param ... methods passed to \code{\link{ggbiplot}}
#'
#' @return signature PCA biplot
#'
#' @import ggbiplot
#' @export

sigBiplot <- function(sigPCA, pcs = c(1,2), groups = NULL, genenames = TRUE, main = "Signature", db = "hgu133plus2.db",
                      palette = "Paired", ...) {

    if (genenames)
        rownames(sigPCA$rotation) <- probe2geneMap(rownames(sigPCA$rotation), db)

    ggbiplot(sigPCA, choices = pcs, groups = groups, ... ) + ggtitle(main) + scale_color_brewer(palette = palette)
}

#' @title Signature PCA plot
#'
#' @description To do..
#'
#' @author Nikos Sidiropoulos
#'
#' @param sigPCA a \code{\link{prcomp}} object, preferrably produced using \code{\link{signaturePCA}}
#' @param pcs which PCs to plot. default: 1,2
#' @param groups optional factor variable indicating the groups that the observations belong to.
#' @param text when TRUE it plots textual annotations instead of points based on \code{groups} parameter.
#' @param main plot title
#' @param palette colorbrewer palette scheme to be used.
#' @param ... methods passed to \code{\link{ggplot}}
#'
#' @import ggplot2
#' @export
#'
sigPlot <- function(sigPCA, pcs = c(1,2), groups = NULL, text = FALSE, main = "Signature",
                    palette = "Paired", ...)  {

    if (is.null(groups)){
        groups <- rep("sample", nrow(sigPCA$x))
    }

    d <- data.frame(groups, sigPCA$x[, pcs])

    xdata <- paste("PC", pcs[1], sep = "")
    ydata <- paste("PC", pcs[2], sep = "")

    p <- ggplot(d, aes(x = get(xdata), y = get(ydata), label = groups, colour = factor(groups)))
    p <- p + xlab(xdata) + ylab(ydata) + ggtitle(main)

    if (text)
        p <- p + geom_text(show_grid = FALSE, size = 4) + theme(legend.position="none")
    else
        p <- p + geom_point()

    p + scale_color_brewer(palette = palette)
}
