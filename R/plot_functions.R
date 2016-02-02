#' @title Signature Biplot for principal compontents
#'
#' @description Produces a bibplot of the selected signature's PCA. Based on
#' \code{\link{ggbiplot}} and \code{ggplot2}.
#'
#' @param data gene expression matrix where rownames correspond to unique gene
#' identifiers in \code{signature} format, and columns correspond to samples.
#' @param signature character vector containing the signature's gene
#' identifiers.
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
#' dummyData <- do.call(rbind, lapply(seq(0.1, 0.9, by = 0.1),
#'                      rnorm, n = 100, m = 6))
#' #add a row with bimodal gene expression
#' dummyData <- rbind(dummyData, c(rnorm(70, 6, 0.1), rnorm(30, 9, 0.1)))
#'
#' rownames(dummyData) <- paste(rep("gene", nrow(dummyData)),
#'                              seq(1, nrow(dummyData)), sep = "")
#'
#' dummySig <- c("gene1", "gene8", "gene9", "gene10")
#'
#' sigBiplot(dummyData, dummySig)
#'
#' groups <- factor(c(rep("GroupA", 70), rep("GroupB", 30)))
#' sigBiplot(dummyData, dummySig, groups, main = "dummySig")
#'
#' @import ggbiplot
#' @export
sigBiplot <- function(data, signature, groups = NULL, labels = NULL,
                      pcs = c(1,2), main = "", obs.size = 2, var.size = 3,
                      var.scaled = FALSE, palette = "Paired", ...){

    pca <- signaturePCA(data, signature)

    if (var.scaled) {
        norms <- .loadingsNorm(pca)
        var.size <- norms*var.size
    }

    if (!is.null(labels)){
        p <- ggbiplot(pca, choices = pcs, groups = groups, labels = labels,
                      labels.size = obs.size, varname.size = var.size, ... )
        p <- p + theme(legend.position="none")
        p <- p + scale_color_brewer(palette = palette)
    }else if (!is.null(groups)){
        p <- ggbiplot(pca, choices = pcs, groups = groups,
                      varname.size = var.size, ... )
        p <- p + geom_point(aes(colour=groups), size = obs.size)
        p <- p + scale_color_brewer(palette = palette)
    }else {
        p <- ggbiplot(pca, choices = pcs, varname.size = var.size, ...)
        p <- p + geom_point(size = obs.size)
    }
    p + ggtitle(main)
}

#' @title Signature PCA plot
#'
#' @description Plot of the selected signature's PCA using
#' \code{ggplot2}.
#'
#' @param data gene expression matrix where rownames correspond to unique gene
#' identifiers in \code{signature} format, and columns correspond to samples.
#' @param signature character vector containing the signature's gene
#' identifiers.
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
#' dummyData <- do.call(rbind, lapply(seq(0.1, 0.9, by = 0.1),
#'                      rnorm, n = 100, m = 6))
#'
#' #add a row with bimodal gene expression
#' dummyData <- rbind(dummyData, c(rnorm(70, 6, 0.1), rnorm(30, 9, 0.1)))
#'
#' rownames(dummyData) <- paste(rep("gene", nrow(dummyData)),
#'                              seq(1, nrow(dummyData)), sep = "")
#'
#' dummySig <- c("gene1", "gene8", "gene9", "gene10")
#'
#' sigPlot(dummyData, dummySig)
#'
#' groups <- factor(c(rep("GroupA", 70), rep("GroupB", 30)))
#' sigPlot(dummyData, dummySig, groups, main = "dummySig")
#'
#' @import ggplot2
#' @export
#'
sigPlot <- function(data, signature, groups = NULL, text = FALSE, pcs = c(1,2),
                    main = "", palette = "Paired", ...)  {

    pca <- signaturePCA(data, signature)

    # Groups flag
    gFlag <- TRUE

    if (is.null(groups)){
        groups <- rep("sample", nrow(pca$x))
        gFlag <- FALSE
    }

    data <- data.frame(groups, pca$x[, pcs])

    colnames(data) <- c("groups", "x", "y")
    xlab <- paste("PC", pcs[1], sep = "")
    ylab <- paste("PC", pcs[2], sep = "")

    #Prevent R CMD check "no global definition.." NOTE
    x <- y <- NULL

    if (gFlag) {
        p <- ggplot(data, aes(x = x, y = y, label = groups,
                              colour = factor(groups)))
        p <- p + scale_color_brewer(palette = palette) + labs(colour = "Groups")
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

.loadingsNorm <- function(pca, varCutoff = 0.75, pc = 2) {

    norms <- c()
    for (i in 1:nrow(pca$rotation)){

        if (is.null(pc))
            pc <- which(summary(pca)$importance[3,] > varCutoff)[1]

        norms <- c(norms, norm(as.matrix(pca$rotation[i, 1:pc])))
    }

    names(norms) <- rownames(pca$rotation)
    norms
}
