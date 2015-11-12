#' @title Score gene signatures based on PCA density
#'
#' @param data data frame of matrix with gene expression values where rows
#' represent genes and columns represent samples.
#' @param parsed logical.
#' @param pcs when set the function will compute the Euclidean norm for
#' selected principal components instead of using the first \code{n}
#' (overrides \code{varCutoff}). Vector of length 2.
#'
#' @return norms a vector with the Euclidean Norm of all the loadings.
#'
#' @export
scoreSigs <- function(data, parsed = FALSE, genes = (if (parsed) NULL),
                      signatures,  varCutoff = 0.75, denCutoffLow = 0.005,
                      denCutoffHigh = 0.02, pc = NULL, show.all = FALSE){

    if (!parsed)
        data <- readSamples(data, genes)

    scores <- do.call(rbind, mclapply(signatures, measureSpread, data,
                                      varCutoff, denCutoffLow, denCutoffHigh,
                                      pc, show.all))


    NAsigs <- scores[,1] == "NA"
    if (sum(NAsigs) > 0) {
        message(paste("Removed", sum(NAsigs), "signature(s)"))
        scores <- scores[ which(!NAsigs), ]
    }

    if (show.all) {
        scores <- as.data.frame(scores, stringsAsFactors = FALSE)
        colnames(scores) <- c("Density", "Size", "PCs")
        s <- sort(scores$Density, decreasing = TRUE, index.return = TRUE)
        scores <- scores[s$ix, ]
    }else {
        s <- sort(scores, decreasing = TRUE, index.return = TRUE)
        scores <- as.data.frame(scores[s$ix], stringsAsFactors = FALSE)
        colnames(scores) <- c("Density")
    }

    if (parsed)
        scores
    else
        list(data = data, scores = scores)
}
