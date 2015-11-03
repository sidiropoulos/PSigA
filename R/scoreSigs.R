#' @export
scoreSigs <- function(data, genes = NULL, signatures, parsed = FALSE, varCutoff = 0.75,
                      denCutoffLow = 0.005, denCutoffHigh = 0.02, pc = NULL, show.all = FALSE){

    if (!parsed)
        data <- readSamples(data, genes)

    scores <- do.call(rbind, mclapply(signatures, measureSpread, data,
                                      varCutoff, denCutoffLow, denCutoffHigh,
                                      pc, show.all))
    if (show.all)
        colnames(scores) <- c("Signature Size", "PCs", "Density")
    else
        colnames(scores) <- c("Density")

    scores <- as.data.frame(scores)

    #sort scores based on density
    s <- sort(scores$Density, decreasing = TRUE, index.return = TRUE)

    scores <- scores[s$ix, ]

    if (parsed)
        scores
    else
        list(data = data, scores = scores)
}
