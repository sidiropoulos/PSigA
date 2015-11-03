#' @export
scoreSigs <- function(data, genes = NULL, signatures, parsed = FALSE, varCutoff = 0.75,
                      denCutoffLow = 0.005, denCutoffHigh = 0.02, pc = NULL, show.all = FALSE){

    if (!parsed)
        data <- readSamples(data, genes)

    r <- mclapply(a, measureSpread, data, varCutoff, denCutoffLow, denCutoffHigh, pc, show.all)
    r <- do.call(rbind, r)

    if (show.all)
        colnames(r) <- c("Signature Size", "PCs", "Density area")
    else
        colnames(r) <- c("Density area")

    if (parsed)
        as.data.frame(r)
    else
        list(data = data, scores = r)
}
