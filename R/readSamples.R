#' Import samples
#'
#' @param data data frame of matrix with gene expression values where rows
#' represent genes and columns represent samples.
#' @param genes vector of length \code{nrow(data)} containing gene names in
#' "HGNC" format.
#'
#' @return A data frame ...
#'
#' @export
readSamples <- function(data, genes) {

    #find and remove NA values
    nonNAgenes <- which(!is.na(genes))
    data <- data[nonNAgenes,]
    genes <- genes[nonNAgenes]

    #select the probe/isoform with the highest median expression
    d <- do.call(rbind, mclapply(unique(genes), .maxMedian, data, genes))

    rownames(d) <- unique(genes)
    d
}

.maxMedian <- function(targetGene, data, genes) {

    match <- genes == targetGene
    if (sum(match) == 1){
        data[match,]
    }else {
        maxGene <- which.max(abs(apply(data[match,], 1, median)))
        data[maxGene, ]
    }

}
