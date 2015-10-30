#' Import samples
#'
#' @param data data frame of matrix with gene expression values where rows
#' represent genes and columns represent samples.
#' @param genes vector of length \code{nrow(data)} containing gene names in
#' "HGNC" format.
#'
#' @return kati gamato
#'
#' @export
readSamples <- function(data, genes) {

    nonNAgenes <- which(!is.na(genes))
    data <- data[nonNAgenes,]

    genes <- genes[nonNAgenes]

    a<- mclapply(unique(genes), .maxMedian, data, genes)

    d <- do.call(rbind, a)

    rownames(d) <- unique(genes)
    d
}

.maxMedian <- function(gene, data, genes) {

    data <- data[genes == gene, ]
    maxGene <- which.max(abs(apply(data, 1, median)))
    data[maxGene, ]

}

