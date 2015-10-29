#' Import samples
#'
#'
#'
readSamples <- function(data, genes, method = c("max", "mean", "median")) {

    match.arg(method, c("max", "mean", "median"))
    print(method)
    nonNAgenes <- which(!is.na(genes))
    data <- data[nonNAgenes,]

    genes <- genes[nonNAgenes]

    a<- mclapply(unique(genes), .testFun, data, genes, method)

    d <- do.call(rbind.data.frame, a)
    rownames(d) <- unique(genes)
    d
}


.testFun <- function(gene, data, genes, method) {

    if (method == "max") {
        maxGene <- which.max(abs(apply(data[genes == gene, ], 1, mean)))
        data[genes == gene,][maxGene, ]
    } else
        apply(data[genes == gene,], 2, method)

}

