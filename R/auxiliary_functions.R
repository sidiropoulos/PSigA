#' @title PCA based on a given gene-signature
#'
#' @description \code{signaturePCA} performs principal component analysis on
#' the given data matrix and gene-signature and returns the results as an object
#' of class \code{prcomp}.
#'
#' @param signature character vector with the signature's gene identifiers.
#' @param data Gene expression matrix where rownames correspond to unique gene
#' identifiers in \code{signature} format and columns correspond to samples.
#' @param center a logical value indicating whether the variables should be
#' shifted to be zero centered. See \code{\link{prcomp}} for more details.
#' @param scale a logical value indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place. See
#' \code{prcomp} for more details.
#' @param ... arguments passed to \code{prcomp}.
#'
#' @return A list with class \code{prcomp}.
#'
#' @examples
#'
#' dummyData <- do.call(rbind, lapply(seq(0.1, 1, by = 0.1),
#'                      rnorm, n = 100, m = 6))
#' rownames(dummyData) <- paste(rep("gene", nrow(dummyData)),
#'                              seq(1, nrow(dummyData)), sep = "")
#' dummySig <- c("gene1", "gene8", "gene9", "gene10")
#'
#' dummyPca <- signaturePCA(dummyData, dummySig)
#'
#' @export signaturePCA
signaturePCA <- function(data, signature, center = TRUE, scale = FALSE, ...){

    inSet <- signature %in% rownames(data)
    if (sum(inSet) < 3){
        stop(strwrap("Less than 3 genes in the signature match the rows in the
                     data. Make sure your data and signature are in the same
                     gene format. Consider using parseData()", prefix = " ",
                     width = getOption("width")))
    }

    genes <- signature[signature %in% rownames(data)]
    prcomp(t(data[genes,]), center = center, scale = scale)
}

#' @title PCA-based Gene Set Enrichment Analysis (GSEA)
#'
#' @description ...
#'
#' @param data Gene expression matrix where rownames correspond to unique gene
#' identifiers in \code{signature} format and columns correspond to samples.
#' @param signatures list where every entry is a character vector of
#' the gene names that correspond to a gene-pathway. If the gene format is
#' different from the one in \code{rownames(data)}, use \code{\link{parseData}}
#' first.
#' @param p.adj p-value correction method. See \code{\link{p.adjust}}.
#' @param pcs principal components to perform GSEA on. Default: c(1,2).
#' @param filtered logical value indicating if genes in the supplied
#' \code{signature} list that are not present in the \code{data} have been
#' filtered out. Default: FALSE.
#' @export
rankedGSEA <- function(data, signatures, p.adj = p.adjust.methods,
                       pcs = c(1,2), filtered = FALSE) {

    if (!filtered) {
        signatures <- lapply(signatures, .filterSignatures, rownames(data))
    }



    ### Check for signature length!

    p.adj = match.arg(p.adj)

    pca <- prcomp(t(data))

    gsea <- list()

    for (i in pcs){

        loadings <- pca$rotation[ ,i]
        ranking <- rank(abs(loadings))

        tmp <- do.call(rbind, mclapply(signatures, .rankTest, loadings,
                                       ranking))
        tmp <- as.data.frame(tmp)
        colnames(tmp) <- c("ks.statistic", "p.value")

        q.value <- p.adjust(tmp$p.value, method = p.adj)
        tmp$q.value <- q.value

        gsea[[paste0("pc", i)]] <- tmp[with(tmp, order(q.value,
                                                       -ks.statistic)), ]

    }

    gsea

}

.rankTest <- function(signature, loadings, ranking) {

    ind <- which(names(loadings) %in% signature)
    geneset <- ranking[ind]
    background <- ranking[-ind]

    ks <- ks.test(geneset, background)
    c(ks$statistic, ks$p.value)
}

.filterSignatures <- function(signature, genes){
    signature[signature %in% genes]
}

