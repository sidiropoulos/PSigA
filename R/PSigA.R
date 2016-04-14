#' PsigA: A gene-signature ranking method based on sample density in PCA space.
#'
#' The PsigA package provides a set of functions to perform Principal Component
#' analysis based on genetic pathways. Moreover PSigA incorporates an algorithm
#' to measure the spread of the data in the PCA space and rank the
#' gene-signatures based on their ability to split the data in distinct
#' entities.
#'
#' @docType package
#' @name PsigA
NULL

#' Rapin et al. signatures
#'
#' 21 gene pathways obtained from Rapin et al. (2014).
#'
#' @usage data(RAPIN)
#' @format An object of type list with 21 entries, each one representing a
#' gene signature. The name of each list entry represents the signature's name.
#'
#' @return A list of 21 character vectors, each one containing gene names in
#' HGNC format.
#'
#' @source \url{http://www.bloodjournal.org/content/123/6/894}
"RAPIN"

#' @title Score gene signatures based on PCA density
#'
#' @description PsigAscore blah blah
#'
#' @param data data frame of matrix with gene expression values where rows
#' represent genes and columns represent samples.
#' @param signatures list where every entry is a character vector of
#' the gene names that correspond to a gene-pathway. If the gene format is
#' different from the one in \code{rownames(data)}, use \code{\link{parseData}}
#' first.
#' @param threshold low density cutoff. Used in \code{\link{peakDistance2d}}.
#' @param n Number of grid points in each direction. Can be scalar or a
#' length-2 integer vector. See \link{kde2d}.
#' @param magnitude When TRUE the score is multiplied by the cluster density.
#' Default: FALSE.
#' @param scale a logical value indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place.
#' See prcomp for more details.
#'
#' @return A data frame with 2 columns; 1) \code{score}, the \code{PsigA}
#' score of a given signature and 2) \code{size}, the number of genes in the
#' \code{signature} that were found in the \code{data}. If
#' \code{parsed = FALSE}, the function returns a list, where the first entry
#' contains the data frame described above and the second entry the parsed
#' \code{data}.
#'
#' @examples
#'
#' require(Biobase)
#' require(breastCancerVDX)
#' data(vdx)
#' data(RAPIN)
#'
#' VDX <- parseData(data = exprs(vdx), geneIds = fData(vdx)$Gene.symbol)
#'
#' scores <- PsigAscore(VDX, RAPIN, threshold = 0.003)
#' head(scores)
#'
#' #plot top signature
#' sigPlot(VDX, RAPIN[[rownames(scores)[1]]])
#' sigBiplot(VDX, RAPIN[[rownames(scores)[1]]], factor(pData(vdx)$grade),
#'           main = rownames(scores)[1])
#'
#' @import parallel
#' @importFrom knitr knit
#' @importFrom knitr pandoc
#' @export
PsigA <- function(data, signatures, threshold = 0.005, n = 200,
                       magnitude = FALSE, scale = FALSE, report = TRUE,
                       p.adj = p.adjust.methods, min.length = 10,
                       max.length = 500){

    #check input
    p.adj <- match.arg(p.adj)

    #Filter signatures
    message("Filtering signatures... ", appendLF = FALSE)
    signatures <- lapply(signatures, .filterSignatures, rownames(data))

    sig.lengths <- sapply(signatures, length)
    keep <- sig.lengths >= min.length & sig.lengths <= max.length

    message("Done! ", appendLF = FALSE)
    message(paste0("(", sum(keep), "/", length(signatures), " passed)"))

    #Call scoring function
    message("Computing PsigA scores... ", appendLF = FALSE)
    scores <- do.call(rbind, mclapply(signatures[keep], peakDistance2d, data,
                                      threshold, n, magnitude, scale,
                                      filtered = TRUE))

    scores <- as.data.frame(scores, stringsAsFactors = FALSE)
    colnames(scores) <- c("Score", "Length")
    scores <- scores[with(scores, order(-Score)), ]

    message("Done!")

    #Perform GSEA
    message("Performing GSEA... ", appendLF = FALSE)
    gsea <- rankedGSEA(data, signatures[keep], p.adj, filtered = TRUE)
    message("Done!")

    if (report) {
        f <- system.file("extdata", "knitr_test.Rmd", package = "PsigA")

        message("Generating HTML report... ", appendLF = FALSE)
        pandoc(knit(f, quiet = TRUE), format = "html")
        message("Done!")

        browseURL("knitr_test.html")
    }

    list(scores = scores, gsea = gsea)
}
