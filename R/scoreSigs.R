#' @title Score gene signatures based on PCA density
#'
#' @param data data frame of matrix with gene expression values where rows
#' represent genes and columns represent samples.
#' @param parsed logical.
#' @param genes vector of length \code{nrow(data)} containing gene names in
#' "HGNC" format. If \code{parsed = TRUE} the parameter is omitted.
#' @param signatures list where every entry is a character vector of
#' the gene names in HGNC format that correspond to a given gene-pathway.
#' @param threshold low density cutoff. Used in \code{\link{peakDistance2d}}.
#' @param n Number of grid points in each direction. Can be scalar or a
#' length-2 integer vector. See \link{kde2d}.
#'
#' @return A data frame with 2 columns; 1) \code{score}, the \code{PSigA}
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
#' VDX <- readSamples(data = exprs(vdx), genes = fData(vdx)$Gene.symbol)
#'
#' scores <- scoreSigs(VDX, parsed = TRUE, signatures = MSigDB[1827:1838],
#'                     threshold = 0.003)
#' head(scores)
#'
#' #plot top signature
#' sigPCA <- signaturePCA(MSigDB[[rownames(scores)[1]]], VDX)
#' sigPlot(sigPCA)
#' sigBiplot(sigPCA, factor(pData(vdx)$grade), main = rownames(scores)[1])
#'
#' @import parallel
#' @export
scoreSigs <- function(data, geneIds, signatures, threshold = 0.005, n = 200){


    #find and remove NA values
    nonNAgenes <- which(!is.na(geneIds))
    data <- data[nonNAgenes,]
    genes <- geneIds[nonNAgenes]

    if (sum(duplicated(genes)) > 0){
        message(strwrap("Duplicate genes found. The row with the highest median
                        expresison will be selected."))
        data <- do.call(rbind, mclapply(unique(genes), .maxMedian, data,
                                        genes))
        rownames(data) <- unique(genes)
    }

    scores <- do.call(rbind, mclapply(signatures, peakDistance2d, data,
                                      threshold, n))

    scores <- as.data.frame(scores, stringsAsFactors = FALSE)
    colnames(scores) <- c("Score", "Genes.found")
    s <- sort(scores$score, decreasing = TRUE, index.return = TRUE)
    scores <- scores[s$ix, ]


    if (length(unique(genes)) == length(geneIds))
        scores
    else
        list(data = data, scores = scores)
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
