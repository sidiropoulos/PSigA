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
#' data("AML")
#' data("RAPIN")
#'
#' scores <- scoreSigs(AML, parsed = TRUE, signatures = RAPIN,
#'                     threshold = 0.003)
#' head(scores)
#'
#' #plot top signature
#' sigPCA <- signaturePCA(RAPIN[[rownames(scores)[1]]], AML)
#' sigBiplot(sigPCA, AML_meta$karyotype, main = rownames(scores)[1])
#'
#' @import parallel
#' @export
scoreSigs <- function(data, parsed = FALSE, genes = (if (parsed) NULL),
                      signatures,  threshold = 0.005, n = 200){

    if (!parsed)
        data <- readSamples(data, genes)

    scores <- do.call(rbind, mclapply(signatures, peakDistance2d, data,
                                      threshold, n, mc.preschedule = FALSE))

    scores <- as.data.frame(scores, stringsAsFactors = FALSE)
    colnames(scores) <- c("score", "size")
    s <- sort(scores$score, decreasing = TRUE, index.return = TRUE)
    scores <- scores[s$ix, ]


    if (parsed)
        scores
    else
        list(data = data, scores = scores)
}
