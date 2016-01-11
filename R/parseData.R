#' Parse data to be used \code{\link{PSigA}} functions
#'
#' @description The function converts the gene identifiers (rownames) of a gene
#' expression matrix to a desired format, supplied by the parameter
#' \code{geneIds}. If any of the gene identifiers correspond to more than one
#' row of the expression matrix (e.g. probes that map to the same gene), the
#' median expression of each patient..
#'
#' @param data a data frame of matrix with gene expression values where rows
#' represent genes and columns represent samples.
#' @param geneIds vector of length \code{nrow(data)} containing gene names in
#' of the desired format. The format must be the same as the format of the
#' signatures that will be used in the \code{PSigA} functions.
#'
#' @return An object of class \code{data} where the rows correspond to
#' \code{unique(geneIds)}.
#'
#' @examples
#'
#' require(Biobase)
#' require(breastCancerVDX)
#'
#' data(vdx)
#'
#' #use fData() to get the gene symbols for vdx. We use only the first 100
#' #probes as an example.
#' VDXparsed <- parseData(data = exprs(vdx)[1:100,],
#'                        geneIds = fData(vdx)$Gene.symbol[1:100])
#'
#'
#' require(ALL)
#' data(ALL)
#'
#' ALLexprs <- exprs(ALL)
#'
#' #featureData ("fData()") for ALL are missing so we need to get the gene
#' #symbols manually.
#'
#' annotation(ALL)
#'
#' require(hgu95av2.db)
#'
#' keys <- AnnotationDbi::select(hgu95av2.db, rownames(ALLexprs)[1:100],
#'                              "SYMBOL", "PROBEID")
#' #remove probe duplicates
#' geneIds <- keys$SYMBOL[ !duplicated(keys$PROBEID) ]
#'
#' ALLparsed <- parseData(ALLexprs, geneIds)
#'
#' @import parallel
#' @export
parseData <- function(data, geneIds) {

    #find and remove NA values
    nonNAgenes <- which(!is.na(geneIds))
    data <- data[nonNAgenes,]
    geneIds <- geneIds[nonNAgenes]

    #select the probe/isoform with the highest median expression
    d <- do.call(rbind, mclapply(unique(geneIds), .maxMedian, data, geneIds))

    rownames(d) <- unique(geneIds)
    d
}

.maxMedian <- function(targetGene, data, genes) {

    match <- genes == targetGene
    if (sum(match) == 1){
        data[match,]
    }else {
        apply(data[match,], 2, median)

    }

}
