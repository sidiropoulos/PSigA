#' Parse data to be used \code{\link{PsigA}} functions
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
#' signatures that will be used in the \code{PsigA} functions.
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
#' #use fData() to get the gene symbols for vdx.
#' VDXparsed <- parseData(data = exprs(vdx), geneIds = fData(vdx)$Gene.symbol)
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
#' keys <- AnnotationDbi::select(hgu95av2.db, rownames(ALLexprs),
#'                              "SYMBOL", "PROBEID")
#' #remove probe duplicates
#' geneIds <- keys$SYMBOL[ !duplicated(keys$PROBEID) ]
#'
#' ALLparsed <- parseData(ALLexprs, geneIds)
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise_each
#' @importFrom dplyr funs
#' @export
parseData <- function(data, geneIds) {

    #check if data is data.frame
    if (!is.data.frame(data))
        data <- data.frame(data)

    data$geneIds <- geneIds

    #find and remove NA values
    data <- data[which(!is.na(geneIds)), ]

    data <- summarise_each(group_by(data, geneIds), funs(mean))
    data <- as.data.frame(data)

    rownames(data) <- data$geneIds
    data <- subset(data, select = - geneIds)

    data
}

