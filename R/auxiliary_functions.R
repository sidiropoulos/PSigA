#' @title PCA based on a given gene-signature
#'
#' @description \code{signaturePCA} performs principal component analysis on the given data matrix and gene-signature and
#' returns the results as an object of class \code{prcomp}.
#'
#' @author Nikos Sidiropoulos
#'
#' @param genes Character vector with the gene identifiers
#' @param data Data matrix with gene expression values
#' @param db Microarray Platform
#' @param center a logical value indicating whether the variables should be shifted to be zero centered. See
#' \code{\link{prcomp}} for more details.
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance before
#' the analysis takes place. See \code{prcomp} for more details.
#' @param ... arguments passed to \code{prcomp}.
#'
#' @export signaturePCA
signaturePCA <- function(genes, data, db, center = TRUE, scale = FALSE, ...) {

  #Convert geneID to probeID
  keys <- reduceProbes( gene2probe(genes, db), data[])
  prcomp(t(data[keys,]), center = center, scale = scale, ...)
}

#' @title  Convert gene identifiers to microarray probeID
#'
#' @description \code{gene2probe} returns probeIDs compatible with several Microarray Platforms
#'
#' @author Nikos Sidiropoulos
#'
#' @param genes List of gene identifiers to be converted
#' @param db Microarray platform. See details for a list of valid inputs
#' @param keytype the keytype that matches the keys used in the \code{gene} parameter. See a list with valid inputs
#' in the details below. (Default = "SYMBOL")
#'
#' @return keys a data frame of possible values
#'
#' @import AnnotationDbi hgu133plus2.db
#' @export gene2probe
gene2probe <- function(genes, db, keytype = "SYMBOL") {

#    keys <- try(suppressWarnings(AnnotationDbi::select( x = get(db), keys = genes, columns = "PROBEID" , keytype)), TRUE)

    keys <- AnnotationDbi::select( x = get(db), keys = genes, columns = "PROBEID" , keytype)

#     #If ALL provided genes don't map to any probes return NULL
#     if ( inherits( keys, "try-error" ) ) {
#         return(NULL)
#     }

    #Remove genes that didn't map to any probes
    keys <- keys[ ! is.na(keys[,2]), ]

    keys

}

#' @title For each gene select the probe with the highest mean expression across the input sample.
#'
#' @description \code{reduceprobes} finds and selects the probe with the highest mean expression among the probes of the
#' same gene.
#'
#' @author Nikos Sidiropoulos
#'
#' @param keys data frame of probe identifiers. Output of \code{\link{gene2probe}} function.
#' @param data gene expression matrix
#'
#' @export reduceProbes
reduceProbes <- function(keys, data) {

  #if ( is.null( keys ) ) { return(NULL) }

  #Initialize vector of probes to keep
  kept <- character(0)

  for (gene in unique(keys[,1]) ) {

    geneProbes <- keys[ keys[,1] == gene, 2]

    # Find the probe with the highest mean expression
    keep <- names( which.max(abs( apply( data[ geneProbes, ], 1, mean) ) ) )

    kept <- c( kept, keep )
  }

  return(kept)

}

#' @title Converts probe id to gene symbols
#'
#' @description Converts probe id to gene symbols
#'
#' @author Nikos Sidiropoulos
#'
#' @param keys a character vector of probe identifiers
#' @param db microarray platform
#' @return character vector of gene symbols
#'
#' @import AnnotationDbi hgu133plus2.db
#' @export probe2geneMap
probe2geneMap <- function(keys, db){

    keys <- AnnotationDbi::select(get(db), keys, columns = "SYMBOL", keytype = "PROBEID")

    #Remove non mapped and get the unique ones

    genes <- c()

    for (probe in unique(keys$PROBEID)) {
        genes <- c(genes, paste( keys$SYMBOL[ keys$PROBEID == probe ], collapse = "/"))
    }

    genes
}


