#' @title PCA based on a given gene-signature
#'
#' @description \code{signaturePCA} performs principal component analysis on the given data matrix and gene-signature and
#' returns the results as an object of class \code{prcomp}.
#'
#' @param genes Character vector with the gene identifiers
#' @param data Data matrix with gene expression values
#' @param db Microarray Platform
#' @param ... arguments passed to \code{\link{prcomp}}.
#'
#' @export signaturePCA
signaturePCA <- function(genes, data, db, ...) {

  #Convert geneID to probeID
  keys <- reduceProbes( gene2probe(genes, db), data[])
  prcomp(t(data[keys,]), center = TRUE, scale = FALSE, ...)
}

#' @title  Convert gene identifiers to microarray probeID
#'
#' @description \code{gene2probe} returns probeIDs compatible with several Microarray Platforms
#'
#' @param genes List of gene identifiers to be converted
#' @param db Microarray platform. See details for a list of valid inputs
#' @param keytype the keytype that matches the keys used in the \code{gene} parameter. See a list with valid inputs
#' in the details below. (Default = "SYMBOL")
#'
#' @return keys a data frame of possible values
#'
#' @import hgu133plus2.db
#' @export
gene2probe <- function(genes, db, keytype = "SYMBOL") {

  keys <- try(suppressWarnings(select( x = get(db), keys = genes, columns = "PROBEID" , keytype)), TRUE)

  #If ALL provided genes don't map to any probes return NULL
  if ( inherits( keys, "try-error" ) ) {
    return(NULL)
  }

  #Remove genes that didn't map to any probes
  keys <- keys[ ! is.na(keys[,2]), ]

  keys

}

#' @title For each gene select the probe with the highest mean expression across the input sample.
#'
#' @description \code{reduceprobes} finds and selects the probe with the highest mean expression among the probes of the
#' same gene.
#'
#' @param keys data frame of probe identifiers. Output of \code{\link{gene2probe}} function.
#' @param data gene expression matrix
#'
#' @export reduceProbes
reduceProbes <- function(keys, data) {

  if ( is.null( keys ) ) { return(NULL) }

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
#' @param keys a character vector of probe identifiers
#' @param db microarray platform
#' @return character vector of gene symbols
#'
#' @import hgu133plus2.db
#' @export probe2geneMap
probe2geneMap <- function(keys, db){

    keys <- select(get(db), keys, "PROBEID", "SYMBOL")
    #Remove non mapped and get the unique ones
    unique(keys[!is.na(keys[,2]),1])
}

