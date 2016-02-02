#' PSigA: A gene-signature ranking method based on sample density in PCA space.
#'
#' The PSigA package provides a set of functions to perform Principal Component
#' analysis based on genetic pathways. Moreover PSigA incorporates an algorithm
#' to measure the spread of the data in the PCA space and rank the
#' gene-signatures based on their ability to split the data in distinct
#' entities.
#'
#' @docType package
#' @name PSigA
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
