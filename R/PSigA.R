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

#' Microarray data of 424 AML patients
#'
#' Gene expression data of 424 AML patients with aberrant karyotyope (AK-AML),
#' retrieved from TCGA, GSE15434, GSE13159 and GSE14468, normalized and batch
#' corrected. Rows represent genes in HGNC format. Each column represents an
#' AML patient.
#' A data frame with clinical information on the karyotype, gender and age of
#' the AK-AML patients is also loaded (\code{\link{AML_meta}}).
#'
#' @usage data(AML)
#' @format
#' A data frame with 21369 rows (genes in HGNC format) and 424
#' columns (patients).
#'
#' @source \url{https://tcga-data.nci.nih.gov/tcga/}
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15434}
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13159}
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14468}
"AML"

#' AML dataset metadata
#'
#' Clinical data for the \code{\link{AML}} dataset.
#'
#' @usage data(AML)
#' @format A data frame with 424 rows (patients) and 3 columns:
#' \itemize{
#'     \item karyotype: Cytogenetic aberration
#'     \item gender: 1 = male, 0 = female
#'     \item age: patient's age in years
#' }
#' @source \url{https://tcga-data.nci.nih.gov/tcga/}
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15434}
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13159}
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14468}
"AML_meta"

#' Molecular Signature DataBase
#'
#' 10,348 gene pathways obtained from Molecular Signature Database
#' (MSigDB v5.0).
#'
#' @usage data(MSigDB)
#' @format An object of type list with 10,348 entries, each one representing a
#' gene signature. The name of each list entry represents the signature's name.
#' Additionally, \code{MSigDB} entries contain the following collection
#' identifiers:
#' \itemize{
#'  \item{H: hallmark gene sets (50)}
#'  \item{C1: positional gene sets (326)}
#'  \item{C2: curated gene sets (4725)}
#'  \item{C3: motif gene sets (836)}
#'  \item{C4: computational gene sets (858)}
#'  \item{C5: GO gene sets (1454)}
#'  \item{C6: oncogenic signatures (189)}
#'  \item{C7: immunologic signatures (1910)}
#' }
#'
#' @source \url{http://www.broadinstitute.org/gsea/msigdb/index.jsp}
"MSigDB"

#' Rapin et al. signatures
#'
#' 24 gene pathways obtained from Rapin et al. (2014).
#'
#' @usage data(RAPIN)
#' @format An object of type list with 24 entries, each one representing a
#' gene signature. The name of each list entry represents the signature's name.
#'
#' @source \url{http://www.bloodjournal.org/content/123/6/894}
"RAPIN"
