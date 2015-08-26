#' SigPCA: A gene-signature based exploratory tool for gene expression data.
#'
#' The SigPCA package provides a set of functions to perform Principal Component analysis based on genetic pathways.
#' Moreover SigPCA incorporates an algorithm to measure the spread of the data in the PCA space and rank the
#' gene-signatures based on their ability to split the data in distinct entities.
#'
#' @docType package
#' @name SigPCA
NULL

#' Microarray data of 183 AML patients
#'
#' Gene expression data of 183 AML patients, retrieved from TCGA, normalized and batch corrected.
#' Rows represent microarray probes of the Affymetrix HGU 133 Plus2 platform. Each column represents
#' an AML patient.
#'
#' @usage data(AML)
#' @format A data frame with 54675 rows (Affy HGU 133 Plus 2 probes) and 183 columns (patients).
#'
#' @source \url{https://tcga-data.nci.nih.gov/tcga/}
"AML"

#' AML dataset metadata
#'
#' Clinical data of the TCGA AML patients.
#'
#' @usage data(AML_meta)
#' @format A data frame with 183 rows (patients) and 5 variables
#' \itemize{
#'     \item karyotype: Cytogenetic aberration
#'     \item gender: 1 = male, 0 = female
#'     \item age: patient's age in years
#'     \item OS: Survival from diagnosis until the last follow up in days
#'     \item EventOS: Patient status by the time of the last follow up. 1 = dead, 0 = alive
#' }
#' @source \url{https://tcga-data.nci.nih.gov/tcga/}
"AML_meta"

#' Bone marrow samples from 60 healthy individuals
#'
#' Gene expression data of 60 FACS sorted blood samples using the Affymetrix HGU 133 plus2 platform.
#'
#' @usage data(blood)
#' @format A data frame with 54675 rows (probes) and 60 columns. The column names describe the blood cell type of
#' the respective sample.
"blood"

#' 12,617 gene signatures
#'
#' Gene pathways were obtained from Molecular Signature Database (MSigDB), Connectivity Map (cmap)
#' and the signatures by Rapin et. al.
#'
#' @usage data(sigdb)
#' @format A list with 12617 entries, each one representing a gene signature. The name of each list entry represents
#' the signature's name.
#' @source \url{http://www.broadinstitute.org/gsea/msigdb/index.jsp}
#' @source \url{https://www.broadinstitute.org/cmap/}
#' @source \url{http://www.bloodjournal.org/cgi/pmidlookup?view=long&pmid=24363398}
"sigdb"
