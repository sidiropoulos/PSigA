#' SigPCA: A gene-signature based exploratory tool for gene expression data.
#'
#' The SigPCA package provides a set of functions to perform Principal Component
#' analysis based on genetic pathways. Moreover SigPCA incorporates an algorithm
#' to measure the spread of the data in the PCA space and rank the
#' gene-signatures based on their ability to split the data in distinct
#' entities.
#'
#' @docType package
#' @name SigPCA
NULL

#' Microarray data of 183 AML patients
#'
#' Gene expression data of 183 AML patients, retrieved from TCGA, normalized
#' and batch corrected. Rows represent microarray probes of the Affymetrix HGU
#' 133 Plus2 platform. Each column represents an AML patient.
#'
#' @usage data(AML)
#' @format A data frame with 54675 rows (Affy HGU 133 Plus 2 probes) and 183
#' columns (patients).
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
#'     \item EventOS: Patient status by the time of the last follow up.
#'           1 = dead, 0 = alive
#' }
#' @source \url{https://tcga-data.nci.nih.gov/tcga/}
"AML_meta"

#' Bone marrow samples from 60 healthy individuals
#'
#' Gene expression data of 60 FACS sorted blood samples using the Affymetrix
#' HGU 133 plus2 platform.
#'
#' @usage data(blood)
#' @format A data frame with 54675 rows (probes) and 60 columns. The column
#' names describe the blood cell type of the respective sample.
"blood"

#' Molecular Signature DataBase
#'
#' 10,348 gene pathways obtained from Molecular Signature Database
#' (MSigDB v5.0).
#'
#' @usage data(MSigDb)
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

#' Drug SIGnature DataBase
#'
#' 22,528 gene pathways obtained from Drug SIGnature DataBase (DSigDB v1.0).
#'
#' @usage data(DSigDb)
#' @format An object of type list with 22,528 entries, each one representing a
#' gene signature. The name of each list entry represents the signature's name.
#' Additionally, \code{MSigDB} entries contain the following collection
#' identifiers:
#' \itemize{
#'  \item{D1: FDA Approved (1202)}
#'  \item{D2: Kinase Inhibitors (1220)}
#'  \item{D3: Perturbagen Signatures (1998)}
#'  \item{D4: Computational Drug Signatures (18107)}
#' }
#'
#' @source \url{http://tanlab.ucdenver.edu/DSigDB/DSigDBv1.0/}
"DSigDB"

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
