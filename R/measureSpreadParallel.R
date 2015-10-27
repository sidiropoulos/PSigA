#' @title Parallel execution of \code{\link{measureSpread}}
#'
#' @description \code{measureSpreadParallel}
#'
#' @author Nikos Sidiropoulos
#'
#' @param data gene expression matrix
#' @param signatures list of gene signatures
#' @param threads number of cpu threads to be used
#' @param ... arugments passed to \code{\link{measureSpread}}
#'
#' @return a data frame with columns: Signature Name, Number of principal
#' components used, Density area.
#'
#' @import foreach
#' @import doMC
#' @export measureSpreadParallel


measureSpreadParallel <- function(data, signatures, db, threads = 4,
                                  ...){

    registerDoMC(threads)

    spread <- foreach(i=1:length(signatures), .combine = rbind) %dopar% {

        .prepare(data, signatures[[i]], db, names(signatures[i]), ...)
        }
    spread
}


.prepare <- function(data, signature, db, name, ...){
    sigpca <- signaturePCA(signature$genes, data, db = db)
    measureSpread(sigpca, sigName = name, categoryID = signature$collectionId,
                  ...)
}
