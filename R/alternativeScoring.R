bicScore <- function(signature, data, G){

    inSet <- signature %in% rownames(data)

    if (sum(inSet) < 3)
        return("NA")

    sigpca <- prcomp(t(data[signature[inSet], ]), center = TRUE,
                     scale = FALSE)

    mod <- Mclust(sigpca$x[,1:2], G)

    c(mod$bic, mod$G, mod$modelName, sum(inSet))
}

kmeansScore <- function(signature, data, k){

    inSet <- signature %in% rownames(data)

    if (sum(inSet) < 3)
        return(c(-1,0))

    sigpca <- prcomp(t(data[signature[inSet], ]), center = TRUE,
                     scale = FALSE)

    a <- kmeans(sigpca$x[,1:2], centers = k)

    c(a$totss, a$betweenss, sum(inSet))
}
