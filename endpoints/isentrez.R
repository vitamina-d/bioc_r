library(AnnotationDbi)
library(org.Hs.eg.db)  

#* isentrez valida que el valor sea un ENTREZ
#* @param entrez
#* @get /
#* @tag entrez
#* @serializer unboxedJSON 
function(entrez, res) {
    is_entrez <- FALSE

    ids <- tryCatch({
        keys(org.Hs.eg.db, keytype = "ENTREZID")    
    }, error = function(e) {
        res$status <- 500
        stop(paste("Fallo la consulta a org.Hs.eg.db: ", e$message), call. = FALSE)
    })

    if (entrez %in% ids) {
        is_entrez = TRUE
    }

    result <- list(
        code = 200,
        message = "Ok.",
        data = list(
            is_entrez = is_entrez
        )
    )
}
