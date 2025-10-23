library(AnnotationDbi)
library(org.Hs.eg.db)  

#* isentrez valida que el valor sea un ENTREZ
#* @param entrez
#* @get /
#* @tag entrez
#* @serializer unboxedJSON 
function(entrez) {
    
    if (is.null(entrez) || entrez == "" ) {
        result <- list(
            code = 400,
            message = "Ingrese un valor.",
            data = NULL
        )
        return(result)
    } 

    ids <- tryCatch({
        keys(org.Hs.eg.db, keytype = "ENTREZID")    
    }, error = function(e) NULL)

    if (is.null(ids)) {
        result <- list(
            code = 500,
            message = "Error del servidor.",
            data = NULL
        )
    }

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
