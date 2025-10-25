library(AnnotationDbi)
library(org.Hs.eg.db)  

#* detail by ENTREZ 
#* @param entrez:int id del gen
#* @get /
#* @tag detail
#* @serializer unboxedJSON 
function(entrez, res) {

    detail <- tryCatch({
        AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("GENETYPE", "GENENAME", "SYMBOL", "ALIAS"), keytype = "ENTREZID")
    }, error = function(e) {
        res$status <- 500
        stop(paste("Error de servicio R: ", e$message), call. = FALSE)
    })

    if (is.null(detail) || length(detail) == 0 || nrow(detail) == 0) {
        res$status <- 404
        stop(paste("No se encontro el entrez: ", e$message), call. = FALSE)
    }

    aliases <- unique(detail$ALIAS)
    print(aliases)
    if (length(aliases) == 0) {
        aliases <- list()
    } else if (length(aliases) == 1) {
        aliases <- list(aliases)
    } 

    #select 1: ok
    #select 1:1 da error
    response <- list(
        code = 200,
        message = "Ok",
        data = list(
            entrez = entrez,
            symbol = unique(detail$SYMBOL),
            genename = unique(detail$GENENAME),
            genetype = unique(detail$GENETYPE),
            alias = aliases
        )
    )
}
