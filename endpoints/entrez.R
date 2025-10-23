#* entrez devuelve el entrez si lo encuentra
#* @param value 
#* @get /
#* @tag entrez
#* @serializer unboxedJSON 
function(value) {

    if (is.null(value) || value == "" ) {
        result <- list(
            code = 400,
            data = list (
                message = "Ingrese un valor.",
                entrez = NULL
            )
        )
        return(result)
    } 

    entrez <- tryCatch({
        AnnotationDbi::select(org.Hs.eg.db, keys = value, columns = "ENTREZID", keytype = "ALIAS")$ENTREZID
    }, error = function(e) NULL)

    if (is.null(entrez)) {
        entrez <- tryCatch({
            AnnotationDbi::select(org.Hs.eg.db, keys = value, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID
        }, error = function(e) NULL)
    } 
    
    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    if (is.null(entrez) || length(entrez) == 0) {
        result <- list(
            code = 404,
            message = paste("no se encontro entrez para ", value),
            data = NULL
        )
    } else {
        result <- list(
            code = 200,
            message = "Ok",
            data = list( 
                entrez = unique(entrez)
            )
        )
    }
    return(result)
}