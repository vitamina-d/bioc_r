#* entrezByAlias devuelve el entrez si lo encuentra
#* @param alias puede ser un alias
#* @get /
#* @tag endpoints
#* @serializer unboxedJSON 
function(alias) {

    start_time <- Sys.time()

    if (is.null(alias) || alias == "" ) {
        result <- list(
            code = 400,
            time_secs = time,
            data = list (
                message = "Ingrese un valor."
            )
        )
        return(result)
    } 

    entrez <- tryCatch({
        AnnotationDbi::select(org.Hs.eg.db, keys = alias, columns = "ENTREZID", keytype = "ALIAS")$ENTREZID
    }, error = function(e) NULL)
    
    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    if (is.null(entrez)) {
        result <- list(
            code = 500,
            time_secs = time,
            data = list (
                message = "try catch"
            )
        )
    } else if (length(entrez) == 0) {
        result <- list(
            code = 404,
            time_secs = time,
            data = list (
                message = paste("no se encontro entrez para ", alias)
            )
        )
    } else {
        result <- list(
            code = 200,
            time_secs = time,
            data = list( 
                entrez = unique(entrez)
            )
        )
    }
    return(result)
}