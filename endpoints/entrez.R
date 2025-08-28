#* entrez devuelve el entrez si lo encuentra
#* @param value 
#* @get /
#* @tag entrez
#* @serializer unboxedJSON 
function(value) {

    start_time <- Sys.time()

    if (is.null(value) || value == "" ) {
        end_time <- Sys.time()
        time <- as.numeric(difftime(end_time, start_time, units = "secs"))

        result <- list(
            code = 400,
            datetime = start_time,
            time_secs = time,
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
            datetime = start_time,
            time_secs = time,
            data = list (
                message = paste("no se encontro entrez para ", value),
                entrez = NULL
            )
        )
    } else {
        result <- list(
            code = 200,
            datetime = start_time,
            time_secs = time,
            data = list( 
                message = "Ok",
                entrez = unique(entrez)
            )
        )
    }
    return(result)
}