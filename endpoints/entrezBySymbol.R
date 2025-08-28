#* entrezBySymbol devuelve el entrez si lo encuentra
#* @param symbol
#* @get /
#* @tag entrez
#* @serializer unboxedJSON 
function(symbol) {

    start_time <- Sys.time()

    if (is.null(symbol) || symbol == "" ) {
        end_time <- Sys.time()
        time <- as.numeric(difftime(end_time, start_time, units = "secs"))

        result <- list(
            code = 400,
            message = "Ingrese un valor.",
            datetime = start_time,
            time_secs = time,
            data = NULL
        )
        return(result)
    } 

    entrez <- tryCatch({

        AnnotationDbi::select(org.Hs.eg.db, keys = symbol, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID

    }, error = function(e) NULL)

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    if (is.null(entrez)) {
        result <- list(
            code = 500,
            message = "Error del servidor.",
            datetime = start_time,
            time_secs = time,
            data = NULL
        )
    } else if (length(entrez) == 0) {
        result <- list(
            code = 404,
            message = paste("no se encontro entrez para ", symbol),
            datetime = start_time,
            time_secs = time,
            data = NULL
        )
    } else {
        result <- list(
            code = 200,
            message = "Ok",
            datetime = start_time,
            time_secs = time,
            data = NULL
        )
    }
    return(result)
}