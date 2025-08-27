#* entrezBySymbol devuelve el entrez si lo encuentra
#* @param alias 
#* @get /
#* @serializer unboxedJSON 
function(alias) {    #################################### no acepta numeros ni minusculas

    start_time <- Sys.time()

    if (is.null(alias) || alias == "" ) {
        end_time <- Sys.time()
        time <- as.numeric(difftime(end_time, start_time, units = "secs"))

        result <- list(
            code = 400,
            datetime = start_time,
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
            datetime = start_time,
            time_secs = time,
            data = list (
                message = "Error del servidor."
            )
        )
    } else if (length(entrez) == 0) {
        result <- list(
            code = 404,
            datetime = start_time,
            time_secs = time,
            data = list (
                message = paste("no se encontro entrez para ", symbol)
            )
        )
    } else {
        result <- list(
            code = 200,
            datetime = start_time,
            time_secs = time,
            data = list( 
                entrez = unique(entrez)
            )
        )
    }
    return(result)
}