#* entrezBySymbol devuelve el entrez si lo encuentra
#* @param symbol
#* @get /
#* @tag entrez
#* @serializer unboxedJSON 
function(symbol) {

    if (is.null(symbol) || symbol == "" ) {
        result <- list(
            code = 400,
            message = "Ingrese un valor.",
            data = NULL
        )
        return(result)
    } 

    entrez <- tryCatch({
        AnnotationDbi::select(org.Hs.eg.db, keys = symbol, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID
    }, error = function(e) NULL)

    if (is.null(entrez)) {
        result <- list(
            code = 500,
            message = "Error del servidor.",
            data = NULL
        )
    } else if (length(entrez) == 0) {
        result <- list(
            code = 404,
            message = paste("no se encontro entrez para ", symbol),
            data = NULL
        )
    } else {
        result <- list(
            code = 200,
            message = "Ok",
            data = NULL
        )
    }
    return(result)
}