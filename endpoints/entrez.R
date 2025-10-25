#* entrez devuelve el entrez de un ALIAS o SYMBOL (enviar mayus)
#* @param value 
#* @get /
#* @tag entrez
#* @serializer unboxedJSON 
function(value, res) {

    if (is.null(value) || value == "" ) {
        res$status <- 400
        stop(paste("Ingrese value: "), call. = FALSE)
    } 

    entrez <- tryCatch({
        value <- toupper(value)
        AnnotationDbi::select(org.Hs.eg.db, keys = value, columns = "ENTREZID", keytype = "ALIAS")$ENTREZID
    }, error = function(e) NULL)

    if (is.null(entrez) || length(entrez) == 0) {
        entrez <- tryCatch({
            AnnotationDbi::select(org.Hs.eg.db, keys = value, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID
        }, error = function(e) NULL)
    }

    if (is.null(entrez) || length(entrez) == 0) {
        res$status <- 404
        stop(paste("No se encontró el ENTREZID para el value: ", value), call. = FALSE)
    }

    list(
        code = 200,
        message = "Ok",
        data = list( 
            entrez = unique(entrez)
        )
    )
}