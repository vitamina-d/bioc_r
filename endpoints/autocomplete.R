library(org.Hs.eg.db)
library(AnnotationDbi)

#* entrez devuelve el entrez si lo encuentra
#* @param input 
#* @get /
#* @tag detail
#* @serializer unboxedJSON 
function(input, res) {

    tryCatch({
        alias <- unique(keys(org.Hs.eg.db, keytype = "ALIAS"))
        matches <- grep(input, alias, ignore.case = TRUE, value = TRUE)

        if (is.null(matches) || length(matches) == 0 ) {
            res$status <- 404
            print("ACA EL 404")
            stop(paste("No se encontro el alias: "), call. = FALSE)
        } else if (length(matches) == 1) {
            matches <- list(matches)
        } 

        list(
            code = 200,
            message = "Ok",
            data = head(matches, 20)
        )

    }, error = function(e) {
        res$status <- 500
        stop(paste("Error de servicio R: ", e$message), call. = FALSE)
    })
}