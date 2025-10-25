library(org.Hs.eg.db)
library(AnnotationDbi)

#* entrez devuelve el entrez si lo encuentra
#* @param input 
#* @get /
#* @tag detail
#* @serializer unboxedJSON 
function(input) {

    tryCatch({
        alias <- unique(keys(org.Hs.eg.db, keytype = "ALIAS"))
        matches <- grep(input, alias, ignore.case = TRUE, value = TRUE)

        if (length(matches) == 0) {
            matches <- list()
        } else if (length(matches) == 1) {
            matches <- list(matches)
        } 
        #print(matches)
        #print(typeof(matches))

        list(
            code = 200,
            message = "Ok",
            data = head(matches, 20)
        )

    }, error = function(e) {
        # exception
        stop(paste("No se obtuvo detalle: ", e$message), call. = FALSE)
    })
}