library(org.Hs.eg.db)
library(AnnotationDbi)

#* entrez devuelve el entrez si lo encuentra
#* @param input 
#* @get /
#* @tag detail
#* @serializer unboxedJSON 
function(input) {

    alias <- unique(keys(org.Hs.eg.db, keytype = "ALIAS"))

    matches <- grep(input, alias, ignore.case = TRUE, value = TRUE)

    if (length(matches) == 0)  {
        matches <- NULL
    } 
    
    if (length(matches) == 0) {
        result <- list(
            code = 404,
            message = paste("no se encontro"),
            data = list()
        )
    } else {
        result <- list(
            code = 200,
            message = "Ok",
            data = head(matches, 20)
        )
    }
    return(result)
}