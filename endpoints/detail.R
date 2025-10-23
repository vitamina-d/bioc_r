library(AnnotationDbi)
library(org.Hs.eg.db)  

#* detail by ENTREZ 
#* @param entrez:int id del gen
#* @get /
#* @tag detail
#* @serializer unboxedJSON 
function(entrez) {

    details <- tryCatch({
        AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("GENETYPE", "GENENAME", "SYMBOL", "ALIAS"), keytype = "ENTREZID")
    }, error = function(e) NULL)

    if (is.null(detail)) {
        response <- list(
            code = 500,
            message = "try catch",
            data = NULL
        )
    } else if (length(detail) == 0) {
        response <- list(
            code = 404,
            message = paste("no se encontraron valores para el entrez ", alias),
            data = NULL
        )
    } else {
      response <- list(
        code = 200,
        message = "Ok"
        data = list(
            entrez = entrez,
            symbol = unique(details$SYMBOL),
            genename = unique(details$GENENAME),
            genetype = unique(details$GENETYPE),
            alias = unique(details$ALIAS)
        )
      )
    }
    return(response)
}
