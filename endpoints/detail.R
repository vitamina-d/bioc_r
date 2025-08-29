library(AnnotationDbi)
library(org.Hs.eg.db)  

#* detail by ENTREZ 
#* @param entrez:int id del gen
#* @get /
#* @tag detail
#* @serializer unboxedJSON 
function(entrez) {

    start_time <- Sys.time()

    details <- tryCatch({
        AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("GENETYPE", "GENENAME", "SYMBOL", "ALIAS"), keytype = "ENTREZID")
    }, error = function(e) NULL)

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))


    if (is.null(detail)) {
        response <- list(
            code = 500,
            message = "try catch",
            datetime = start_time,
            time_secs = time,
            data = NULL
        )
    } else if (length(detail) == 0) {
        response <- list(
            code = 404,
            message = paste("no se encontraron valores para el entrez ", alias),
            datetime = start_time,
            time_secs = time,
            data = NULL
        )
    } else {
      response <- list(
        code = 200,
        message = "Ok",
        datetime = start_time,
        time_secs = time,
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
