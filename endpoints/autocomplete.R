library(org.Hs.eg.db)
library(AnnotationDbi)

#* entrez devuelve el entrez si lo encuentra
#* @param input 
#* @get /
#* @tag detail
#* @serializer unboxedJSON 
function(input) {

    start_time <- Sys.time()

    alias <- unique(keys(org.Hs.eg.db, keytype = "ALIAS"))

    matches <- grep(input, alias, ignore.case = TRUE, value = TRUE)

    if (length(matches) == 0)  {
        matches <- NULL
    } 
    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    if (length(matches) == 0) {
        result <- list(
            code = 404,
            message = paste("no se encontro"),
            datetime = start_time,
            time_secs = time,
            data = list()
        )
    } else {
        result <- list(
            code = 200,
            message = "Ok",
            datetime = start_time,
            time_secs = time,
            data = head(matches, 20)
        )
    }
    return(result)
}