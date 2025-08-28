library(AnnotationDbi)
library(org.Hs.eg.db)  

#* isentrez valida que el valor sea un ENTREZ
#* @param entrez
#* @get /
#* @tag entrez
#* @serializer unboxedJSON 
function(entrez) {

    start_time <- Sys.time()
    is_entrez = FALSE
    
    if (is.null(entrez) || entrez == "" ) {
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

    ids <- tryCatch({
        keys(org.Hs.eg.db, keytype = "ENTREZID")    
    }, error = function(e) NULL)

    if (is.null(ids)) {
        result <- list(
            code = 500,
            message = "Error del servidor.",
            datetime = start_time,
            time_secs = time,
            data = NULL
        )
    }

    if (entrez %in% ids) {
        is_entrez = TRUE
    }

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    result <- list(
        code = 200,
        message = "Ok.",
        datetime = start_time,
        time_secs = time,
        data = list(
            is_entrez = is_entrez
        )
    )
}
