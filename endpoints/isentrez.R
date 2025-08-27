library(AnnotationDbi)
library(org.Hs.eg.db)  

#* isentrez valida que el valor sea un ENTREZ
#* @param entrez
#* @get /
#* @tag endpoints
#* @serializer unboxedJSON 
function(entrez) {

    start_time <- Sys.time()
    is_entrez = FALSE
    
    if (is.null(entrez) || entrez == "" ) {
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

    ids <- tryCatch({
        keys(org.Hs.eg.db, keytype = "ENTREZID")    
    }, error = function(e) NULL)

    if (is.null(ids)) {
        result <- list(
            code = 500,
            datetime = start_time,
            time_secs = time,
            data = list (
                message = "Error del servidor."
            )
        )
    }

    if (entrez %in% ids) {
        is_entrez = TRUE
    }

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    result <- list(
        code = 200,
        datetime = start_time,
        time_secs = time,
        data = list(
            is_entrez = is_entrez
        )
    )
}
