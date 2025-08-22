#* entrez devuelve el entrez del value
#* @param value es un symbol o alias
#* @get /
#* @serializer unboxedJSON 
function(value = "DHCR7") {

    start_time <- Sys.time()
    label <- NULL
    entrez <- tryCatch({
        AnnotationDbi::select(org.Hs.eg.db, keys = value, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID
    }, error = function(e) NULL)
    

    if (length(entrez) != 0) {
        label <- "SYMBOL" 
        status <- "success"
    } else {
        entrez <- tryCatch({
            AnnotationDbi::select(org.Hs.eg.db, keys = value, columns = "ENTREZID", keytype = "ALIAS")$ENTREZID
        }, error = function(e) NULL)
    
        status <- "error"

        if (length(entrez) != 0) {
            label <- "ALIAS" 
            status <- "success"
        }
    }

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    result <- list(
        status = status, 
        time_secs = time,
        data = list(
            value = value,
            label = label,
            entrez = entrez
        )
    )
}