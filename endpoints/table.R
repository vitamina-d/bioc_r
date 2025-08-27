library(AnnotationDbi)
library(org.Hs.eg.db)

#* table json
#* @get /
#* @serializer unboxedJSON 
function() {

    start_time <- Sys.time()

    #table <- as.list(org.Hs.egSYMBOL)

    select <- AnnotationDbi::select(
        org.Hs.eg.db,
        keys = keys(org.Hs.eg.db), #keytype = "ENTREZID"),
        columns = c("SYMBOL", "GENENAME"),
        keytype = "ENTREZID"
    )

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    return(list(
        code = 200,
        datetime = start_time,
        time_secs = time,
        data = list (
            count = nrow(select),
            table = select
            )
        )
    )
}
