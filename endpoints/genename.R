library(AnnotationDbi)
library(org.Hs.eg.db)

#* genename json
#* @get /
#* @tag tables
#* @serializer unboxedJSON 
function() {

    start_time <- Sys.time()

    list <- as.list(org.Hs.egGENENAME)

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    return(list(
        code = 200,
        datetime = start_time,
        time_secs = time,
        data = list (
            message = "Ok.",
            count = nrow(list),
            table = list
            )
        )
    )
}
