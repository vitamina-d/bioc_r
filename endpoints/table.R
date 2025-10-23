library(AnnotationDbi)
library(org.Hs.eg.db)

#* table json
#* @get /
#* @tag tables
#* @serializer unboxedJSON 
function() {

    #table <- as.list(org.Hs.egSYMBOL)
    #list <- as.list(org.Hs.egGENENAME)

    select <- AnnotationDbi::select(
        org.Hs.eg.db,
        keys = keys(org.Hs.eg.db), #keytype = "ENTREZID"),
        columns = c("SYMBOL", "GENENAME"),
        keytype = "ENTREZID"
    )

    return(list(
        code = 200,
        message = "Ok.",
        data = list (
            count = nrow(select),
            table = select
            )
        )
    )
}
