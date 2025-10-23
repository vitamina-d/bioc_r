library(AnnotationDbi)
library(org.Hs.eg.db)

#* genename json
#* @get /
#* @tag tables
#* @serializer unboxedJSON 
function() {

    list <- org.Hs.egGENENAME

    return(list(
        code = 200,
        data = list (
            message = "Ok.",
            count = nrow(list),
            table = list
            )
        )
    )
}
