library(AnnotationDbi)
library(org.Hs.eg.db)  

#* detail by ENTREZ 
#* @param entrez:int id del gen
#* @get /
#* @tag endpoints
#* @serializer unboxedJSON 
function(entrez) {

  start_time <- Sys.time()

  details <- tryCatch({
    AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("GENETYPE", "SYMBOL", "ALIAS"), keytype = "ENTREZID")
  }, error = function(e) NULL)

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))


  if (is.null(detail)) {
      result <- list(
          code = 500,
          datetime = start_time,
          time_secs = time,
          data = list (
              message = "try catch"
          )
      )
  } else if (length(detail) == 0) {
      result <- list(
          code = 404,
          datetime = start_time,
          time_secs = time,
          data = list (
              message = paste("no se encontraron valores para el entrez ", alias)
          )
      )
  } else {
    result <- list(
      code = 200,
      datetime = start_time,
      time_secs = time,
      data = list(
        entrez = entrez,
        symbol = unique(details$SYMBOL),
        genetype = unique(details$GENETYPE),
        alias = unique(details$ALIAS)
      )
    )
  }
  return(result)
}
