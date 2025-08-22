library(AnnotationDbi)
library(org.Hs.eg.db)  

#* isentrez valida que el valor sea un ENTREZ
#* @param num entrez
#* @get /
#* @tag endpoints
#* @serializer unboxedJSON 
function(entrez = "1717") {

  start_time <- Sys.time()
  
  is_entrez = FALSE

  #if (entrez %in% keys(org.Hs.eg.db, keytype = "ENTREZID")) {
   # is_entrez = TRUE
  #}

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))


  result <- list(
    status = "success", 
    time_secs = time,
    data = list(
      is_entrez = is_entrez,
    )
  )
}