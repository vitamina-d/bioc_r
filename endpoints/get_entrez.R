#* get_entrez back the input
#* @param msg The message to echo
#* @get /
#* @serializer unboxedJSON 
function(msg = "msg") {

  start_time <- Sys.time()
  
  entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = msg, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID

  if (length(entrez) == 0) {
    entrez <- AnnotationDbi::select( org.Hs.eg.db, keys = msg, columns = "ENTREZID", keytype = "ALIAS")$ENTREZID
  }

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))


  result <- list(
    status = "success", 
    time_secs = time,
    data = list(
      entrez = entrez
    )
  )
}