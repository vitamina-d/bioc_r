# library(AnnotationDbi)
# library(org.Hs.eg.db)  

# #* entrez_num retorna el ENTREZ de un SYMBOL o ALIAS
# #* @param input SYMBOL o ALIAS
# #* @get /
# #* @tag endpoints
# #* @serializer unboxedJSON 
# function(input = "DHCR7") {

#   start_time <- Sys.time()
  
#   #entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = input, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID
# entrez  <- "entrez"
#   #if (length(entrez) == 0) {
#    # entrez <- AnnotationDbi::select( org.Hs.eg.db, keys = input, columns = "ENTREZID", keytype = "ALIAS")$ENTREZID
#   #}

#   end_time <- Sys.time()
#   time <- as.numeric(difftime(end_time, start_time, units = "secs"))


#   result <- list(
#     status = "success", 
#     time_secs = time,
#     data = list(
#       entrez = entrez
#     )
#   )
# }