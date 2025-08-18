library(AnnotationDbi)
library(org.Hs.eg.db)  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#* detail muestra info de un gen, dado su symbol
#* @param symbol Nombre del gen
#* @get /
#* @tag endpoints
#* @serializer unboxedJSON 
function(symbol = "DHCR7") {

  start_time <- Sys.time()

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

  entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = symbol, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID
  details <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("ENSEMBL", "ENSEMBLPROT", "UNIPROT", "ENTREZID", "GENETYPE", "MAP", "SYMBOL"), keytype = "ENTREZID")

  # Obtener rangos, sin filtrar genes que están en ambas cadenas .. granges / listgranges
  #range_list <- genes(txdb, single.strand.genes.only = FALSE)[entrez]
  range <- genes(txdb)[entrez]
  range_df <- as.data.frame(range)

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))


  result <- list(
    status = "success", 
    time_secs = time,
    data = list(
      entrezID = entrez,
      symbol = symbol,
      type = unique(details$GENETYPE),
      location_chr = unique(details$MAP),
      chr = as.character(range_df$seqnames),
      start = range_df$start,
      end = range_df$end,
      length = range_df$width,
      strand = as.character(range_df$strand),
      ensembl_id_gene = unique(details$ENSEMBL),
      ensembl_id_protein = unique(details$ENSEMBLPROT),
      uniprot_id = unique(details$UNIPROT)
    )
  )
}