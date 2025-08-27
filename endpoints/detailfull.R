library(AnnotationDbi)
library(org.Hs.eg.db)  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#* detail muestra info de un gen, dado su ENTREZ 
#* @param entrez id del gen
#* @get /
#* @tag endpoints
#* @serializer unboxedJSON 
function(entrez = "1717") {

  start_time <- Sys.time()

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

  details <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("ENSEMBL", "ENSEMBLPROT", "UNIPROT", "ENTREZID", "GENETYPE", "MAP", "SYMBOL", "GO", "ONTOLOGY"), keytype = "ENTREZID")

  range <- genes(txdb)[entrez]
  range_df <- as.data.frame(range)

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  response <- list(
    code = 200,
    datetime = start_time,
    time_secs = time,
    data = list(
      entrezID = entrez,
      symbol = unique(details$SYMBOL),
      type = unique(details$GENETYPE),
      location = list(
        citogenetic = unique(details$MAP),
        strand = as.character(range_df$strand),
        chr = as.character(range_df$seqnames),
        start = range_df$start,
        end = range_df$end,
        length = range_df$width
      ),
      ensembl_id_gene = unique(details$ENSEMBL),
      ensembl_id_protein = unique(details$ENSEMBLPROT),
      uniprot_ids = unique(details$UNIPROT)
      #go = unique(details$go),
      #ontology = unique(details$ontology)
    )
  )
}