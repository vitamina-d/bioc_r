library(plumber)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationDbi)
library(org.Hs.eg.db)  
library(GenomicFeatures)   # genes(), exonsBy()
library(Biostrings)        # getSeq(), xscat()

#* @get /sequence
#* @param gene_symbol Nombre del gen
#* @param complete:boolean Solo exones (CDS):FALSE
#* @serializer unboxedJSON 
function(gene_symbol="DHCR7", complete = TRUE) {
  
  # inicio contador
  start_time <- Sys.time()

  #genoma y coord
  human_genome <- BSgenome.Hsapiens.UCSC.hg38
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

  ### identificar el gen 
  #mapIds(org.Hs.eg.db, keys = gene_entrez, column = "SYMBOL", keytype = "ENTREZID")
  #mapIds(org.Hs.eg.db, keys = gene_symbol, column = "CHR", keytype = "SYMBOL")
  gene_entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbol, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID

  ####VALIDAR

  #coordenadas: objeto GRanges
  #secuencia: objeto DNAStringSet de biostrings
  if(complete){
    
    coord_gene <- genes(txdb)[gene_entrez]
    sequence <- getSeq(human_genome, coord_gene)

  } else {
      coord_exones <- exonsBy(txdb, by = "gene")[[gene_entrez]] 
      seq_exones <- getSeq(human_genome, coord_exones)

      #CDS: secuencia de la region codificante
      sequence <- do.call(xscat, as.list(seq_exones))
  }

  type <- ifelse(complete, "complete", "cds")

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))
 
  list(
    time_secs = time,
    type = type,
    sequence_length = nchar(sequence),
    sequence = as.character(sequence)
  )
}