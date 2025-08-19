library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

library(AnnotationDbi)
library(org.Hs.eg.db)  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#* seq_by_symbol devuelve la secuencia completa o de exones, dado el symbol de un gen
#* @param gene_symbol Nombre del gen
#* @param complete:boolean Secuencia completa (TRUE) o solo exones (FALSE)
#* @get /
#* @tag endpoints
#* @serializer unboxedJSON 
function(gene_symbol="DHCR7", complete = TRUE) {

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

  type <- ifelse(complete, "complete", "exons")

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  result <- list(
    status = "success", 
    time_secs = time,
    data = list(
      complete = complete,
      sequence_length = nchar(sequence),
      sequence = as.character(sequence)
    )
  )
}
