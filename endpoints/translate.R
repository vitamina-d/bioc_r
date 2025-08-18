library(plumber)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

library(AnnotationDbi)
library(org.Hs.eg.db)  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#* Translate traduce la secuencia a proteina
#* @param seq Secuencia
#* @get /translate
#* @tag endpoints
#* @serializer unboxedJSON 
function(seq = "ACGT") {

  start_time <- Sys.time()

  # biostrings
  DNA_str <- DNAString(seq)
  protein <- translate(dna)

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  result <- list(
    status = "success", 
    time_secs = time,
    data = protein
  )  
}
