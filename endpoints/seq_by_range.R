library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

#* seq_by_range devuelve la secuencia dado el cromosoma y el rango
#* @param chrom Cromosoma (ej: "dhcr7")
#* @param start Inicio
#* @param end Fin
#* @get /
#* @serializer unboxedJSON 
function(chrom = "chr11", start = 71428193, end = 71452868) {
  start_time <- Sys.time()

  start <- as.integer(start)
  end <- as.integer(end)

  seq <- BSgenome.Hsapiens.UCSC.hg38[[chrom]][start:end]
  sequence <- as.character(seq)
  sequence_length <- nchar(sequence)
  
  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  result <- list(
    status = "success", 
    time_secs = time,
    data = list(
      sequence_length = sequence_length,
      sequence = sequence,
      complete = TRUE
    )
  )
}
