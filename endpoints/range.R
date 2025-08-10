library(plumber)
library(BSgenome.Hsapiens.UCSC.hg38)

#* Devuelve la secuencia de BSGenome, dado un cromosoma y los valores inicio y fin.

#* @get /range
#* @param chrom Cromosoma (ej: "dhcr7")
#* @param start Inicio
#* @param end Fin
#* @serializer unboxedJSON 
function(chrom = "chr11", start = 100000, end = 100100) {
  # inicio contador
  start_time <- Sys.time()

  start <- as.integer(start)
  end <- as.integer(end)

  seq <- BSgenome.Hsapiens.UCSC.hg38[[chrom]][start:end]
  sequence <- as.character(seq)
  sequence_length <- nchar(sequence)
  
  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  list(
    time_secs = time,
    chrom = chrom,
    start = start,
    end = end,
    sequence_length = sequence_length,
    sequence = sequence
  )
}
