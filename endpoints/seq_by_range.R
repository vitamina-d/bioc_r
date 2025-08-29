library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

#* seq_by_range devuelve la secuencia dado el cromosoma y el rango
#* @param chrom Cromosoma (ej: chr11)
#* @param start Inicio
#* @param end Fin
#* @get /
#* @tag sequence
#* @serializer unboxedJSON 
function(chrom, start, end) {
    start_time <- Sys.time()

    start <- as.integer(start)
    end <- as.integer(end)

    seq <- BSgenome.Hsapiens.UCSC.hg38[[chrom]][start:end]

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    result <- list(
        code = 200,
        message = "Ok.",
        datetime = start_time,
        time_secs = time,
        data = list(
            sequence_length = nchar(sequence),
            sequence = as.character(seq),
            complete = TRUE
        )
    )
}
