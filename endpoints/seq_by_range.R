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
    start <- as.integer(start)
    end <- as.integer(end)

    seq <- BSgenome.Hsapiens.UCSC.hg38[[chrom]][start:end]

    result <- list(
        code = 200,
        message = "Ok.",
        data = list(
            sequence_length = nchar(seq),
            sequence = as.character(seq),
            complete = TRUE
        )
    )
}
