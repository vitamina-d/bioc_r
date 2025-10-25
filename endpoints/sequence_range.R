library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

#* seq_by_range devuelve la secuencia dado el cromosoma y el rango
#* @param chrom Cromosoma (ej: chr11)
#* @param start Inicio
#* @param end Fin
#* @get /
#* @tag sequence
#* @serializer unboxedJSON 
function(chrom, start, end, res) {
    tryCatch({
        if (is.null(chrom) || chrom == "" || is.null(start) || is.null(end)) {
            res$status <- 400
            stop("Parametros invalidos.", call. = FALSE)
        }

        start <- as.integer(start)
        end <- as.integer(end)

        if (!(chrom %in% names(BSgenome.Hsapiens.UCSC.hg38))) {
            res$status <- 404
            stop(paste(chrom, " no es un cromosoma valido."), call. = FALSE)
        }
        chr_len <- length(BSgenome.Hsapiens.UCSC.hg38[[chrom]])
        if (start <= 0 || end <= 0 || start > chr_len || end > chr_len || start > end) {
            res$status <- 400
            stop("Ingrese un rango valido.", call. = FALSE)
        }
        
        seq <- BSgenome.Hsapiens.UCSC.hg38[[chrom]][start:end]

        list(
            code = 200,
            message = "Ok.",
            data = list(
                sequence_length = nchar(seq),
                sequence = as.character(seq),
                complete = TRUE
            )
        )
    }, error = function(e) {
        res$status <- 500
        stop(paste("Error de servicio R: ", e$message), call. = FALSE)
    })
}