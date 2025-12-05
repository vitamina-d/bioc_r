library(Biostrings)

#* align devuelve el alineamiento global o local de dos secuecias
#* @param pattern Lectura
#* @param subject Genoma de referencia
#* @param type "global", "local", "overlap"
#* @param gapOpening:number
#* @param gapExtension:number
#* @post /
#* @tag sequence
#* @serializer unboxedJSON 
function(pattern, subject, type, gapOpening, gapExtension, res) { 
    gapOpening <- as.numeric(gapOpening)
    gapExtension <- as.numeric(gapExtension)

    if (is.null(pattern) || pattern == "") {
        res$status <- 400
        stop("Ingrese pattern", call. = FALSE) 
    }
    if (is.null(subject) || subject == "") {
        res$status <- 400
        stop("Ingrese subject", call. = FALSE) 
    }
    if (!(type %in% c("global","local","overlap"))) {
        res$status <- 400
        stop("Ingrese type valido: 'global' | 'local' | 'overlap'", call. = FALSE) 
    }
    if (is.na(gapOpening) || is.na(gapExtension) || gapOpening < 0 || gapExtension < 0) {
        res$status <- 400
        stop("gapOpening y gapExtension deben ser numeros positivos", call. = FALSE)
    }

    tryCatch({
        seqA <- DNAString(pattern)
        seqB <- DNAString(subject)
        align <- pairwiseAlignment(	seqA, seqB, type = type, substitutionMatrix = NULL, gapOpening = gapOpening, gapExtension = gapExtension)

        list(
            code = 200,
            message = "Ok.",
            data = list(
                score = score(align),
                pattern_align = as.character(pattern(align)),
                subject_align = as.character(subject(align))
                )
            )
    }, error = function(e){
        # exception
        res$status <- 500
        stop(paste("Error de servicio R: ", e$message), call. = FALSE)
    })
}