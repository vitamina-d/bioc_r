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
function(pattern, subject, type, gapOpening, gapExtension) { ## allocate 733.5 Mb of memory (i.e. length(pattern) * length(subject) * 3 bytes)
    gapOpening <- as.numeric(gapOpening)
    gapExtension <- as.numeric(gapExtension)
    # Validación de entrada
    if (is.null(pattern) || pattern == "") {
        stop("Ingrese pattern", call. = FALSE) # HTTP 400
    }
    if (is.null(subject) || subject == "") {
        stop("Ingrese subject", call. = FALSE) # HTTP 400
    }
    if (!(type %in% c("global","local","overlap"))) {
        stop("Ingrese type valido: 'global' | 'local' | 'overlap'", call. = FALSE) # HTTP 400
    }
    tryCatch({
        seqA <- DNAString(pattern)
        seqB <- DNAString(subject)
        #matrix <- nucleotideSubstitutionMatrix(match = match, mismatch = mismatch, baseOnly=TRUE)
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
        stop(paste("Internal server error:", e$message), call. = FALSE)
    })
}