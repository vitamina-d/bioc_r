library(Biostrings)

#* align devuelve el alineamiento global o local de dos secuecias
#* @param pattern Lectura
#* @param subject Genoma de referencia
#* @param type "global", "local", "overlap"
#* @param match number
#* @param mismatch number
#* @param gapOpening number
#* @param gapExtension number
#* @post /
#* @tag sequence
#* @serializer unboxedJSON 
function(pattern, subject, type, match, mismatch, gapOpening, gapExtension) { ## allocate 733.5 Mb of memory (i.e. length(pattern) * length(subject) * 3 bytes)
    start_time <- Sys.time()

    seqA <- DNAString(pattern)
    seqB <- DNAString(subject)

    matrix <- nucleotideSubstitutionMatrix(match=match, mismatch=mismatch, baseOnly=TRUE)
    align <- pairwiseAlignment(	seqA, seqB, type = type, substitutionMatrix = matrix, gapOpening = gapOpening, gapExtension = gapExtension)

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    response <- list(
        code = 200,
        datetime = start_time,
        time_secs = time,
        data = list(
          message = "Ok",
          score = score(align),
          pattern_align = as.character(pattern(align)),
          subject_align = as.character(subject(align))
        )
    )
}
