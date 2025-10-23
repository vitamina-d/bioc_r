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
    
    seqA <- DNAString(pattern)
    seqB <- DNAString(subject)

    #matrix <- nucleotideSubstitutionMatrix(match = match, mismatch = mismatch, baseOnly=TRUE)
    align <- pairwiseAlignment(	seqA, seqB, type = type, substitutionMatrix = NULL, gapOpening = gapOpening, gapExtension = gapExtension)

    response <- list(
        code = 200,
        message = "Ok",
        data = list(
          score = score(align),
          pattern_align = as.character(pattern(align)),
          subject_align = as.character(subject(align))
        )
    )
}
