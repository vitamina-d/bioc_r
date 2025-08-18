library(Biostrings)

#* pairwise devuelve el alineamiento global o local de dos secuecias
#* @param pattern Lectura
#* @param subject Genoma de referencia
#* @param global:boolean Alineamiento global (TRUE) o local (FALSE)
#* @get /
#* @tag endpoints
#* @serializer unboxedJSON 
function(pattern = "", subject = "", global = TRUE) {
  start_time <- Sys.time()

  seqA <- DNAString(pattern)
  seqB <- DNAString(subject)

  type <- ifelse(global, "global", "local")

  # Alineamiento global (tipo Needleman-Wunsch)
  align <- pairwiseAlignment(	seqA, seqB, substitutionMatrix = NULL, gapOpening = -2, gapExtension = -1, type = type)

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  result <- list(
    status = "success", 
    time_secs = time,
    data = list(
      score = score(align),
      type = type,
      pattern = pattern,
      subject = subject,
      pattern_align = as.character(pattern(align)),
      subject_align = as.character(subject(align))
    )
  )
}
