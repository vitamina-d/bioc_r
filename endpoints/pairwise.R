library(plumber)
library(Biostrings)

#* @get /paiswise
#* @param pattern Lectura
#* @param subject Genoma de referencia
#* @param global:boolean Alineamiento global (Needleman-Wunsch) o local (Smith-Waterman)
#* @serializer unboxedJSON 
function(pattern = "", subject = "", global = TRUE) {
  # inicio contador
  start_time <- Sys.time()

  seqA <- DNAString(pattern)
  seqB <- DNAString(subject)

  type <- ifelse(global, "global", "local")

  # Alineamiento global (tipo Needleman-Wunsch)
  align <- pairwiseAlignment(	seqA, seqB, substitutionMatrix = NULL, gapOpening = -2, gapExtension = -1, type = type)

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  list(
    time_secs = time,
    score = score(align), # ver matriz
    type = type,
    pattern = pattern,
    subject = subject,
    pattern_align = as.character(pattern(align)),
    subject_align = as.character(subject(align))
  )
}