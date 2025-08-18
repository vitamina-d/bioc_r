library(Biostrings)

#* translate traduce la secuencia a proteina
#* @param seq Secuencia
#* @get /
#* @tag endpoints
#* @serializer unboxedJSON 
function(seq = "ACGT") {

  start_time <- Sys.time()

  # biostrings
  DNA_str <- DNAString(seq)
  protein <- translate(DNA_str)

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  result <- list(
    status = "success", 
    time_secs = time,
    data = as.character(protein) # AAString a texto
  )  

  return(result)
}
