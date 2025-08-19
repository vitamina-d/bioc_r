library(Biostrings)

#* percent muestra el porcentaje de bases, A/T y C/G de una secuencia
#* @param seq Secuencia
#* @get /
#* @tag endpoints
#* @serializer unboxedJSON 
function(seq = "ACGT") {

  start_time <- Sys.time()

  # biostrings

  DNA_str <- DNAString(seq)
  
  counter_base <- alphabetFrequency(DNA_str, baseOnly = TRUE)

  AT <- counter_base["A"] + counter_base["T"]
  CG <- counter_base["C"] + counter_base["G"] #termicamente mas estables

  ATCG <- AT + CG
  
  AT_percent <- (AT / ATCG) * 100
  CG_percent <- (CG / ATCG) * 100


  ####
  pattern_CpG <- "CG"
  counter_CpG <- countPattern(pattern_CpG, DNA_str)
  match_CpG <- matchPattern(pattern_CpG, DNA_str)

  CpG_ranges <- as.data.frame(match_CpG@ranges)

  ######## Convertir cada fila a lista: [{start:.., end:.., width:..}, ...]
  CpG_ranges_list <- apply(CpG_ranges, 1, function(row) {
    list(
      start = as.integer(row["start"]),
      end = as.integer(row["end"]),
      width = as.integer(row["width"])
    )
  })
   
  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  result <- list(
    status = "success", 
    time_secs = time,
    data = list(
      nucleotides = list(
        labels = c("A", "T", "C", "G"),
        counts = c(
          as.integer(counter_base["A"]),
          as.integer(counter_base["T"]),
          as.integer(counter_base["C"]),
          as.integer(counter_base["G"])
        )
      ),
      A = counter_base["A"],
      T = counter_base["T"],
      C = counter_base["C"],
      G = counter_base["G"],
      ATCG = ATCG,
      AT_percent = AT_percent,
      CG_percent = CG_percent,
      counter_CpG = counter_CpG,
      CpG_ranges = CpG_ranges
    )
  )
}