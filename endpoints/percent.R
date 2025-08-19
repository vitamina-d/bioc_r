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

  total <- AT + CG
  
  at_percent <- (AT / total) * 100
  cg_percent <- (CG / total) * 100


  ####
  pattern_CpG <- "CG"
  counter_CpG <- countPattern(pattern_CpG, DNA_str)
  match_CpG <- matchPattern(pattern_CpG, DNA_str)

  cpg_info <- as.data.frame(match_CpG@ranges)
  
  ######## Convertir cada fila a lista: [{start:.., end:.., width:..}, ...]
  cpg_info_list <- apply(cpg_info, 1, function(row) {
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
      composition = list(
        nucleotides = list(
          labels = c("A", "T", "C", "G"),
          counts = c(
            as.integer(counter_base["A"]),
            as.integer(counter_base["T"]),
            as.integer(counter_base["C"]),
            as.integer(counter_base["G"])
          )
        ),
        total = total,
        at_percent = at_percent,
        cg_percent = cg_percent
      ),
      cpg_islands = list(
        length = counter_CpG,
        ranges = cpg_info
      )
    )
  )
}

