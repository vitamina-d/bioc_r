library(Biostrings)

#* percent muestra el porcentaje de bases, A/T y C/G de una secuencia
#* @post /
#* @parser json
#* @param seq
#* @tag sequence
#* @serializer unboxedJSON 
function(seq) {

    start_time <- Sys.time()

    if (is.null(seq) || seq == "" ) {
        end_time <- Sys.time()
        time <- as.numeric(difftime(end_time, start_time, units = "secs"))

        result <- list(
            code = 400,
            message = "Ingrese un valor.",
            datetime = start_time,
            time_secs = time,
            data = NULL
        )
        return(result)
    }  
  
    DNA_str <- tryCatch({
        DNAString(seq)
        #as.character(DNA_str) seq
    }, error = function(e) NULL)
      
    if (is.null(DNA_str)) {
        end_time <- Sys.time()
        time <- as.numeric(difftime(end_time, start_time, units = "secs"))

        result <- list(
            code = 400,
            message = "Ingrese una secuencia valida.",
            datetime = start_time,
            time_secs = time,
            data = NULL
        )
        return(response)

    } else {
        counter_base <- alphabetFrequency(DNA_str, baseOnly = TRUE)
        pattern_CpG <- "CG"
        counter_CpG <- countPattern(pattern_CpG, DNA_str)
        match_CpG <- matchPattern(pattern_CpG, DNA_str)
        cpg_info <- as.list(match_CpG@ranges@start)

        end_time <- Sys.time()
        time <- as.numeric(difftime(end_time, start_time, units = "secs"))

        result <- list(
            code = 200,
            message = "Ok.",
            datetime = start_time,
            time_secs = time,
            data = list(
                composition = list(
                    length = sum(counter_base),
                    nucleotides = as.list(counter_base)
                ),
                cpg_islands = list(
                    count = counter_CpG,
                    start = cpg_info
                )
            )
        )
    }
    return(result)
}


        #nucleotides = list(
         # labels = c("A", "T", "C", "G"),
          #counts = c(
           # as.integer(counter_base["A"]),
            #as.integer(counter_base["T"]),
            #as.integer(counter_base["C"]),
            #as.integer(counter_base["G"])
          #)
        #),

        #dataframe
        # "start": [
        # {
        #   "match_CpG@ranges@start": 11
        # },
        # {
        #   "match_CpG@ranges@start": 13
        # },
        # {
        #   "match_CpG@ranges@start": 20
        # },
        # {
        #list