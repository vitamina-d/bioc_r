library(Biostrings)

#* percent muestra el porcentaje de bases, A/T y C/G de una secuencia
#* @post /
#* @parser json
#* @param seq
#* @tag sequence
#* @serializer unboxedJSON 
function(seq) {

    if (is.null(seq) || seq == "" ) {
        result <- list(
            code = 400,
            message = "Ingrese un valor.",
            data = NULL
        )
        return(result)
    }  
  
    DNA_str <- tryCatch({
        DNAString(seq)
        #as.character(DNA_str) seq
    }, error = function(e) NULL)
      
    if (is.null(DNA_str)) {
        result <- list(
            code = 400,
            message = "Ingrese una secuencia valida.",
            data = NULL
        )
        return(response)

    } else {
        counter_base <- alphabetFrequency(DNA_str, baseOnly = TRUE)
        pattern_CpG <- "CG"
        counter_CpG <- countPattern(pattern_CpG, DNA_str)
        match_CpG <- matchPattern(pattern_CpG, DNA_str)
        cpg_info <- as.list(match_CpG@ranges@start)

        result <- list(
            code = 200,
            message = "Ok.",
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
