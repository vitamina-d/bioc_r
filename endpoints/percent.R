library(Biostrings)

#* percent muestra el porcentaje de bases, A/T y C/G de una secuencia
#* @post /
#* @parser json
#* @param seq
#* @tag sequence
#* @serializer unboxedJSON 
function(seq, res) {

    if (is.null(seq) || seq == "") {
        res$status <- 400
        stop("Ingrese un valor.", call. = FALSE)
    }
  
    DNA_str <- tryCatch({
        DNAString(seq)
    }, error = function(e){
        res$status <- 400
        stop(paste("No se puede convertir a DNAString:", e$message), call. = FALSE)
    })
    
    counter_base <- alphabetFrequency(DNA_str, baseOnly = TRUE)
    pattern_CpG <- DNAString("CG")
    counter_CpG <- countPattern(pattern_CpG, DNA_str)
    match_CpG <- matchPattern(pattern_CpG, DNA_str)
    cpg_info <- as.list(match_CpG@ranges@start)

    result <- list(
        code = 200,
        message = "Ok.",
        data = list(
            sequence_length = nchar(DNA_str),
            nucleotides = as.list(counter_base),
            cpg_islands = list(
                count = counter_CpG,
                start = cpg_info
            ),
            sequence = as.character(DNA_str)
        )
    )
}
