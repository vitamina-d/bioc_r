library(Biostrings)

#* complement devuelve la secuencia complementaria e inversa
#* @param seq Secuencia
#* @param is_reverse:boolean 
#* @param is_complement:boolean 
#* @get /
#* @tag endpoints
#* @serializer unboxedJSON 
function(seq, is_reverse, is_complement) {

    start_time <- Sys.time()
sequence <- seq
    DNA_str <- DNAString(seq)

    if (is_reverse) {
        sequence <- reverse(DNA_str)
    } else if (is_complement) {
        sequence <- complement(DNA_str)
    } 

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    result <- list(
        code = 200,
        datetime = start_time,
        time_secs = time,
        data = list(
            sequence = as.character(sequence) # AAString a texto
        )
    )

    return(result)
}
