library(Biostrings)

#* complement devuelve la secuencia complementaria e inversa
#* @param seq Secuencia
#* @param to_reverse:boolean 
#* @param to_complement:boolean 
#* @post /
#* @tag sequence
#* @serializer unboxedJSON 
function(seq, to_reverse, to_complement) {

    start_time <- Sys.time()
    DNA_str <- DNAString(seq)

    if (to_reverse) {
        DNA_str <- reverse(DNA_str)
    } 
    
    if (to_complement) {
        DNA_str <- complement(DNA_str)
    } 

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    result <- list(
        code = 200,
        message = "Ok",
        datetime = start_time,
        time_secs = time,
        data = list(
            sequence = as.character(DNA_str) # AAString a texto
        )
    )

    return(result)
}
