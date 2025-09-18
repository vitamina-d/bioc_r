library(Biostrings)

#* translate devuelve el marco pedido 
#* @param sequence Secuencia
#* @param frame:number marco +1 +2 +3 -1 -2 -3
#* @post /
#* @tag sequence
#* @serializer unboxedJSON 
function(sequence, frame) {
frame <- as.numeric(frame)
    print(frame)
    print(sequence)

    start_time <- Sys.time()
    DNA_str <- DNAString(sequence)

    if (frame < 0) {
        DNA_str <- reverseComplement(DNA_str)
    } 
    frame <- abs(frame)

    amino <- translate(subseq(DNA_str, start=frame))

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    result <- list(
        code = 200,
        message = "Ok",
        datetime = start_time,
        time_secs = time,
        data = list(
            sequence = as.character(amino)
        )
    )

    return(result)
}
