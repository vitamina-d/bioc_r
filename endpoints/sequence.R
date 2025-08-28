library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

library(AnnotationDbi)
library(org.Hs.eg.db)  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#* sequence devuelve la secuencia completa o de exones a partir de su entrez
#* @param entrez EntrezID
#* @param complete:boolean Secuencia completa (TRUE) o solo exones (FALSE)
#* @get /
#* @tag sequence
#* @serializer unboxedJSON 
function(entrez, complete = TRUE) {

    start_time <- Sys.time()

    human_genome <- BSgenome.Hsapiens.UCSC.hg38
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    #coordenadas: objeto GRanges
    #secuencia: objeto DNAStringSet de biostrings
    if(complete){
      
        coord_gene <- genes(txdb)[entrez]
        sequence <- getSeq(human_genome, coord_gene) #devuelve la codificante

    } else {
        coord_exones <- exonsBy(txdb, by = "gene")[[entrez]] 
        seq_exones <- getSeq(human_genome, coord_exones)

        #exones
        sequence <- do.call(xscat, as.list(seq_exones))
    }

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    result <- list(
        code = 200,
        datetime = start_time,
        time_secs = time,
        data = list(
            message = "Ok.",
            complete = as.logical(complete),
            sequence_length = nchar(sequence),
            sequence = as.character(sequence)
        )
    )
}
