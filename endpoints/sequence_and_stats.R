library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

library(AnnotationDbi)
library(org.Hs.eg.db)  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#* sequence_and_stats devuelve la secuencia completa o de exones a partir de su entrez y percent
#* @param entrez EntrezID
#* @param complete:boolean Secuencia completa (TRUE) o solo exones (FALSE)
#* @get /
#* @tag sequence
#* @serializer unboxedJSON 
function(entrez, complete = TRUE) {

    human_genome <- BSgenome.Hsapiens.UCSC.hg38
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    #coordenadas: objeto GRanges
    #secuencia: objeto DNAStringSet de biostrings
    if(complete){
        #coord_gene <- genes(txdb)[entrez]
        
        coord_gene <- tryCatch({
            genes(txdb)[entrez]
        }, error = function(e) NULL)

        if (is.null(coord_gene)) {
            result <- list(
                code = 404,
                message = paste("no se encontro secuencia para el entrez: ", entrez),
                data = NULL
            )
            return(result)
        } 

        sequence <- getSeq(human_genome, coord_gene) 
        DNA_str <- unlist(sequence)

    } else {

        #coord_exones <- exonsBy(txdb, by = "gene")[[entrez]] 
        #seq_exones <- getSeq(human_genome, coord_exones)
        #DNA_str <- do.call(xscat, as.list(seq_exones)) # concatenar 

        coord_exones <- tryCatch({
            exonsBy(txdb, by = "gene")[[entrez]] 
        }, error = function(e) NULL)

        if (is.null(coord_exones)) {
            result <- list(
                code = 404,
                message = paste("no se encontro secuencia para el entrez: ", entrez),
                data = NULL
            )
            return(result)
        } 
        coord_exones <- exonsBy(txdb, by = "gene")[[entrez]] 
        seq_exones <- getSeq(human_genome, coord_exones)
        DNA_str <- do.call(xscat, as.list(seq_exones)) # concatenar         
    }
    
    counter_base <- alphabetFrequency(DNA_str, baseOnly = TRUE)
    pattern_CpG <- DNAString("CG")
    counter_CpG <- countPattern(pattern_CpG, DNA_str)
    match_CpG <- matchPattern(pattern_CpG, DNA_str)
    cpg_info <- as.list(match_CpG@ranges@start)

    result <- list(
        code = 200,
        message = "Ok.",
        data = list(
            complete = as.logical(complete),
            sequence = as.character(DNA_str),
            length = nchar(DNA_str),
            nucleotides = as.list(counter_base),
            cpg_islands = list(
                count = counter_CpG,
                start = cpg_info
            )
        )
    )
}
