library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

library(AnnotationDbi)
library(org.Hs.eg.db)  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#* sequence_and_stats devuelve la secuencia completa o de exones a partir de su entrez y percent
#* @param entrez EntrezID
#* @get /
#* @tag sequence
#* @serializer unboxedJSON 
function(entrez, res) {

    human_genome <- BSgenome.Hsapiens.UCSC.hg38
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    coord_gene <- tryCatch({
        genes(txdb)[entrez]
    }, error = function(e) {
        res$status <- 404
        stop(paste("No se encontro el entrez: ", entrez), call. = FALSE)
    })

    if (length(coord_gene) == 0) {
        res$status <- 404
        stop(paste("No se encontro el entrez: ", entrez), call. = FALSE)
    }
    seq_gene <- getSeq(human_genome, coord_gene)
    DNA_str <- unlist(seq_gene) 
    
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
