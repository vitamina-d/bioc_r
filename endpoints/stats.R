library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

library(AnnotationDbi)
library(org.Hs.eg.db)  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#* stats devuelve la secuencia completa a partir de su entrez
#* @param entrez EntrezID
#* @get /
#* @tag sequence
#* @serializer unboxedJSON 
function(entrez, res) {
    tryCatch({
        human_genome <- BSgenome.Hsapiens.UCSC.hg38
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

        if (is.null(entrez) || entrez == "") {
            res$status <- 400
            stop("Ingrese entrez.", call. = FALSE)
        }
        #coordenadas: objeto GRanges
        #secuencia: objeto DNAStringSet de biostrings
        coord_gene <- genes(txdb)[entrez]
        if (length(coord_gene) == 0) {
            res$status <- 404
            stop(paste("No se encontro el entrez: ", entrez), call. = FALSE)
        }
        sequence <- getSeq(human_genome, coord_gene) #devuelve la codificante
        DNA_str <- unlist(sequence) ##STATS

        counter_base <- alphabetFrequency(DNA_str, baseOnly = TRUE)
        pattern_CpG <- DNAString("CG")
        counter_CpG <- countPattern(pattern_CpG, DNA_str)
        match_CpG <- matchPattern(pattern_CpG, DNA_str)
        cpg_info <- as.list(match_CpG@ranges@start)

        list(
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
    }, error = function(e) {
        res$status <- 500
        stop(paste("Error de servicio R: ", e$message), call. = FALSE)
    })
}
