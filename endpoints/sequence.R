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
function(entrez, complete = TRUE, res) {
    tryCatch({
        human_genome <- BSgenome.Hsapiens.UCSC.hg38
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

        if (is.null(entrez) || entrez == "") {
            res$status <- 400
            stop("Ingrese entrez.", call. = FALSE)
        }
        if(complete){      
            
            coord_gene <- tryCatch({
                genes(txdb)[entrez]
            }, error = function(e) NULL)
            
            if (length(coord_gene) == 0) {
                res$status <- 404
                stop(paste("No se encontro el entrez: ", entrez), call. = FALSE)
            }
            seq_gene <- getSeq(human_genome, coord_gene)

            result <- list(
                list(
                    index = 1,
                    start = start(coord_gene)[1],
                    end   = end(coord_gene)[1],
                    sequence_length = nchar(seq_gene),
                    sequence = as.character(seq_gene)
                )
            )
        } else {
        
            coord_exones <- tryCatch({
                exonsBy(txdb, by = "gene")[[entrez]] 
            }, error = function(e) NULL)
            
            if (is.null(coord_exones)) {
                res$status <- 404
                stop(paste("No se encontro secuencia para el entrez: ", entrez), call. = FALSE)
            } 
            coord_exones <- exonsBy(txdb, by = "gene")[[entrez]] 
            seq_exones <- getSeq(human_genome, coord_exones)

            result <- lapply(seq_along(seq_exones), function(i) {
                list(
                    index = i,
                    start = start(coord_exones)[i],
                    end = end(coord_exones)[i],
                    sequence_length = nchar(seq_exones[i]),
                    sequence = as.character(seq_exones[i])
                )
            })
        }

        list(
            code = 200,
            message = "Ok.",
            data = result
        )

    }, error = function(e) {
        res$status <- 500
        stop(paste("Error de servicio R: ", e$message), call. = FALSE)
    })
}