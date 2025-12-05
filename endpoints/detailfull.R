library(AnnotationDbi)
library(org.Hs.eg.db)  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#* detail muestra info de un gen, dado su ENTREZ 
#* @param entrez id del gen
#* @get /
#* @tag detail
#* @serializer unboxedJSON 
function(entrez, res) {

    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    detail <- tryCatch({
        AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("ENTREZID", "GENETYPE", "GENENAME", "MAP", "SYMBOL"), keytype = "ENTREZID")
    }, error = function(e) {
        res$status <- 500
        stop(paste("Error de servicio R: ", e$message), call. = FALSE)
    })

    if (is.null(detail) || length(detail) == 0 || nrow(detail) == 0) {
        res$status <- 404
        stop(paste("No se encontro el entrez: ", entrez), call. = FALSE)
    }

    aliases <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("ALIAS"), keytype = "ENTREZID")
    aliases <- unique(aliases$ALIAS)
    aliases <- aliases[!is.na(aliases)]
    if (length(aliases) == 0) {
        aliases <- list()
    } else if (length(aliases) == 1) {
        aliases <- list(aliases)
    }

    prot_ensembl <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("ENSEMBLPROT"), keytype = "ENTREZID")
    prot_ensembl <- unique(prot_ensembl$ENSEMBLPROT)
    prot_ensembl <- prot_ensembl[!is.na(prot_ensembl)]
    if (length(prot_ensembl) == 0) {
        prot_ensembl <- list()
    } else if (length(prot_ensembl) == 1) {
        prot_ensembl <- list(prot_ensembl)
    }

    gene_ensembl <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("ENSEMBL"), keytype = "ENTREZID")
    gene_ensembl <- unique(gene_ensembl$ENSEMBL)
    gene_ensembl <- gene_ensembl[!is.na(gene_ensembl)]
    if (length(gene_ensembl) == 0) {
        gene_ensembl <- list()
    } else if (length(gene_ensembl) == 1) {
        gene_ensembl <- list(gene_ensembl)
    }
    
    uniprot <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("UNIPROT"), keytype = "ENTREZID")
    uniprot <- unique(uniprot$UNIPROT)
    uniprot <- uniprot[!is.na(uniprot)]
    if (length(uniprot) == 0) {
        uniprot <- list()
    } else if (length(uniprot) == 1) {
        uniprot <- list(uniprot)
    }
    
    grangeslist <- tryCatch({
        genes(txdb, single.strand.genes.only = FALSE)[entrez]
    }, error = function(e) {NULL})

    locations <- list()
    if (!is.null(grangeslist) && length(grangeslist) > 0) {
        granges <- grangeslist[[1]]
        for (i in seq_along(granges)) {
            locations[[i]] <- list(
                strand = as.character(strand(granges[i])),
                seqnames = as.character(seqnames(granges[i])),
                start = start(granges[i]),
                end = end(granges[i]),
                length = width(granges[i])
            )
        }
    } 

    if (is.null(detail) || length(detail) == 0) {
        res$status <- 404
        stop(paste("No se obtuvo detalle: "), call. = FALSE)
    }

    response <- list(
        code = 200,
        message = "Ok",
        data = list(
            entrez = entrez,
            symbol = unique(detail$SYMBOL),
            genetype = unique(detail$GENETYPE),
            genename = unique(detail$GENENAME),
            citogenetic = unique(detail$MAP), #principal
            location = locations, #list
            alias = aliases, #list
            ensembl_id_gene = gene_ensembl, #list
            ensembl_id_protein = prot_ensembl, #list
            uniprot_ids = uniprot #list
          )
    )
}