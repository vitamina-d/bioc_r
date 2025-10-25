library(AnnotationDbi)
library(org.Hs.eg.db)  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#* detail muestra info de un gen, dado su ENTREZ 
#* @param entrez id del gen
#* @get /
#* @tag detail
#* @serializer unboxedJSON 
function(entrez) {

    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    detail <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("ENTREZID", "GENETYPE", "GENENAME", "MAP", "SYMBOL"), keytype = "ENTREZID")

    #manejo listas: [], ["item"], ["item1", "item2"]
    aliases <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("ALIAS"), keytype = "ENTREZID")
    print("aliases")
    print(aliases$ALIAS)
    aliases <- unique(aliases$ALIAS)
    aliases <- aliases[!is.na(aliases)]
    if (length(aliases) == 0) {
        aliases <- list()
    } else if (length(aliases) == 1) {
        aliases <- list(aliases)
    }

    prot_ensembl <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("ENSEMBLPROT"), keytype = "ENTREZID")
    print("prot_ensembl")
    print(prot_ensembl$ENSEMBLPROT)
    prot_ensembl <- unique(prot_ensembl$ENSEMBLPROT)
    prot_ensembl <- prot_ensembl[!is.na(prot_ensembl)]
    if (length(prot_ensembl) == 0) {
        prot_ensembl <- list()
    } else if (length(prot_ensembl) == 1) {
        prot_ensembl <- list(prot_ensembl)
    }

    gene_ensembl <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("ENSEMBL"), keytype = "ENTREZID")
    print("gene_ensembl")
    print(gene_ensembl$ENSEMBL)
    gene_ensembl <- unique(gene_ensembl$ENSEMBL)
    gene_ensembl <- gene_ensembl[!is.na(gene_ensembl)]
    if (length(gene_ensembl) == 0) {
        gene_ensembl <- list()
    } else if (length(gene_ensembl) == 1) {
        gene_ensembl <- list(gene_ensembl)
    }
    
    uniprot <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("UNIPROT"), keytype = "ENTREZID")
    print("uniprot")
    print(uniprot$UNIPROT)
    uniprot <- unique(uniprot$UNIPROT)
    uniprot <- uniprot[!is.na(uniprot)]
    if (length(uniprot) == 0) {
        uniprot <- list()
    } else if (length(uniprot) == 1) {
        uniprot <- list(uniprot)
    }
    
    # locations
    grangeslist <- tryCatch({
        genes(txdb, single.strand.genes.only = FALSE)[entrez]
    }, error = function(e) NULL)
    print("grangeslist")
    print(grangeslist)

    locations <- list()
    if (!is.null(grangeslist) && length(grangeslist) > 0) {
        granges <- grangeslist[[1]]
        print("granges")
        print(granges)
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
        # exception
        stop(paste("No se obtuvo detalle: ", e$message), call. = FALSE)
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

    #range <- genes(txdb)[entrez] #devuelve GRanges si son single strand (rango simple)
    #range_df <- as.data.frame(range)

    #grangeslist <- genes(txdb, single.strand.genes.only = FALSE)[entrez] ##GRangesList: el gen tiene varios rangos.
    #granges <- grangeslist[[1]]
}