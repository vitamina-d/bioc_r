library(AnnotationDbi)
library(org.Hs.eg.db)  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#* detail muestra info de un gen, dado su ENTREZ 
#* @param entrez id del gen
#* @get /
#* @tag detail
#* @serializer unboxedJSON 
function(entrez = "1717") {

    start_time <- Sys.time()

    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

    details <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("ENSEMBL", "ENSEMBLPROT", "UNIPROT", "ENTREZID", "GENETYPE", "MAP", "SYMBOL", "ALIAS"), keytype = "ENTREZID")

    #range <- genes(txdb)[entrez] #devuelve GRanges si son single strand (rango simple)
    #range_df <- as.data.frame(range)

    grangeslist <- genes(txdb, single.strand.genes.only = FALSE)[entrez] ##GRangesList: el gen tiene varios rangos.
    granges <- grangeslist[[1]]

    locations <- list()
    for (i in seq_along(granges)) {
        locations[[i]] <- list(
            strand = as.character(strand(granges[i])),
            seqnames = as.character(seqnames(granges[i])),
            start = start(granges[i]),
            end = end(granges[i]),
            length = width(granges[i])
        )
    }

    ensembl_id_gene <- unique(details$ENSEMBL)
    if (length(ensembl_id_gene) == 1)  {
        ensembl_id_gene <- list(ensembl_id_gene)
    } 

    end_time <- Sys.time()
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    response <- list(
        code = 200,
        message = "Ok",
        datetime = start_time,
        time_secs = time,
        data = list(
            entrez = entrez,
            alias = unique(details$ALIAS),
            symbol = unique(details$SYMBOL),
            genetype = unique(details$GENETYPE),
            citogenetic = unique(details$MAP), #principal
            location = locations,
            ensembl_id_gene = ensembl_id_gene,
            ensembl_id_protein = unique(details$ENSEMBLPROT),
            uniprot_ids = unique(details$UNIPROT)
          )
    )
}