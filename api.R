library(plumber)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

library(AnnotationDbi)
library(org.Hs.eg.db)  
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#* @apiTitle Vitamina D
#* @apiDescription API para cargar archivos FASTQ, alinearlos con el genoma humano (hg38), extraer fragmentos de genes relacionados con la vitamina D y analizar variantes genomicas.

#* Echo back the input
#* @param msg The message to echo
#* @get /echo
#* @tag endpoints
function(msg = "") {
  list(
    msg = paste0("The message is: '", msg, "'")
  )
}

#* seqByRange devuelve la secuencia dado el cromosoma y el rango
#* @param chrom Cromosoma (ej: "dhcr7")
#* @param start Inicio
#* @param end Fin
#* @get /seqByRange
#* @tag endpoints
#* @serializer unboxedJSON 
function(chrom = "chr11", start = 71428193, end = 71452868) {
  start_time <- Sys.time()

  start <- as.integer(start)
  end <- as.integer(end)

  seq <- BSgenome.Hsapiens.UCSC.hg38[[chrom]][start:end]
  sequence <- as.character(seq)
  sequence_length <- nchar(sequence)
  
  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  result <- list(
    status = "success", 
    time_secs = time,
    data = list(
      sequence_length = sequence_length,
      sequence = sequence
    )
  )
}

#* Pairwise devuelve el alineamiento global o local de dos secuecias
#* @param pattern Lectura
#* @param subject Genoma de referencia
#* @param global:boolean Alineamiento global (TRUE) o local (FALSE)
#* @get /paiswise
#* @tag endpoints
#* @serializer unboxedJSON 
function(pattern = "", subject = "", global = TRUE) {
  start_time <- Sys.time()

  seqA <- DNAString(pattern)
  seqB <- DNAString(subject)

  type <- ifelse(global, "global", "local")

  # Alineamiento global (tipo Needleman-Wunsch)
  align <- pairwiseAlignment(	seqA, seqB, substitutionMatrix = NULL, gapOpening = -2, gapExtension = -1, type = type)

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  result <- list(
    status = "success", 
    time_secs = time,
    data = list(
      score = score(align),
      type = type,
      pattern = pattern,
      subject = subject,
      pattern_align = as.character(pattern(align)),
      subject_align = as.character(subject(align))
    )
  )
}



#* SeqBySymbol devuelve la secuencia completa o de exones, dado el symbol de un gen
#* @param gene_symbol Nombre del gen
#* @param complete:boolean Secuencia completa (TRUE) o solo exones (FALSE)
#* @get /seqBySymbol
#* @tag endpoints
#* @serializer unboxedJSON 
function(gene_symbol="DHCR7", complete = TRUE) {

  start_time <- Sys.time()

  #genoma y coord
  human_genome <- BSgenome.Hsapiens.UCSC.hg38
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

  ### identificar el gen 
  #mapIds(org.Hs.eg.db, keys = gene_entrez, column = "SYMBOL", keytype = "ENTREZID")
  #mapIds(org.Hs.eg.db, keys = gene_symbol, column = "CHR", keytype = "SYMBOL")
  gene_entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbol, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID

  ####VALIDAR

  #coordenadas: objeto GRanges
  #secuencia: objeto DNAStringSet de biostrings
  if(complete){
    
    coord_gene <- genes(txdb)[gene_entrez]
    sequence <- getSeq(human_genome, coord_gene)

  } else {
      coord_exones <- exonsBy(txdb, by = "gene")[[gene_entrez]] 
      seq_exones <- getSeq(human_genome, coord_exones)

      #CDS: secuencia de la region codificante
      sequence <- do.call(xscat, as.list(seq_exones))
  }

  type <- ifelse(complete, "complete", "exons")

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  result <- list(
    status = "success", 
    time_secs = time,
    data = list(
      type = type,
      sequence_length = nchar(sequence),
      sequence = as.character(sequence)
    )
  )
}

#* Detail muestra info de un gen, dado su symbol
#* @param symbol Nombre del gen
#* @get /detail
#* @tag endpoints
#* @serializer unboxedJSON 
function(symbol = "DHCR7") {

  start_time <- Sys.time()

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

  entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = symbol, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID
  details <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("ENSEMBL", "ENSEMBLPROT", "UNIPROT", "ENTREZID", "GENETYPE", "MAP", "SYMBOL"), keytype = "ENTREZID")

  # Obtener rangos, sin filtrar genes que están en ambas cadenas .. granges / listgranges
  #range_list <- genes(txdb, single.strand.genes.only = FALSE)[entrez]
  range <- genes(txdb)[entrez]
  range_df <- as.data.frame(range)

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))


  result <- list(
    status = "success", 
    time_secs = time,
    data = list(
      entrezID = entrez,
      symbol = symbol,
      type = unique(details$GENETYPE),
      location_chr = unique(details$MAP),
      chr = as.character(range_df$seqnames),
      start = range_df$start,
      end = range_df$end,
      length = range_df$width,
      strand = as.character(range_df$strand),
      ensembl_id_gene = unique(details$ENSEMBL),
      ensembl_id_protein = unique(details$ENSEMBLPROT),
      uniprot_id = unique(details$UNIPROT)
    )
  )
}

#* RangeBySymbol X muestra el rango de un gen, dado su symbol
#* @param symbol Nombre del gen
#* @get /rangeBySymbol
#* @tag endpoints
#* @serializer unboxedJSON 
function(symbol = "DHCR7") {

  start_time <- Sys.time()

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

  entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = symbol, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID
  details <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez, columns = c("MAP"), keytype = "ENTREZID")

  # Obtener rangos, sin filtrar genes que están en ambas cadenas .. granges / listgranges
  #range_list <- genes(txdb, single.strand.genes.only = FALSE)[entrez]
  range <- genes(txdb)[entrez]
  range_df <- as.data.frame(range)

  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))


  result <- list(
    status = "success", 
    time_secs = time,
    data = list(
      entrezID = entrez,
      location_chr = unique(details$MAP),
      chr = as.character(range_df$seqnames),
      start = range_df$start,
      end = range_df$end,
      length = range_df$width,
      strand = as.character(range_df$strand),
    )
  )
}


#* Percent muestra el porcentaje de bases, A/T y C/G de una secuencia
#* @param seq Secuencia
#* @get /percent
#* @tag endpoints
#* @serializer unboxedJSON 
function(seq = "ACGT") {

  start_time <- Sys.time()

  # biostrings

  DNA_str <- DNAString(seq)
  
  counter_base <- alphabetFrequency(DNA_str, baseOnly = TRUE)

  AT <- counter_base["A"] + counter_base["T"]
  CG <- counter_base["C"] + counter_base["G"] #termicamente mas estables

  ATCG <- AT + CG
  
  AT_percent <- (AT / ATCG) * 100
  CG_percent <- (CG / ATCG) * 100


  ####
  pattern_CpG <- "CG"
  counter_CpG <- countPattern(pattern_CpG, DNA_str)
  match_CpG <- matchPattern(pattern_CpG, DNA_str)

  CpG_ranges <- as.data.frame(match_CpG@ranges)

  ######## Convertir cada fila a lista: [{start:.., end:.., width:..}, ...]
  CpG_ranges_list <- apply(CpG_ranges, 1, function(row) {
    list(
      start = as.integer(row["start"]),
      end = as.integer(row["end"]),
      width = as.integer(row["width"])
    )
  })
   
  end_time <- Sys.time()
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  result <- list(
    status = "success", 
    time_secs = time,
    data = list(
      nucleotides = list(
        labels = c("A", "T", "C", "G"),
        counts = c(
          as.integer(counter_base["A"]),
          as.integer(counter_base["T"]),
          as.integer(counter_base["C"]),
          as.integer(counter_base["G"])
        )
      ),
      #A = counter_base["A"],
      #T = counter_base["T"],
      #C = counter_base["C"],
      #G = counter_base["G"],
      #ATCG = ATCG,
      #AT_percent = AT_percent,
      #CG_percent = CG_percent,
      counter_CpG = counter_CpG,
      CpG_ranges = CpG_ranges
    )
  )
}



#ver secuencias completas
#cat(as.character(gene_seq), "\n")
#cat(as.character(concat_exon), "\n")

 
  #obtener symbol
  #mapIds(org.Hs.eg.db, keys = gene_entrez, column = "SYMBOL", keytype = "ENTREZID")
  #obtener chr



  #entrezDHCR7 <- AnnotationDbi::select(org.Hs.eg.db, keys = "DHCR7", columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID
  #entrezTP53 <- AnnotationDbi::select(org.Hs.eg.db, keys = "TP53", columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID

  #coord_DHCR7 <- exonsBy(txdb, by = "gene")[["1717"]] 
  #coord_TP53 <- exonsBy(txdb, by = "gene")[["7157"]] 
  
  #seq_DHCR7  <- getSeq(human_genome, coord_DHCR7)
  #seq_TP53  <- getSeq(human_genome, coord_TP53)

  #CDS: secuencia de la region codificante
  #sequence_DHCR7 <- do.call(xscat, as.list(seq_DHCR7))
  #sequence_TP53 <- do.call(xscat, as.list(seq_TP53))

  #mapIds(org.Hs.eg.db, keys = gene_symbol, column = "CHR", keytype = "SYMBOL")