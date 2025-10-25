library(plumber)

#* @apiTitle Vitamina D
#* @apiDescription API para cargar archivos FASTQ, alinearlos con el genoma humano (hg38), extraer fragmentos de genes relacionados con la vitamina D y analizar variantes genomicas.

### AGReGAR COMPARTIDAS
#   response <- list(
#     code = 200,
#     datetime = start_time,
#     time_secs = time,
#     data = list( ... )
#   )

api <- Plumber$new()
api$mount("/autocomplete", Plumber$new("endpoints/autocomplete.R"))#ok
api$mount("/align", Plumber$new("endpoints/align.R")) #ok
api$mount("/detail", Plumber$new("endpoints/detail.R")) #ok

api$mount("/detailfull", Plumber$new("endpoints/detailfull.R"))
api$mount("/entrez", Plumber$new("endpoints/entrez.R"))
api$mount("/isentrez", Plumber$new("endpoints/isentrez.R"))
api$mount("/percent", Plumber$new("endpoints/percent.R"))
api$mount("/seq_by_range", Plumber$new("endpoints/seq_by_range.R"))
api$mount("/sequence", Plumber$new("endpoints/sequence.R"))
api$mount("/stats", Plumber$new("endpoints/sequence_and_stats.R"))

api$run(host = "0.0.0.0", port = 8000)

#RUN: api <- Plumber$new("bioc_r/plumber.R")
#"ClustalW", "ClustalOmega" o "Muscle". Todos son algoritmos clásicos de MSA
