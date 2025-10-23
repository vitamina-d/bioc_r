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
api$mount("/autocomplete", Plumber$new("endpoints/autocomplete.R"))


api$mount("/align", Plumber$new("endpoints/align.R"))
api$mount("/complement", Plumber$new("endpoints/complement.R"))
api$mount("/detail", Plumber$new("endpoints/detail.R"))
api$mount("/detailfull", Plumber$new("endpoints/detailfull.R"))
api$mount("/isentrez", Plumber$new("endpoints/isentrez.R"))
api$mount("/sequence", Plumber$new("endpoints/sequence.R"))
api$mount("/percent", Plumber$new("endpoints/percent.R"))
api$mount("/seq_by_range", Plumber$new("endpoints/seq_by_range.R"))
api$mount("/stats", Plumber$new("endpoints/sequence_and_stats.R"))

api$mount("/echo", Plumber$new("endpoints/echo.R"))

api$mount("/genename", Plumber$new("endpoints/genename.R"))
api$mount("/table", Plumber$new("endpoints/table.R"))

api$mount("/entrez", Plumber$new("endpoints/entrez.R"))
api$mount("/entrezByAlias", Plumber$new("endpoints/entrezByAlias.R"))
api$mount("/entrezBySymbol", Plumber$new("endpoints/entrezBySymbol.R"))

api$run(host = "0.0.0.0", port = 8000)

#RUN: api <- Plumber$new("bioc_r/plumber.R")
#"ClustalW", "ClustalOmega" o "Muscle". Todos son algoritmos clásicos de MSA
