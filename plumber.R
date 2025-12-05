library(plumber)

#* @apiTitle Vitamina D
#* @apiDescription API para cargar archivos FASTQ, alinearlos con el genoma humano (hg38), extraer fragmentos de genes relacionados con la vitamina D y analizar variantes genomicas.

api <- Plumber$new()

api$mount("/autocomplete", Plumber$new("endpoints/autocomplete.R"))
api$mount("/detail", Plumber$new("endpoints/detail.R")) 
api$mount("/detailfull", Plumber$new("endpoints/detailfull.R")) 
api$mount("/entrez", Plumber$new("endpoints/entrez.R"))
api$mount("/isentrez", Plumber$new("endpoints/isentrez.R"))
api$mount("/stats", Plumber$new("endpoints/stats.R"))
api$mount("/sequence_range", Plumber$new("endpoints/sequence_range.R"))

api$mount("/align", Plumber$new("endpoints/align.R")) 
api$mount("/percent", Plumber$new("endpoints/percent.R"))
api$mount("/sequence", Plumber$new("endpoints/sequence.R"))
api$mount("/sequence_and_stats", Plumber$new("endpoints/sequence_and_stats.R"))

api$run(host = "0.0.0.0", port = 8000)