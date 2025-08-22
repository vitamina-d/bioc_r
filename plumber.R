library(plumber)

#* @apiTitle Vitamina D
#* @apiDescription API para cargar archivos FASTQ, alinearlos con el genoma humano (hg38), extraer fragmentos de genes relacionados con la vitamina D y analizar variantes genomicas.

api <- Plumber$new()
api$mount("/isentrez", Plumber$new("endpoints/isentrez.R"))
api$mount("/detail", Plumber$new("endpoints/detail.R"))
api$mount("/seq_by_symbol", Plumber$new("endpoints/seq_by_symbol.R"))
api$mount("/percent", Plumber$new("endpoints/percent.R"))
api$mount("/align", Plumber$new("endpoints/align.R"))
api$mount("/seq_by_range", Plumber$new("endpoints/seq_by_range.R"))
api$mount("/echo", Plumber$new("endpoints/echo.R"))

api$run(host = "0.0.0.0", port = 8000)

#RUN: api <- Plumber$new("bioc_r/plumber.R")
#"ClustalW", "ClustalOmega" o "Muscle". Todos son algoritmos clásicos de MSA