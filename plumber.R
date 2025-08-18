library(plumber)

#* @apiTitle Vitamina D
#* @apiDescription API para cargar archivos FASTQ, alinearlos con el genoma humano (hg38), extraer fragmentos de genes relacionados con la vitamina D y analizar variantes genomicas.

api <- Plumber$new()
api$mount("/echo", Plumber$new("endpoints/echo.R"))
api$mount("/seq_by_symbol", Plumber$new("endpoints/seq_by_symbol.R"))
api$mount("/percent", Plumber$new("endpoints/percent.R"))
api$mount("/detail", Plumber$new("endpoints/detail.R"))
api$mount("/seq_by_range", Plumber$new("endpoints/seq_by_range.R"))
api$mount("/pairwise", Plumber$new("endpoints/pairwise.R"))
api$mount("/translate", Plumber$new("endpoints/translate.R"))

api$run(host = "0.0.0.0", port = 8000)

#RUN: api <- Plumber$new("bioc_r/plumber.R")