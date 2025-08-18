library(plumber)

#* @apiTitle Vitamina D
#* @apiDescription API para cargar archivos FASTQ, alinearlos con el genoma humano (hg38), extraer fragmentos de genes relacionados con la vitamina D y analizar variantes genomicas.

api <- Plumber$new()
api$mount("/echo", Plumber$new("endpoints/echo.R"))
api$run(host = "0.0.0.0", port = 8000)

#RUN: api <- Plumber$new("bioc_r/plumber.R")