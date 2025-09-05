
FROM veroyols/myapp_bioc:librerias

WORKDIR /rservice
COPY . .
EXPOSE 8000 8787

CMD ["R", "-e", "library(plumber); api <- Plumber$new('plumber.R'); api$run(host='0.0.0.0', port=8000)"]

#RUN: api <- Plumber$new("bioc_r/plumber.R")
