# 🧬 bioc_r

Servicio de Bioconductor en R del proyecto **Vitamina-D** para análisis genómico y procesamiento de secuencias.

## 📋 Descripción

**bioc_r** es el motor de análisis genómico construido sobre Bioconductor y expuesto como API REST mediante Plumber. Proporciona acceso a genomas de referencia, análisis de secuencias, alineamientos y búsqueda de genes relacionados con la vitamina D.

Este servicio forma parte del **proyecto Vitamina** y se integra con:

- **[bioc_front](https://github.com/vitamina-d/bioc_front)** - Frontend web para visualización
- **[bioc_back](https://github.com/vitamina-d/bioc_back)** - API intermedia en ASP.NET Core

## ✨ Funcionalidades Principales

### 🔍 Búsqueda y Consulta de Secuencias

- Extracción de secuencias genómicas por cromosoma y rango
- Autocompletado de genes por símbolo
- Validación de IDs Entrez
- Búsqueda de información detallada de genes

### 🧪 Análisis de Secuencias

- Alineamiento múltiple de secuencias
- Cálculo de composición nucleotídica (porcentajes de A, T, G, C)
- Análisis de fragmentos específicos del genoma

### 🧬 Genomas y Anotaciones

- Acceso a genomas de referencia vía BSGenome
- Genoma humano hg38 (GRCh38)
- Anotaciones genómicas de genes relacionados con vitamina D
- Integración con bases de datos NCBI

### 📊 RStudio Server

- Interfaz web de RStudio disponible en puerto 8787
- Acceso interactivo al entorno R
- Desarrollo y testing de endpoints

## 🛠️ Stack Tecnológico

### Core

- **R** - Lenguaje de programación estadístico
- **Bioconductor** - Suite de paquetes para análisis bioinformático
- **Plumber** - Framework para APIs REST en R

### Paquetes Bioconductor Principales

- **BSgenome** - Genomas de referencia
- **BSgenome.Hsapiens.UCSC.hg38** - Genoma humano GRCh38
- **GenomicRanges** - Manipulación de rangos genómicos
- **Biostrings** - Manejo de secuencias biológicas
- **AnnotationDbi** - Acceso a bases de datos de anotación

### Infraestructura

- **Docker** - Contenedorización
- **RStudio Server** - IDE web (puerto 8787)
- **Imagen base**: `veroyols/myapp_bioc:librerias`

## 📁 Estructura del Proyecto

```
bioc_r/
├── endpoints/
│   ├── sequence_range.R    # Extracción de secuencias por rango
│   ├── align.R             # Alineamiento de secuencias
│   ├── detail.R            # Información básica de genes
│   ├── detailfull.R        # Información completa de genes
│   ├── percent.R           # Composición nucleotídica
│   ├── autocomplete.R      # Búsqueda autocompletada
│   ├── entrez.R            # Búsqueda por Entrez ID
│   └── isentrez.R          # Validación de Entrez ID
├── plumber.R               # API principal
├── Dockerfile              # Configuración Docker
└── README.md               # Este archivo
```

## 🚀 Inicio Rápido

### Con Docker

```bash
git clone https://github.com/vitamina-d/bioc_r.git
cd bioc_r

docker build -t bioc_r .

docker run -p 8000:8000 -p 8787:8787 \
  -e PASSWORD=bioc \
  -e USER=rstudio \
  bioc_r
```

### Acceso a los Servicios

- **API REST**: http://localhost:8000
- **RStudio Server**: http://localhost:8787
  - Usuario: `rstudio`
  - Password: `bioc`

## 📡 API Reference

### 1. Extracción de Secuencias

**GET** `/sequence_range/`

Extrae una secuencia genómica del genoma humano hg38.

#### Parámetros

| Parámetro | Tipo    | Descripción                    |
| --------- | ------- | ------------------------------ |
| chrom     | string  | Cromosoma (ej: "chr1", "chrX") |
| start     | integer | Posición inicial               |
| end       | integer | Posición final                 |

#### Ejemplo

```bash
curl "http://localhost:8000/sequence_range/?chrom=chr1&start=100000&end=100100"
```

#### Response

```json
{
  "code": 200,
  "message": "Ok",
  "data": {
    "sequence": "ATGGCTAGCTAG...",
    "length": 100
  }
}
```

---

### 2. Alineamiento de Secuencias

**POST** `/align/`

Realiza alineamiento de secuencias.

#### Request Body

```json
{
  "pattern": "ATGGCTAGCTAG",
  "subject": "ATGACTACCTAG",
  "type": "local",
  "gapOpening": 1,
  "gapExtension": 1
}
```

#### Response

```json
{
  "code": 200,
  "message": "Ok",
  "data": {
    "score ": 1,
    "pattern_align": "...",
    "subject_align": "..."
  }
}
```

---

### 3. Información de Genes

**GET** `/detail/`

Obtiene información básica de un gen por símbolo.

#### Parámetros

| Parámetro | Tipo   | Descripción                            |
| --------- | ------ | -------------------------------------- |
| symbol    | string | Símbolo del gen (ej: "VDR", "CYP27B1") |

#### Ejemplo

```bash
curl "http://localhost:8000/detail/?symbol=VDR"
```

---

**GET** `/detailfull/`

Obtiene información completa de un gen incluyendo anotaciones.

#### Parámetros

| Parámetro | Tipo   | Descripción     |
| --------- | ------ | --------------- |
| symbol    | string | Símbolo del gen |

---

### 4. Composición Nucleotídica

**GET** `/percent/`

Calcula el porcentaje de cada nucleótido en una secuencia.

#### Parámetros

| Parámetro | Tipo   | Descripción              |
| --------- | ------ | ------------------------ |
| sequence  | string | Secuencia de nucleótidos |

#### Ejemplo

```bash
curl "http://localhost:8000/percent/?sequence=ATGGCTAGCTAG"
```

#### Response

```json
{
  "code": 200,
  "message": "Ok",
  "data": {
    "sequence_length": 1,
    "nucleotides": {
      "A": 25.0,
      "T": 25.0,
      "G": 33.3,
      "C": 16.7
    },
    "cpg_islands":{
      "count": 1,
      "start": [1, 2, 3]
    },
    "sequence": "ATGGCTAGCTAG"
  }
}
```

---

### 5. Autocompletado

**GET** `/autocomplete/`

Busca los alias que coincidan con un input para el autocompletado.

#### Parámetros

| Parámetro | Tipo    | Descripción                     |
| --------- | ------- | ------------------------------- |
| input     | string  | Prefijo del alias del gen     |

#### Ejemplo

```bash
curl "http://localhost:8000/autocomplete/?input=VD"
```

#### Response

```json
{
  "code": 200,
  "message": "Ok",
  "data": [ "VDR", ... , "VDB"]
}
```

---

### 6. Búsqueda por Entrez ID

**GET** `/entrez/`

Obtiene información de un gen por su Entrez ID.

#### Parámetros

| Parámetro | Tipo    | Descripción    |
| --------- | ------- | -------------- |
| id        | integer | Entrez Gene ID |

#### Ejemplo

```bash
curl "http://localhost:8000/entrez/?id=7421"
```

---

**GET** `/isentrez/`

Valida si un ID es un Entrez ID válido.

#### Parámetros

| Parámetro | Tipo    | Descripción  |
| --------- | ------- | ------------ |
| id        | integer | ID a validar |

#### Response

```json
{
  "code": 200,
  "message": "Ok",
  "data": {
    "valid": true,
    "symbol": "VDR"
  }
}
```

## 🔧 Configuración

### Dockerfile

```dockerfile
FROM veroyols/myapp_bioc:librerias

WORKDIR /rservice
COPY . .
EXPOSE 8000 8787

CMD ["R", "-e", "library(plumber); api <- Plumber$new('plumber.R'); api$run(host='0.0.0.0', port=8000)"]
```

## 🔗 Integración 

### Flujo de Datos

```
Frontend (React)
      ↓
   Backend (ASP.NET)
      ↓
   bioc_r (R/Plumber) ←→ BSGenome (hg38)
      ↓
   Bioconductor Packages
```

### Servicios Relacionados

- **[bioc_back](https://github.com/vitamina-d/bioc_back)** - Consume los endpoints de este servicio
- **[bioc_front](https://github.com/vitamina-d/bioc_front)** - Visualiza los resultados
- **[bioc_blast](https://github.com/vitamina-d/bioc_blast)** - Análisis BLAST complementario

## 📊 Datos Genómicos

### Genoma de Referencia

- **Versión**: GRCh38 (hg38)
- **Organismo**: Homo sapiens
- **Fuente**: UCSC Genome Browser
- **Paquete**: BSgenome.Hsapiens.UCSC.hg38

## 📝 Licencia

Este proyecto tiene fines educativos y forma parte del Proyecto Integrador Profesional (PIP).
