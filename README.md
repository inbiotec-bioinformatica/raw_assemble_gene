# raw_assemble_gene

`raw_assemble_gene` es una herramienta bioinformática desarrollada por **Inbiotec.bioinformática** para realizar el ensamblado y extensión de genes a partir de lecturas pareadas (paired-end) usando **SPAdes**, **BWA**, **Minimap2** y **Samtools**.  
Permite identificar contigs candidatos, extender regiones UTR y calcular estadísticas básicas y cobertura de ensamblado.

---

## 🚀 Instalación

### 🔹 Requisitos del sistema

- Python ≥ 3.10  
- Linux / macOS / Windows (soporte parcial en Windows)  
- Dependencias externas instaladas en `$PATH`:
  - minimap2  
  - samtools  
  - bwa  
  - spades.py  

---

### 🔹 Instalación del entorno Python (método clásico)

```bash
git clone https://github.com/Inbiotec-bioinformatica/raw_assemble_gene.git
cd raw_assemble_gene

# Crear entorno virtual (opcional)
python -m venv venv
source venv/bin/activate    # En Linux / macOS
venv\Scripts\activate       # En Windows

# Instalar dependencias
pip install -r requirements.txt

---

Si utilizás Conda, podés crear el entorno de forma reproducible con el archivo environment.yml:

conda env create -f environment.yml
conda activate raw_assemble_gene_env

---

También podés ejecutar la herramienta sin instalar nada localmente, usando el contenedor Docker incluido:

# Construir la imagen
docker build -t raw_assemble_gene .

# Ejecutar un ejemplo
docker run --rm -v $(pwd)/benchmark_example:/data raw_assemble_gene \
    bash /data/scripts/run_benchmark.sh

Esto ejecutará el benchmark completo dentro del contenedor y dejará los resultados en
benchmark_example/results/.

---

🧬 Ejemplo de uso

python raw_assemble_gene.py \
  --query genes.fasta \
  --reads1 reads_R1.fastq \
  --reads2 reads_R2.fastq \
  --extension 100 \
  --output_prefix resultado

| Parámetro         | Descripción                                        |
| ----------------- | -------------------------------------------------- |
| `--query`         | Archivo FASTA con genes de referencia              |
| `--reads1`        | Lecturas forward (FASTQ)                           |
| `--reads2`        | Lecturas reverse (FASTQ)                           |
| `--minimap2_args` | Argumentos personalizados para Minimap2 (opcional) |
| `--extension`     | Longitud de extensión de UTRs en nt (default: 0)   |
| `--output_prefix` | Prefijo para archivos de salida                    |

📁 Archivos de salida

| Archivo                             | Descripción                      |
| ----------------------------------- | -------------------------------- |
| `*.putative_gene_with_utrs.fasta`   | Secuencia ensamblada y extendida |
| `spades_output/contigs.fasta`       | Contigs generados por SPAdes     |
| `spades_output/coverage_plot.png`   | Gráfico de cobertura             |
| `spades_output/assembly_stats.tsv`  | Estadísticas del ensamblado      |
| `<output_prefix>.status_report.txt` | Informe resumen de cada gen      |

---

🧪 Test rápido

Podés verificar que todo funcione correctamente ejecutando el test automático incluido:

bash tests/test_raw_assemble_gene.sh

---

Si la instalación es correcta, deberías ver:
✅ Archivo generado correctamente: ../results/test_result.putative_gene_with_utrs.fasta
✅ Test completado con éxito

---

📊 Benchmark incluido

El repositorio incluye un ejemplo de benchmark en benchmark_example/ que compara
raw_assemble_gene frente a SPAdes, midiendo tiempo, memoria y calidad del ensamblado.

Para ejecutarlo:
cd benchmark_example/scripts
bash run_benchmark.sh

Los resultados se guardarán en benchmark_example/results/.

---

🧩 Estructura del repositorio
.
├── benchmark_example/
│   ├── data/
│   ├── results/
│   └── scripts/
├── tests/
│   └── test_raw_assemble_gene.sh
├── raw_assemble_gene.py
├── requirements.txt
├── environment.yml
├── Dockerfile
├── LICENSE.txt
└── README.md

---

🧠 Créditos

Desarrollado por el equipo de Inbiotec.bioinformática
bioinfo.inbiotec@gmail.com

Para consultas o contribuciones, crear un issue o pull request en este repositorio.
