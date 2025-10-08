# raw_assemble_gene

**raw_assemble_gene** es una herramienta bioinformática desarrollada por **Inbiotec.bioinformática** para realizar el **ensamblado y extensión de genes** a partir de lecturas pareadas (paired-end) usando *SPAdes*, *BWA*, *Minimap2* y *Samtools*.  
Permite identificar contigs candidatos, extender regiones UTR y calcular estadísticas básicas y cobertura de ensamblado.

---

## 🚀 Instalación

### 🔹 Requisitos del sistema
- Python ≥ 3.10  
- Linux / macOS / Windows (soporte parcial en Windows)
- Dependencias externas instaladas en `$PATH`:
  - `minimap2`
  - `samtools`
  - `bwa`
  - `spades.py`

### 🔹 Instalación del entorno Python
```bash
git clone https://github.com/Inbiotec-bioinformatica/raw_assemble_gene.git
cd raw_assemble_gene
python -m venv venv
venv\Scripts\activate  # En Windows
pip install -r requirements.txt

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


| Archivo                             | Descripción                      |
| ----------------------------------- | -------------------------------- |
| `*.putative_gene_with_utrs.fasta`   | Secuencia ensamblada y extendida |
| `spades_output/contigs.fasta`       | Contigs generados por SPAdes     |
| `spades_output/coverage_plot.png`   | Gráfico de cobertura             |
| `spades_output/assembly_stats.tsv`  | Estadísticas del ensamblado      |
| `<output_prefix>.status_report.txt` | Informe resumen de cada gen      |

