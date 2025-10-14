# Benchmark Example

Esta carpeta contiene un ejemplo reproducible para evaluar el desempeño de `raw_assemble_gene` frente a SPAdes estándar.
# Benchmark Example: raw_assemble_gene.py vs SPAdes

Este ejemplo permite comparar el desempeño de:

- **raw_assemble_gene.py**   
- **SPAdes** (ensamblador de referencia)

usando datos simulados de RNA-Seq y el gen de referencia `actin.fasta`.

---

## Estructura del directorio

benchmark_example/
├── data/
│ ├── actin_simulated1.fq
│ ├── actin_simulated2.fq
│ └── actin.fasta
├── scripts/
│ ├── raw_assemble_gene.py
│ └── run_benchmark.sh
├── results/
│ ├── assemble_gene.putative_gene_with_utrs.fasta
│ ├── spades_out/
│ ├── *.log
│ ├── *_blast.xml
│ └── benchmark_report.md
└── README.md


---

## Requisitos

- Linux  
- Python 3.x  
- SPAdes instalado y accesible desde la línea de comando  
- BLAST+ (`blastn`) instalado

---

## Cómo ejecutar el benchmark

Desde la carpeta `scripts/`:

```bash
bash run_benchmark.sh


---

## 3️⃣ `results/benchmark_report.md`

```markdown
# Benchmark Report: raw_assemble_gene.py vs SPAdes

Este reporte describe el benchmark usando datos simulados de RNA-Seq (`actin_simulated1.fq` y `actin_simulated2.fq`) y el gen de referencia `actin.fasta`.

---

## Objetivo

Comparar:

1. **raw_assemble_gene.py**  
2. **SPAdes**

en términos de:

- Tiempo de ejecución  
- Memoria utilizada  
- Calidad del ensamblado (medida mediante BLASTN contra el gen de referencia)

---

## Pasos realizados

1. Ejecutar `raw_assemble_gene.py` con los reads simulados.  
2. Ejecutar SPAdes con los mismos reads.  
3. Medir tiempo y memoria con `/usr/bin/time -v`.  
4. Evaluar los ensamblados con `blastn` contra `actin.fasta`.  
5. Guardar todos los resultados en la carpeta `results/`.

---

## Archivos generados

- `assemble_gene.putative_gene_with_utrs.fasta`: ensamblado de raw_assemble_gene.py  
- `spades_out/contigs.fasta`: ensamblado de SPAdes  
- `*.log`: logs de tiempo y memoria  
- `*_blast.xml`: resultados de BLASTN  
- `benchmark_report.md`: descripción de este benchmark

---

## Cómo interpretar los resultados

- Revisar los logs (`*.log`) para ver tiempo de ejecución y memoria máxima usada.  
- Analizar los archivos BLAST XML (`*_blast.xml`) para verificar la calidad del ensamblado.  
- Comparar los contigs generados para ver cuál ensamblador reproduce mejor la secuencia de referencia.

---

## Conclusión

Este benchmark está diseñado para alguien que nunca usó SPAdes o raw_assemble_gene.py, y permite:

- Evaluar velocidad y uso de memoria  
- Comparar calidad de ensamblado  
- Tener un flujo reproducible paso a paso
