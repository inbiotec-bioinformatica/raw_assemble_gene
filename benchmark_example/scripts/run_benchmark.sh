#!/bin/bash

# ============================================
# Benchmark: raw_assemble_gene.py vs SPAdes
# ============================================

set -e

# -------------------------
# Archivos de entrada
# -------------------------
READ1="actin_simulated1.fq"
READ2="actin_simulated2.fq"
GENE="actin.fasta"

# -------------------------
# Carpeta de resultados
# -------------------------
RESULTS_DIR="../results"

# Crear carpeta de resultados si no existe
mkdir -p "$RESULTS_DIR"

# -------------------------
# Ejecutar raw_assemble_gene.py
# -------------------------
echo "=== Ejecutando raw_assemble_gene.py ==="
/usr/bin/time -v python3 raw_assemble_gene.py \
    -i "$READ1" "$READ2" \
    -g "$GENE" \
    -o "$RESULTS_DIR/assemble_gene.putative_gene_with_utrs.fasta" \
    > "$RESULTS_DIR/assemble_gene.log" 2>&1

# -------------------------
# Ejecutar SPAdes
# -------------------------
echo "=== Ejecutando SPAdes ==="
/usr/bin/time -v spades.py \
    -1 "$READ1" -2 "$READ2" \
    -o "$RESULTS_DIR/spades_out" \
    > "$RESULTS_DIR/spades.log" 2>&1

# -------------------------
# Ejecutar BLASTN
# -------------------------
echo "=== Ejecutando BLASTN ==="

# Ensamblado raw_assemble_gene.py
blastn -query "$RESULTS_DIR/assemble_gene.putative_gene_with_utrs.fasta" \
       -subject "$GENE" \
       -outfmt 5 \
       -out "$RESULTS_DIR/assemble_gene_blast.xml"

# Ensamblado SPAdes
blastn -query "$RESULTS_DIR/spades_out/contigs.fasta" \
       -subject "$GENE" \
       -outfmt 5 \
       -out "$RESULTS_DIR/spades_blast.xml"

echo "=== Benchmark finalizado ==="
echo "Todos los resultados están en $RESULTS_DIR"
