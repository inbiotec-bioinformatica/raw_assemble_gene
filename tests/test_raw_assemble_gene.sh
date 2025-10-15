#!/bin/bash

# ============================================
# Test automático de raw_assemble_gene.py
# ============================================

set -e  # detener ejecución ante cualquier error

# Ir al directorio del benchmark de ejemplo
cd "$(dirname "$0")/../benchmark_example/scripts"

# Archivos de entrada
READ1="../data/actin_simulated1.fq"
READ2="../data/actin_simulated2.fq"
QUERY="../data/actin.fasta"
OUTPUT_PREFIX="../results/test_result"
PLATFORM="illumina"
EXTENSION=200

echo "=== Ejecutando prueba de raw_assemble_gene.py ==="

python3 raw_assemble_gene.py \
    --reads1 "$READ1" \
    --reads2 "$READ2" \
    --query "$QUERY" \
    --output_prefix "$OUTPUT_PREFIX" \
    --extension "$EXTENSION" \
    --platform "$PLATFORM"

# Verificar que se haya generado el archivo esperado
OUTPUT_FILE="${OUTPUT_PREFIX}.putative_gene_with_utrs.fasta"

if [ -f "$OUTPUT_FILE" ]; then
    echo "✅ Archivo generado correctamente: $OUTPUT_FILE"
else
    echo "❌ Error: no se generó el archivo de salida esperado ($OUTPUT_FILE)"
    exit 1
fi

# Verificar que el archivo no esté vacío
if [ ! -s "$OUTPUT_FILE" ]; then
    echo "❌ Error: el archivo de salida está vacío"
    exit 1
fi

echo "✅ Test completado con éxito"
