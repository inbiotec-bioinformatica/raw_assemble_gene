# Imagen base
FROM python:3.10-slim

LABEL maintainer="Inbiotec.bioinformática <soporte@inbiotec.bioinformática>"

# Instalar dependencias del sistema
RUN apt-get update && apt-get install -y \
    bwa \
    samtools \
    minimap2 \
    spades \
    && rm -rf /var/lib/apt/lists/*

# Crear directorio de trabajo
WORKDIR /app

# Copiar archivos
COPY . /app

# Instalar dependencias de Python
RUN pip install --no-cache-dir -r requirements.txt

# Comando por defecto
ENTRYPOINT ["python", "raw_assemble_gene.py"]
