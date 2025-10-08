#!/usr/bin/env python3
"""
Asamblea y extensión de genes desde lecturas pareadas
Autor: [Tu Nombre]
Versión: 1.1 (2025-10)
"""

import os
import sys
import shutil
import argparse
import subprocess
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import matplotlib.pyplot as plt


# -------------------------------------------------
# Utilidades generales
# -------------------------------------------------

def check_dependencies(tools):
    """Verifica que las herramientas externas estén disponibles."""
    for tool in tools:
        if shutil.which(tool) is None:
            sys.exit(f"[ERROR] No se encontró '{tool}' en el PATH. Instálalo o agrega su ruta.")


def run_command(command, description):
    """Ejecuta un comando y detiene si falla."""
    print(f"[INFO] Ejecutando: {description}")
    print(f"        {' '.join(command)}")
    result = subprocess.run(command)
    if result.returncode != 0:
        sys.exit(f"[ERROR] Falló la ejecución: {description}")
    print(f"[OK] {description} completado.")


def validate_file(path, label):
    """Valida existencia y tamaño >0."""
    if not os.path.exists(path):
        sys.exit(f"[ERROR] No se encontró {label}: {path}")
    if os.path.getsize(path) == 0:
        sys.exit(f"[ERROR] El archivo {label} está vacío: {path}")


# -------------------------------------------------
# Funciones principales
# -------------------------------------------------

def get_aligned_reads(query, reads1, reads2, minimap2_args, output_prefix):
    """Alinea las lecturas al gen y genera FASTQ pareadas."""
    sam_file = f"{output_prefix}.sam"
    cmd = ["minimap2", "-ax", "sr", *minimap2_args.split(), query, reads1, reads2]
    with open(sam_file, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)

    # Generar fastq pareadas
    fastq1 = f"{output_prefix}_R1.fastq"
    fastq2 = f"{output_prefix}_R2.fastq"
    cmd = [
        "samtools", "fastq",
        "-f", "0x2",
        "-1", fastq1,
        "-2", fastq2,
        sam_file
    ]
    subprocess.run(cmd, check=True)
    validate_file(fastq1, "FASTQ R1")
    validate_file(fastq2, "FASTQ R2")
    return fastq1, fastq2


def run_spades(reads1, reads2, outdir):
    """Ejecuta SPAdes con lecturas pareadas."""
    print(f"[INFO] Ensamblando con SPAdes en {outdir}")
    os.makedirs(outdir, exist_ok=True)
    command = ["spades.py", "-1", reads1, "-2", reads2, "-o", outdir]
    run_command(command, "SPAdes")
    contigs_path = os.path.join(outdir, "contigs.fasta")
    validate_file(contigs_path, "contigs.fasta")
    return contigs_path


def find_best_contig(query_seq, contigs_path):
    """Selecciona el contig con mayor score local vs query."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    min_score = len(query_seq) * 0.75
    best_score = 0
    best_contig = None

    for contig in SeqIO.parse(contigs_path, "fasta"):
        score = aligner.score(contig.seq, query_seq.seq)
        if score > best_score:
            best_score = score
            best_contig = contig

    if best_contig is None:
        raise ValueError(f"No se encontró contig para {query_seq.id}")

    if best_score < min_score:
        print(f"[WARNING] {query_seq.id}: solo similitud parcial (score={best_score:.1f}).")

    return best_contig


def extend_contig_with_utrs(best_contig, query_seq, extension_length):
    """Extiende el contig con UTRs en minúsculas."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    alignment = aligner.align(str(best_contig.seq), str(query_seq.seq))[0]

    if alignment.score == 0 or len(alignment.aligned) == 0:
        raise ValueError("El gen no se alinea con el contig.")

    contig_aligned_regions, _ = alignment.aligned
    start = contig_aligned_regions[0][0]
    end = contig_aligned_regions[-1][1]

    extended_start = max(0, start - extension_length)
    extended_end = min(len(best_contig), end + extension_length)

    contig_seq = str(best_contig.seq)
    extended_seq = contig_seq[extended_start:extended_end]

    gene_start = start - extended_start
    gene_end = gene_start + (end - start)

    annotated_seq = (
        extended_seq[:gene_start].lower()
        + extended_seq[gene_start:gene_end].upper()
        + extended_seq[gene_end:].lower()
    )

    if gene_start == 0:
        print(f"[WARNING] {query_seq.id}: no se pudo extender el 5'-UTR.")
    if gene_end == len(annotated_seq):
        print(f"[WARNING] {query_seq.id}: no se pudo extender el 3'-UTR.")

    return SeqRecord(Seq(annotated_seq), id=query_seq.id, description="gene_with_utrs")


def analyze_coverage(contigs_fasta, reads1, reads2, outdir):
    """Calcula cobertura por contig y genera gráfico."""
    print("[INFO] Analizando cobertura del ensamblado...")
    os.makedirs(outdir, exist_ok=True)

    bam_file = os.path.join(outdir, "mapped.bam")
    sorted_bam = os.path.join(outdir, "mapped.sorted.bam")
    cov_file = os.path.join(outdir, "coverage_summary.tsv")
    plot_file = os.path.join(outdir, "coverage_plot.png")

    # Mapeo
    subprocess.run(["bwa", "index", contigs_fasta], check=True)
    p1 = subprocess.Popen(["bwa", "mem", contigs_fasta, reads1, reads2], stdout=subprocess.PIPE)
    with open(bam_file, "wb") as bam_out:
        p2 = subprocess.Popen(["samtools", "view", "-bS", "-"], stdin=p1.stdout, stdout=bam_out)
        p1.stdout.close()
        p2.communicate()

    subprocess.run(["samtools", "sort", bam_file, "-o", sorted_bam], check=True)
    subprocess.run(["samtools", "index", sorted_bam], check=True)

    # Cobertura
    depth_out = subprocess.check_output(["samtools", "depth", sorted_bam])
    depth_lines = depth_out.decode().strip().split("\n")
    data = {}

    for line in depth_lines:
        contig, pos, depth = line.split("\t")
        data.setdefault(contig, []).append(int(depth))

    summary = []
    plt.figure(figsize=(10, 6))
    for contig, depths in data.items():
        summary.append({
            "contig": contig,
            "mean_coverage": sum(depths) / len(depths),
            "max_coverage": max(depths),
            "min_coverage": min(depths),
        })
        plt.plot(range(len(depths)), depths, label=contig)

    pd.DataFrame(summary).to_csv(cov_file, sep="\t", index=False)
    plt.xlabel("Posición")
    plt.ylabel("Cobertura")
    plt.title("Cobertura por contig")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_file)
    plt.close()


def assembly_stats(contigs_fasta, outdir):
    """Calcula estadísticas básicas del ensamblado."""
    lengths = [len(rec.seq) for rec in SeqIO.parse(contigs_fasta, "fasta")]
    if not lengths:
        raise ValueError("No se encontraron contigs.")

    total_length = sum(lengths)
    max_length = max(lengths)
    lengths.sort(reverse=True)
    acc = 0
    n50 = 0
    for l in lengths:
        acc += l
        if acc >= total_length / 2:
            n50 = l
            break

    stats_file = os.path.join(outdir, "assembly_stats.tsv")
    with open(stats_file, "w") as f:
        f.write("total_contigs\ttotal_length\tmax_length\tN50\n")
        f.write(f"{len(lengths)}\t{total_length}\t{max_length}\t{n50}\n")


# -------------------------------------------------
# MAIN
# -------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Asamblea y extensión de gen desde lecturas pareadas")
    parser.add_argument("--query", required=True)
    parser.add_argument("--reads1", required=True)
    parser.add_argument("--reads2", required=True)
    parser.add_argument("--minimap2_args", default="")
    parser.add_argument("--extension", type=int, default=0)
    parser.add_argument("--output_prefix", default="genes")
    args = parser.parse_args()

    # Validaciones iniciales
    for path, label in [
        (args.query, "Archivo FASTA de genes"),
        (args.reads1, "Lecturas forward"),
        (args.reads2, "Lecturas reverse"),
    ]:
        validate_file(path, label)

    check_dependencies(["minimap2", "samtools", "bwa", "spades.py"])

    status_lines = []

    for query_seq in SeqIO.parse(args.query, "fasta"):
        gene_id = query_seq.id
        print(f"\n[INFO] Procesando gen: {gene_id}")
        try:
            reads1, reads2 = get_aligned_reads(args.query, args.reads1, args.reads2, args.minimap2_args, gene_id)
            spades_dir = os.path.join(f"{gene_id}.output", "spades_output")
            contigs_path = run_spades(reads1, reads2, spades_dir)
            best_contig = find_best_contig(query_seq, contigs_path)
            annotated = extend_contig_with_utrs(best_contig, query_seq, args.extension)
            output_fasta = f"{gene_id}.putative_gene_with_utrs.fasta"
            SeqIO.write(annotated, output_fasta, "fasta")

            analyze_coverage(contigs_path, reads1, reads2, spades_dir)
            assembly_stats(contigs_path, spades_dir)
            print(f"[OK] {gene_id} completado.\n")
            status_lines.append(f"{gene_id}\tOK")

        except Exception as e:
            print(f"[ERROR] {gene_id}: {e}")
            status_lines.append(f"{gene_id}\tERROR\t{e}")

    with open(f"{args.output_prefix}.status_report.txt", "w") as f:
        f.write("GeneID\tStatus\tDetails\n")
        for line in status_lines:
            f.write(line + "\n")

    print(f"[INFO] Informe final guardado en: {args.output_prefix}.status_report.txt")


if __name__ == "__main__":
    main()
