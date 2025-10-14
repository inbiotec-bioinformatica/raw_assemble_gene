#!/usr/bin/python3
import os
import sys
import argparse
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Align

def run_command(command):
    print(f"Running: {command}")
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        sys.exit(f"Command failed: {command}")

def get_aligned_reads(query, reads1, reads2, minimap2_args, output_prefix):
    sam_file = f"{output_prefix}.sam"
    command = ["minimap2", "-ax", "sr", *minimap2_args.split(), query, reads1, reads2]
    with open(sam_file, "w") as out:
        subprocess.run(command, stdout=out, check=True)

    fastq_file = f"{output_prefix}_reads.fastq"
    command = ["samtools", "fastq", "-f", "0x2", sam_file]
    with open(fastq_file, "w") as out:
        subprocess.run(command, stdout=out, check=True)

    if not os.path.exists(fastq_file) or os.path.getsize(fastq_file) == 0:
        raise ValueError(f"El archivo {fastq_file} no se generó o está vacío.")

    return fastq_file

def run_spades(aligned_reads, outdir):
    print(f"Archivo de entrada para SPAdes: {aligned_reads}")
    command = f"spades.py -s {aligned_reads} -o {outdir}"
    run_command(command)

    contigs_path = os.path.join(outdir, "contigs.fasta")
    if not os.path.exists(contigs_path):
        raise FileNotFoundError("ERROR: contigs.fasta no fue generado por SPAdes.")
    return contigs_path

def find_best_contig(query_seq, contigs_path):
    from Bio import Align
    from Bio import SeqIO

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    min_score = len(query_seq) * 0.75  # Umbral más flexible

    best_score = 0
    best_contig = None

    for contig in SeqIO.parse(contigs_path, "fasta"):
        score = aligner.score(contig.seq, query_seq.seq)
        if score > best_score:
            best_score = score
            best_contig = contig

    if best_score < min_score:
        print(f"[WARNING] {query_seq.id}: No se encontró un contig con similitud completa. Se utilizará el mejor fragmento parcial (score={best_score:.1f}, requerido>{min_score:.1f}).")
        if best_contig is None:
            raise ValueError(f"{query_seq.id}: No se encontró ningún contig con similitud significativa.")

    return best_contig


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import Counter
import os
import subprocess
import tempfile

def consensus_sequence(sequences):
    """Genera una secuencia consenso a partir de múltiples secuencias."""
    if not sequences:
        return ""
    max_len = max(len(seq) for seq in sequences)
    consensus = []
    for i in range(max_len):
        bases = [seq[i] for seq in sequences if i < len(seq)]
        if bases:
            most_common = Counter(bases).most_common(1)[0][0]
            consensus.append(most_common)
    return "".join(consensus)


#!/usr/bin/python3
import os
import sys
import re
import gzip
import argparse
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Align
import matplotlib.pyplot as plt
import pandas as pd

# -----------------------
# Helpers
# -----------------------
def run_command(command):
    print(f"Running: {command}")
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        sys.exit(f"Command failed: {command}")

def is_gz(path):
    return path.endswith(".gz")

def open_text(path, mode="rt"):
    if path is None:
        return None
    if is_gz(path):
        return gzip.open(path, mode)
    return open(path, mode)

def concat_reads(out_path, reads1, reads2=None):
    """Concatena reads1 y reads2 en out_path (soporta .gz)."""
    with open(out_path, "wt") as out:
        if reads1:
            with open_text(reads1, "rt") as f:
                for line in f:
                    out.write(line)
        if reads2:
            with open_text(reads2, "rt") as f:
                for line in f:
                    out.write(line)

def revcomp(s):
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return s.translate(comp)[::-1]

def parse_cigar(cigar):
    """Devuelve lista de (length,int_op) del CIGAR."""
    return [(int(l), op) for l, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]

# -----------------------
# Core pipeline functions
# -----------------------
def get_aligned_reads(query, reads1, reads2, minimap2_args, output_prefix, platform):
    """
    Alinea lecturas contra el query (gen) y extrae las lecturas alineadas.
    Devuelve path a FASTQ con lecturas alineadas (un solo archivo).
    """
    sam_file = f"{output_prefix}.sam"

    if platform == "illumina":
        cmd = ["minimap2", "-ax", "sr", *minimap2_args.split(), query, reads1, reads2]
    elif platform == "nanopore":
        cmd = ["minimap2", "-ax", "map-ont", *minimap2_args.split(), query, reads1]
    else:
        raise ValueError("Plataforma no soportada: use 'illumina' o 'nanopore'")

    with open(sam_file, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)

    fastq_file = f"{output_prefix}_reads.fastq"
    if platform == "illumina":
        cmd = ["samtools", "fastq", "-f", "0x2", sam_file]
    else:
        cmd = ["samtools", "fastq", sam_file]

    with open(fastq_file, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)

    if not os.path.exists(fastq_file) or os.path.getsize(fastq_file) == 0:
        raise ValueError(f"El archivo {fastq_file} no se generó o está vacío.")

    return fastq_file

def run_spades(aligned_reads, outdir, platform):
    print(f"Archivo de entrada para SPAdes: {aligned_reads}")
    if platform == "illumina":
        command = f"spades.py -s {aligned_reads} -o {outdir}"
    elif platform == "nanopore":
        command = f"spades.py --nanopore {aligned_reads} -o {outdir}"
    else:
        raise ValueError("Plataforma no soportada para SPAdes")

    run_command(command)

    contigs_path = os.path.join(outdir, "contigs.fasta")
    if not os.path.exists(contigs_path):
        raise FileNotFoundError("ERROR: contigs.fasta no fue generado por SPAdes.")
    return contigs_path

def find_best_contig(query_seq, contigs_path):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
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

    if best_score < min_score:
        print(f"[WARNING] {query_seq.id}: mejor fragmento parcial (score={best_score:.1f}, requerido>{min_score:.1f}).")
        if best_contig is None:
            raise ValueError(f"{query_seq.id}: No se encontró ningún contig con similitud significativa.")

    return best_contig

# -----------------------
# La función principal que extiende contig con lecturas (iterativa)
# -----------------------
def extend_contig_with_utrs(best_contig, query_seq, extension_len, output_prefix, reads1, reads2, platform="illumina"):
    """
    Extiende el gen en el contig con UTRs usando lecturas por superposición iterativa.
    - best_contig: SeqRecord del contig seleccionado
    - query_seq: SeqRecord del gen (query)
    - extension_len: long total deseada para cada UTR (5' y 3')
    - output_prefix: prefijo para archivos temporales (ej: spades_dir/geneid)
    - reads1/reads2: paths a fastq (reads2 puede ser None para nanopore)
    - platform: 'illumina' o 'nanopore'
    Retorna: SeqRecord con secuencia: left_utr (minusc) + GEN (MAYUSC) + right_utr (minusc)
    """

    import time
    from Bio.SeqRecord import SeqRecord as BSeqRecord  # import local para evitar NameError

    # presets minimap2
    if platform == "illumina":
        minimap_preset = "sr"
    elif platform == "nanopore":
        minimap_preset = "map-ont"
    else:
        raise ValueError("Plataforma desconocida")

    # secuencias base
    contig_seq = str(best_contig.seq)
    gene_seq_query = str(query_seq.seq).upper()

    # localizar gen dentro del contig (alineamiento local)
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(contig_seq, gene_seq_query)
    if len(alignments) == 0:
        raise ValueError("El gen no se alinea con el contig (no se encuentra region).")
    aln = alignments[0]
    contig_aligned_regions, gene_aligned_regions = aln.aligned
    # regiones en contig (0-based, end-exclusive)
    gene_contig_start = contig_aligned_regions[0][0]
    gene_contig_end = contig_aligned_regions[-1][1]  # end-exclusive

    # concatenar lecturas (descomprimir .gz si hace falta) -> archivo temporal
    reads_concat = f"{output_prefix}_all_reads.fastq"
    concat_reads(reads_concat, reads1, reads2)

    # variables de extensión
    left_ext = ""   # secuencia upstream (5') en orientación de referencia (5'->3')
    right_ext = ""  # secuencia downstream (3') en orientación de referencia

    max_iters = 40
    iter_count = 0
    changed = True

    # helper para merge (left y right)
    def merge_left(existing, candidate):
        """Greedy merge: candidate es prefijo nuevo (más a la izquierda). Evita duplicar solapamiento."""
        if not existing:
            return candidate
        maxk = min(len(candidate), len(existing))
        for k in range(maxk, 0, -1):
            if candidate[-k:] == existing[:k]:
                return candidate + existing[k:]
        return candidate + existing

    def merge_right(existing, candidate):
        """Greedy merge: candidate es sufijo nuevo (más a la derecha)."""
        if not existing:
            return candidate
        maxk = min(len(candidate), len(existing))
        for k in range(maxk, 0, -1):
            if existing[-k:] == candidate[:k]:
                return existing + candidate[k:]
        return existing + candidate

    # Iteratively extend left and right
    while (len(left_ext) < extension_len or len(right_ext) < extension_len) and iter_count < max_iters and changed:
        iter_count += 1
        changed = False
        # crear referencia temporal = left_ext + contig + right_ext
        ref_seq = left_ext + contig_seq + right_ext
        ref_fa = f"{output_prefix}_iter{iter_count}_ref.fasta"
        with open(ref_fa, "w") as f:
            f.write(f">ref\n{ref_seq}\n")

        # mapear lecturas contra la referencia temporal
        mapped_sam = f"{output_prefix}_iter{iter_count}.sam"
        cmd_align = [ 
                "minimap2", "-a", "-x", minimap_preset, "-p", "0.8", "--secondary=no", "-N", "100",
                ref_fa, reads_concat
            ], 
        with open (mapped_sam, "W") as out:
                subprocess.run (cmd_align, stdout=out, check=True)

        # coordenadas del gen dentro de la referencia actual
        gene_ref_start = len(left_ext) + gene_contig_start  # 0-based
        gene_ref_end_excl = len(left_ext) + gene_contig_end  # exclusive

        # anchor (primeras y ultimas 15 bases del gen, en coords ref)
        anchor_left_start = gene_ref_start
        anchor_left_end = gene_ref_start + 15 - 1
        anchor_right_start = max(gene_ref_start, gene_ref_end_excl - 15)
        anchor_right_end = gene_ref_end_excl - 1

        # track best candidate found this iteration
        best_left_candidate = ""
        best_right_candidate = ""

        # parse SAM
        with open(mapped_sam, "r") as sam:
            for line in sam:
                if line.startswith("@"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 11:
                    continue
                qname = fields[0]
                flag = int(fields[1])
                rname = fields[2]
                pos = int(fields[3])  # 1-based
                cigar = fields[5]
                seq = fields[9]
                # We only care alignments to our 'ref' record
                if rname != "ref":
                    continue

                # parse CIGAR to build aligned blocks
                ops = parse_cigar(cigar)
                ref_pos = pos - 1  # 0-based
                query_pos = 0
                blocks = []
                # detect leading/trailing soft clip lengths
                leading_softclip = ops[0][0] if ops and ops[0][1] == "S" else 0
                trailing_softclip = ops[-1][0] if ops and ops[-1][1] == "S" else 0

                # build blocks M/=/X
                for (l, op) in ops:
                    if op in ("M", "=", "X"):
                        block_ref_start = ref_pos
                        block_ref_end = ref_pos + l - 1
                        block_q_start = query_pos
                        block_q_end = query_pos + l - 1
                        blocks.append((block_ref_start, block_ref_end, block_q_start, block_q_end))
                        ref_pos += l
                        query_pos += l
                    elif op == "I" or op == "S":
                        query_pos += l
                    elif op == "D" or op == "N":
                        ref_pos += l
                    elif op == "H" or op == "P":
                        # ignore
                        pass

                # for each block, check overlap with anchors
                for i, (b_ref_start, b_ref_end, b_q_start, b_q_end) in enumerate(blocks):
                    # --- LEFT anchor: if this block overlaps the left anchor, and there are ref positions < anchor_left_start mapped by this block
                    if not (b_ref_end < anchor_left_start or b_ref_start > anchor_left_end):
                        # block overlaps left anchor
                        if b_ref_start < anchor_left_start:
                            # number of ref bases in block that are left of anchor_left_start
                            left_ref_end = min(b_ref_end, anchor_left_start - 1)
                            left_ref_len = left_ref_end - b_ref_start + 1
                            if left_ref_len > 0:
                                qstart = b_q_start
                                qend = b_q_start + left_ref_len - 1
                                candidate = ""
                                # include leading softclip if this is first aligned block and leading_softclip>0
                                if i == 0 and leading_softclip > 0:
                                    candidate = seq[0:leading_softclip] + seq[qstart:qend+1]
                                else:
                                    candidate = seq[qstart:qend+1]
                                # if read mapped to reverse strand, reverse complement candidate to put in ref orientation
                                if (flag & 16) != 0:
                                    candidate = revcomp(candidate)
                                # keep longest
                                if len(candidate) > len(best_left_candidate):
                                    best_left_candidate = candidate

                    # --- RIGHT anchor: if this block overlaps right anchor, and there are ref positions > anchor_right_end mapped by this block
                    if not (b_ref_end < anchor_right_start or b_ref_start > anchor_right_end):
                        if b_ref_end > anchor_right_end:
                            # number of ref bases in block that are right of anchor_right_end
                            right_ref_start = max(b_ref_start, anchor_right_end + 1)
                            right_ref_len = b_ref_end - right_ref_start + 1
                            if right_ref_len > 0:
                                # query positions corresponding to that ref tail
                                qoffset = right_ref_start - b_ref_start
                                qstart = b_q_start + qoffset
                                qend = qstart + right_ref_len - 1
                                candidate = ""
                                # include trailing softclip if this is last aligned block and trailing_softclip>0
                                if i == len(blocks) - 1 and trailing_softclip > 0:
                                    # trailing softclip located at end of seq
                                    candidate = seq[qstart:qend+1] + seq[len(seq) - trailing_softclip : len(seq)]
                                else:
                                    candidate = seq[qstart:qend+1]
                                if (flag & 16) != 0:
                                    candidate = revcomp(candidate)
                                if len(candidate) > len(best_right_candidate):
                                    best_right_candidate = candidate

        # now merge best candidates if they add new info
        # trim candidates to requested extension_len
        if best_left_candidate:
            best_left_candidate = best_left_candidate[-extension_len:]
            merged = merge_left(left_ext, best_left_candidate)
            # if merged is strictly longer than left_ext -> we gained bases
            if len(merged) > len(left_ext):
                # keep only up to extension_len
                left_ext = merged[-extension_len:]
                changed = True

        if best_right_candidate:
            best_right_candidate = best_right_candidate[:extension_len]
            merged_r = merge_right(right_ext, best_right_candidate)
            if len(merged_r) > len(right_ext):
                right_ext = merged_r[:extension_len]
                changed = True

        # safety: prevent infinite loops
        # if no candidate found in this iteration for both sides, break automatically (changed==False)
        # continue next iteration if changed True
        # (temporary files are left for debugging; optional: remove them)

    # Final assembled putative gene: left_ext (lowercase) + gene_seq (uppercase) + right_ext (lowercase)
    left_trim = left_ext[-extension_len:] if left_ext else ""
    right_trim = right_ext[:extension_len] if right_ext else ""
    putative_seq = left_trim.lower() + gene_seq_query.upper() + right_trim.lower()

    # Write log with extension sizes
    log_file = f"{output_prefix}_utr_extension.log"
    with open(log_file, "w") as log:
        log.write(f"Requested extension: {extension_len} bp each side\n")
        log.write(f"Left extension obtained: {len(left_trim)} bp\n")
        log.write(f"Right extension obtained: {len(right_trim)} bp\n")
        log.write(f"Iterations: {iter_count}\n")

    # Return SeqRecord
    return BSeqRecord(Seq(putative_seq), id=f"{query_seq.id}_putative", description="putative_gene_with_utrs")

# -----------------------
# Coverage & stats (unchanged / small fixes)
# -----------------------
def analyze_coverage(contigs_fasta, reads1, reads2, outdir, platform):
    print("Analizando cobertura del ensamblado...")

    if not os.path.exists(contigs_fasta):
        raise FileNotFoundError(f"No se encontró el archivo de contigs: {contigs_fasta}")

    if not reads1 or not os.path.exists(reads1):
        raise FileNotFoundError(f"No se encontró el archivo reads1: {reads1}")

    os.makedirs(outdir, exist_ok=True)
    bam_file = os.path.join(outdir, "mapped.bam")
    sorted_bam = os.path.join(outdir, "mapped.sorted.bam")
    cov_file = os.path.join(outdir, "coverage_summary.tsv")
    plot_file = os.path.join(outdir, "coverage_plot.png")
    stats_file = os.path.join(outdir, "assembly_stats.tsv")

    if platform == "illumina":
        subprocess.run(["bwa", "index", contigs_fasta], check=True)
        p1 = subprocess.Popen(["bwa", "mem", contigs_fasta, reads1, reads2], stdout=subprocess.PIPE)
    else:
        # nanopore: use minimap2 for mapping long reads
        p1 = subprocess.Popen(["minimap2", "-ax", "map-ont", contigs_fasta, reads1], stdout=subprocess.PIPE)

    with open(bam_file, "wb") as bam_out:
        p2 = subprocess.Popen(["samtools", "view", "-bS", "-"], stdin=p1.stdout, stdout=bam_out)
        p1.stdout.close()
        p2.communicate()

    subprocess.run(["samtools", "sort", bam_file, "-o", sorted_bam], check=True)
    subprocess.run(["samtools", "index", sorted_bam], check=True)

    depth_out = subprocess.check_output(["samtools", "depth", sorted_bam])
    if not depth_out:
        print("[WARNING] samtools depth returned empty output.")
        return

    depth_lines = depth_out.decode().strip().split("\n")
    data = {}
    for line in depth_lines:
        if not line:
            continue
        contig, pos, depth = line.split("\t")
        data.setdefault(contig, []).append(int(depth))

    summary = []
    for contig, depths in data.items():
        summary.append({
            "contig": contig,
            "mean_coverage": sum(depths) / len(depths),
            "max_coverage": max(depths),
            "min_coverage": min(depths)
        })
        plt.plot(range(len(depths)), depths, label=contig)

    pd.DataFrame(summary).to_csv(cov_file, sep="\t", index=False)

    with open(stats_file, "a") as f:
        f.write("\n# Cobertura por contig\n")
        f.write("contig\tmean_coverage\tmax_coverage\tmin_coverage\n")
        for row in summary:
            f.write(f"{row['contig']}\t{row['mean_coverage']:.2f}\t{row['max_coverage']}\t{row['min_coverage']}\n")

    plt.xlabel("Position")
    plt.ylabel("Coverage")
    plt.title("Contig Coverage")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_file)
    plt.close()

def assembly_stats(contigs_fasta, outdir):
    print("Calculando estadísticas del ensamblado...")
    lengths = [len(rec.seq) for rec in SeqIO.parse(contigs_fasta, "fasta")]
    if not lengths:
        raise ValueError("No se encontraron contigs en el ensamblado.")

    total_contigs = len(lengths)
    total_length = sum(lengths)
    max_length = max(lengths)
    lengths.sort(reverse=True)
    acc = 0
    half = total_length / 2
    for l in lengths:
        acc += l
        if acc >= half:
            n50 = l
            break

    stats_file = os.path.join(outdir, "assembly_stats.tsv")
    with open(stats_file, "w") as f:
        f.write("total_contigs\ttotal_length\tmax_contig_length\tN50\n")
        f.write(f"{total_contigs}\t{total_length}\t{max_length}\t{n50}\n")

    print(f"Estadísticas escritas en {stats_file}")

# -----------------------
# Main
# -----------------------
def main():
    parser = argparse.ArgumentParser(description="Asamblea y extensión de gen desde lecturas")
    parser.add_argument("--query", required=True, help="Archivo FASTA con genes objetivo")
    parser.add_argument("--reads1", required=True, help="Lecturas forward o Nanopore")
    parser.add_argument("--reads2", help="Lecturas reverse (solo Illumina)")
    parser.add_argument("--original_reads1", help="Lecturas forward originales (opcional)")
    parser.add_argument("--original_reads2", help="Lecturas reverse originales (opcional)")
    parser.add_argument("--platform", choices=["illumina","nanopore"], required=True, help="Plataforma de secuenciación")
    parser.add_argument("--minimap2_args", default="", help="Argumentos extra para minimap2")
    parser.add_argument("--extension", type=int, default=0, help="Tamaño de extensión UTR")
    parser.add_argument("--output_prefix", default="genes", help="Prefijo general para reporte final")
    args = parser.parse_args()

    original_reads1 = args.original_reads1 or args.reads1
    original_reads2 = args.original_reads2 or args.reads2

    status_lines = []

    # procesar cada gen del fasta query
    for query_seq in SeqIO.parse(args.query, "fasta"):
        gene_id = query_seq.id
        try:
            # prefijo por gen (archivos intermedios dentro de carpeta del gen)
            spades_dir = os.path.join(f"{gene_id}.output", "spades_output")
            os.makedirs(spades_dir, exist_ok=True)
            output_prefix = os.path.join(spades_dir, gene_id)

            aligned_reads = get_aligned_reads(
                args.query, args.reads1, args.reads2,
                args.minimap2_args, output_prefix, args.platform
            )

            contigs_path = run_spades(aligned_reads, spades_dir, args.platform)
            best_contig = find_best_contig(query_seq, contigs_path)

            # extender usando lecturas originales (no las submuestradas si existieran)
            annotated = extend_contig_with_utrs(
                best_contig, query_seq, args.extension,
                output_prefix, original_reads1, original_reads2, args.platform
            )

            output_fasta = f"{gene_id}.putative_gene_with_utrs.fasta"
            SeqIO.write(annotated, output_fasta, "fasta")
            print(f"[OK] {gene_id} guardado en {output_fasta}")
            status_lines.append(f"{gene_id}\tOK")

            analyze_coverage(contigs_path, original_reads1, original_reads2, spades_dir, args.platform)
            assembly_stats(contigs_path, spades_dir)

        except Exception as e:
            print(f"[ERROR] {gene_id}: {e}")
            status_lines.append(f"{gene_id}\tERROR\t{e}")

    with open(f"{args.output_prefix}.status_report.txt", "w") as f:
        f.write("GeneID\tStatus\tDetails\n")
        for line in status_lines:
            f.write(line + "\n")
    print(f"Informe final guardado en: {args.output_prefix}.status_report.txt")

if __name__ == "__main__":
    main()


