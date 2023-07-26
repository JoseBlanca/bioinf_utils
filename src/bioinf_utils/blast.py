import subprocess
from pathlib import Path
import tempfile
from collections import defaultdict

import pandas


def prepare_blast_db(fasta_file: Path, db_type, out_dir=None, skip_if_exists=False):
    fasta_file = Path(fasta_file)

    if out_dir is None:
        out_dir = fasta_file.parent

    fname = fasta_file.stem
    out_base_path = out_dir / fname

    if db_type == "nucl":
        out_path = Path(str(out_base_path) + ".nhr")
    if db_type == "prot":
        out_path = Path(str(out_base_path) + ".phr")

    if skip_if_exists and out_path.exists():
        return {"db_path": out_base_path}

    cmd = [
        "makeblastdb",
        "-in",
        str(fasta_file),
        "-out",
        str(out_base_path),
        "-dbtype",
        db_type,
    ]
    subprocess.run(cmd, check=True)
    return {"db_path": out_base_path}


def create_fasta_file(seqs, tmp_dir):
    fasta_fhand = tempfile.NamedTemporaryFile(suffix=".fasta", dir=tmp_dir, mode="wt")

    for seq in seqs:
        fasta_fhand.write(f'>{seq["name"]}\n{seq["seq"]}\n')

    fasta_fhand.flush()
    return fasta_fhand


TABBLAST_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe"


def blast_seqs(seqs, db_path, blast_program, tmp_dir=None, evalue_threshold=1e-5):
    fasta_fhand = create_fasta_file(seqs, tmp_dir)

    return blast_file(
        fasta_path=fasta_fhand.name,
        db_path=db_path,
        blast_program=blast_program,
        evalue_threshold=evalue_threshold,
    )


COLUMN_MAPPING = {
    "qseqid": "query_id",
    "sseqid": "subject_id",
    "pident": "percent_ident",
    "length": "alignment_length",
    "mismatch": "num_mismatches",
    "gapopen": "num_gapopenings",
    "qstart": "query_start",
    "qend": "query_end",
    "sstart": "subject_start",
    "send": "subject_end",
    "evalue": "evalue",
    "bitscore": "bitscore",
    "qframe": "query_frame",
    "sframe": "subject_frame",
}


def blast_file(
    fasta_path,
    db_path,
    blast_program,
    evalue_threshold=1e-5,
    out_fmt=TABBLAST_OUTFMT,
    tmp_dir=None,
):
    assert blast_program in ["blastn", "blastp", "blastx", "tblastn"]

    with tempfile.NamedTemporaryFile(suffix=".blast.csv", dir=tmp_dir) as blast_fhand:
        cmd = [
            blast_program,
            "-query",
            str(fasta_path),
            "-db",
            str(db_path),
            "-evalue",
            str(evalue_threshold),
            "-outfmt",
            out_fmt,
            "-out",
            blast_fhand.name,
        ]

        subprocess.run(cmd, check=True, capture_output=True)

        blast_result = pandas.read_csv(blast_fhand.name, sep="\t")
        columns = [COLUMN_MAPPING.get(col, col) for col in TABBLAST_OUTFMT.split()[1:]]
        blast_result.columns = columns
    return blast_result
