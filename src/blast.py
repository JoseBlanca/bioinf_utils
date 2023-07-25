import config

import subprocess
from pathlib import Path
import tempfile
from collections import defaultdict
import hashlib, gzip, pickle


def prepare_blast_db(fasta_file: Path, db_type, out_dir=None, skip_if_exists=False):
    fasta_file = Path(fasta_file)

    if out_dir is None:
        out_dir = fasta_file.parent
    print(fasta_file)
    print(out_dir)
    1 / 0

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
    assert blast_program in ["blastn", "blastp", "blastx", "tblastn"]

    fasta_fhand = create_fasta_file(seqs, tmp_dir)

    cmd = [
        blast_program,
        "-query",
        fasta_fhand.name,
        "-db",
        str(db_path),
        "-evalue",
        str(evalue_threshold),
        "-outfmt",
        TABBLAST_OUTFMT,
    ]

    process = subprocess.run(cmd, check=False, capture_output=True)

    lines = process.stdout.decode().splitlines()
    result = defaultdict(dict)
    for line in lines:
        items = line.strip().split()
        if len(items) == 14:
            (
                query,
                subject,
                identity,
                ali_len,
                mis,
                gap_opens,
                query_start,
                query_end,
                subject_start,
                subject_end,
                expect,
                score,
                qstrand,
                sstrand,
            ) = items
        else:
            raise RuntimeError("Wrong blast outuput")

        hsp = {
            "identity": float(identity),
            "ali_len": int(ali_len),
            "mis": int(mis),
            "gap_opens": int(gap_opens),
            "query_start": int(query_start),
            "query_end": int(query_end),
            "subject_start": int(subject_start),
            "subject_end": int(subject_end),
            "evalue": float(expect),
            "score": float(score),
            "query_strand": int(qstrand),
            "subject_strand": int(sstrand),
        }
        try:
            hsps = result[query][subject]
        except KeyError:
            hsps = []
            result[query][subject] = hsps
        hsps.append(hsp)
    return result


def filter_hsps_by_align_len(hsps, len_threshold=None):
    if len_threshold is None:
        return hsps

    filtered_hsps = []
    for hsp in hsps:
        if hsp["ali_len"] >= len_threshold:
            filtered_hsps.append(hsp)
    return filtered_hsps


def filter_hsps_by_identity(hsps, threshold=None):
    if threshold is None:
        return hsps

    filtered_hsps = []
    for hsp in hsps:
        if hsp["identity"] >= threshold:
            filtered_hsps.append(hsp)
    return filtered_hsps


def get_cdna_ids_by_blasting(
    seq, evalue_threshold=1e-20, blast_program="blastn", cache_dir=None
):
    if cache_dir:
        key = str(seq)
        key += str(evalue_threshold)
        key = hashlib.md5(str(key).encode()).hexdigest()
        cache_path = cache_dir / ("blasted_cdnas" + key + ".pickle")
        if cache_path.exists():
            return pickle.load(gzip.open(cache_path, "rb"))

    blast_db_dir = config.CACHE_DIR / "tomato_blast_db"

    if blast_program in ("blastn", "blastx"):
        db_type = "nucl"
    elif blast_program in ("blastp", "tblastn"):
        db_type = "prot"

    res = prepare_blast_db(
        config.CDNA_FASTA, out_dir=blast_db_dir, db_type=db_type, skip_if_exists=True
    )

    db_path = res["db_path"]

    seq = {"name": "cdna", "seq": seq}
    res = blast_seqs(
        [seq],
        db_path,
        blast_program=blast_program,
        tmp_dir=None,
        evalue_threshold=evalue_threshold,
    )

    cdnas = [cdna_name for cdna_name in res.get("cdna", {}).keys()]

    if cache_dir:
        pickle.dump(cdnas, gzip.open(cache_path, "wb"))

    return cdnas


if __name__ == "__main__":
    blast_db_dir = config.CACHE_DIR / "tomato_blast_db"

    res = prepare_blast_db(
        config.TOMATO_GENOME_FASTA,
        out_dir=blast_db_dir,
        db_type="nucl",
        skip_if_exists=True,
    )

    seq = "gatctcctggcagcaatggctggaaaagcttctgccattgatgtgccaggccctgaggttgatctcctggcagcaatggctggaaaatacaaggtgtacttggtgatgggtgtaattgagagagatggatacacgctatattgcacatacaaggtgtacttggtgatgggtgtaattgagagagatggatacacgctatattgcacatacaaggtgtacttggtgatgggtgtaattgagagagatggatacacgctatattgcacagtgcttttcttcgactctcagggtcactaccttgggaagcatcggaagataatgccaacagtgcttttcttcgactctcagggtcactaccttgggaagcatcggaagataatgccaacagtgcttttcttcgactctcagggtcactaccttgggaagcatcggaagataatgccaacagcgttagagcggataatctggggttttggggatggatcaacaattccagtttatgacactgcgttagagcggataatctggggttttggggatggatcaacaattccagtttatgacactgcgttagagcggataatctggggttttggggatggatcaacaattccagtttatgacactcctgttggaaaaataggtgctgcaatatgttgggagaacagaatgccacttctaaggacccctgttggaaaaataggtgctgcaatatgttgggagaacagaatgccacttctaaggacccctgttggaaaaataggtgctgcaatatgttgggagaacagaatgccacttctaaggaccgcaatgtatgctaaaggcattgagatatattgtgcacctacagctgatgctagggaagtggcaatgtatgctaaaggcattgagatatattgtgcacctacagctgatgctagggaagtggcaatgtatgctaaaggcattgagatatattgtgcacctacagctgatgctagggaagtgtgg"
    cdnas = get_cdna_ids_by_blasting(seq, cache_dir=config.CACHE_DIR)

    print(cdnas)

    cdnas = get_cdna_ids_by_blasting(seq, cache_dir=config.CACHE_DIR)

    print(cdnas)

    seq = {
        "name": "seq1",
        "seq": "AGACAAGTGGTGAAGAAKAAGATGATATGCAGCAATGCATTTCACCACTTTATATAGCATGGAGTGGATTTCTCCACCTCATTTAATAGTATGAAGTGGAGGCAGCCCCCCTCTACACCTGTCCACTAAGGCCAGCCCACAATCTGATCCCTTTTAATTTTTGCCTTGAGTGGTGGGGCCCATTGGATTAAATCAATCCAAATTAGCCAC",
    }
    res = blast_seqs([seq], res["db_path"], "blastn")
