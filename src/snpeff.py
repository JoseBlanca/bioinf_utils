import config

from pathlib import Path
import os
import subprocess

import utils_subprocess

SNPEFF_JAR = config.SNPEFF_JAR

HIGH_IMPACT = "HIGH"
MODERATE_IMPACT = "MODERATE"
LOW_IMPACT = "LOW"
HIGH_OR_MODERATE_IMPACTS = (HIGH_IMPACT, MODERATE_IMPACT)


def create_snpeff_db(
    genome_version: str,
    genome_species: str,
    genome_fasta_path: Path,
    genome_gff_path: Path,
    genome_protein_fasta_path,
    genome_cds_fasta_path,
    snpeff_dir: Path,
    snpeff_config: Path,
):
    data_dir = snpeff_dir / "data"
    if not snpeff_dir.exists():
        snpeff_dir.mkdir(exist_ok=True)
        data_dir.mkdir()

    if not snpeff_config.exists():
        fhand = snpeff_config.open("wt")
        fhand.close()

    with snpeff_config.open("rt") as fhand:
        config_content = fhand.read()

    genome_line = f"{genome_version}.genome : {genome_species}\n"
    if genome_line not in config_content:
        with snpeff_config.open("wt") as fhand:
            fhand.write("\n")
            fhand.write(genome_line)
            fhand.flush()
    genome_dir = data_dir / genome_version
    genome_dir.mkdir(exist_ok=True)

    snpeff_gff_path = genome_dir / "genes.gff"
    if not snpeff_gff_path.exists():
        os.symlink(genome_gff_path, snpeff_gff_path)
    snpeff_fasta_path = genome_dir / f"sequences.fa"
    if not snpeff_fasta_path.exists():
        os.symlink(genome_fasta_path, snpeff_fasta_path)
    protein_fasta_path = genome_dir / f"protein.fa"
    if not protein_fasta_path.exists():
        os.symlink(genome_protein_fasta_path, protein_fasta_path)
    cds_fasta_path = genome_dir / f"cds.fa"
    if not cds_fasta_path.exists():
        os.symlink(genome_cds_fasta_path, cds_fasta_path)

    cmd = ["java", "-jar", str(SNPEFF_JAR), "build", "-gff3", "-v", genome_version]
    subprocess.run(cmd, cwd=snpeff_dir, check=True)


def annotate_vcf(
    in_vcf, out_vcf, genome_species_version, config_path, cwd=None, capture_output=False
):
    cmd = [
        "java",
        "-jar",
        str(SNPEFF_JAR),
        "-c",
        str(config_path),
        "-v",
        genome_species_version,
        str(in_vcf),
        "-noStats",
        "-noLog",
        "|",
        "bgzip" ">",
        str(out_vcf),
    ]
    utils_subprocess.run_in_sh_script(cmd, cwd=cwd, capture_output=capture_output)


if __name__ == "__main__":
    create_snpeff_db(
        genome_fasta_path=config.GENOME_TOMATO_FASTA,
        genome_gff_path=config.GENOME_TOMATO_GFF,
        genome_protein_fasta_path=config.GENOME_TOMATO_PROTEIN_FASTA,
        genome_cds_fasta_path=config.GENOME_TOMATO_CDS_FASTA,
        snpeff_dir=config.SNPEFF_DIR,
        snpeff_config=config.SNPEFF_CONFIG,
        genome_version=config.GENOME_TOMATO_SNPEFF_VERSION,
        genome_species=config.GENOME_TOMATO_SNPEFF_SPECIES,
    )
