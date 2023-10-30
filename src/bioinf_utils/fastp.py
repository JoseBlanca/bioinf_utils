from pathlib import Path
import json
from tempfile import NamedTemporaryFile

import pandas

FASTP_BIN = "fastp"
DEFAULT_READ_ILLUMINA_ADAPTER_SEQ_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
DEFAULT_READ_ILLUMINA_ADAPTER_SEQ_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


class FastpReport:
    def __init__(self, report_fhand):
        self.report = json.load(report_fhand)

    @property
    def num_reads1_before_filtering(self):
        return self.report["read1_before_filtering"]["total_reads"]

    @property
    def num_reads2_before_filtering(self):
        try:
            result = self.report["read2_before_filtering"]["total_reads"]
        except KeyError:
            raise ValueError("No r2 reads")
        return result

    @property
    def num_reads1_after_filtering(self):
        return self.report["read1_after_filtering"]["total_reads"]

    @property
    def num_reads2_after_filtering(self):
        try:
            result = self.report["read2_after_filtering"]["total_reads"]
        except KeyError:
            raise ValueError("No r2 reads")
        return result

    @property
    def num_merged_reads(self):
        try:
            result = self.report["merged_and_filtered"]["total_reads"]
        except KeyError:
            raise ValueError("No merged reads")
        return result

    @property
    def reads1_mean_qual_per_position_before_filtering(self):
        quals = self.report["read1_before_filtering"]["quality_curves"]["mean"]
        quals = pandas.Series(quals, index=range(1, len(quals) + 1))
        return quals

    @property
    def reads2_mean_qual_per_position_before_filtering(self):
        try:
            quals = self.report["read2_before_filtering"]["quality_curves"]["mean"]
        except KeyError:
            raise ValueError("No r2 reads")
        quals = pandas.Series(quals, index=range(1, len(quals) + 1))
        return quals

    @property
    def merged_mean_qual_per_position(self):
        try:
            quals = self.report["merged_and_filtered"]["quality_curves"]["mean"]
        except KeyError:
            raise ValueError("No merged reads")
        quals = pandas.Series(quals, index=range(1, len(quals) + 1))
        return quals

    @property
    def reads1_mean_qual_per_position_after_filtering(self):
        quals = self.report["read1_after_filtering"]["quality_curves"]["mean"]
        quals = pandas.Series(quals, index=range(1, len(quals) + 1))
        return quals

    @property
    def reads2_mean_qual_per_position_after_filtering(self):
        try:
            quals = self.report["read2_after_filtering"]["quality_curves"]["mean"]
        except KeyError:
            raise ValueError("No r2 reads")
        quals = pandas.Series(quals, index=range(1, len(quals) + 1))
        return quals

    @property
    def reads1_nucleotide_content_per_position_before_filtering(self):
        content = self.report["read1_before_filtering"]["content_curves"]
        content = pandas.DataFrame(content, index=range(1, len(content["A"]) + 1))
        return content

    @property
    def reads2_nucleotide_content_per_position_before_filtering(self):
        try:
            content = self.report["read2_before_filtering"]["content_curves"]
        except KeyError:
            raise ValueError("No r2 reads")
        content = pandas.DataFrame(content, index=range(1, len(content["A"]) + 1))
        return content

    @property
    def merged_nucleotide_content_per_position(self):
        try:
            content = self.report["merged_and_filtered"]["content_curves"]
        except KeyError:
            raise ValueError("No merged reads")
        content = pandas.DataFrame(content, index=range(1, len(content["A"]) + 1))
        return content

    @property
    def reads1_nucleotide_content_per_position_after_filtering(self):
        content = self.report["read1_after_filtering"]["content_curves"]
        content = pandas.DataFrame(content, index=range(1, len(content["A"]) + 1))
        return content

    @property
    def reads2_nucleotide_content_per_position_after_filtering(self):
        try:
            content = self.report["read2_after_filtering"]["content_curves"]
        except KeyError:
            raise ValueError("No r2 reads")
        content = pandas.DataFrame(content, index=range(1, len(content["A"]) + 1))
        return content

    @property
    def reads_1_mean_length_before_filtering(self):
        return self.report["summary"]["before_filtering"]["read1_mean_length"]

    @property
    def reads_1_mean_length_after_filtering(self):
        return self.report["summary"]["after_filtering"]["read1_mean_length"]

    @property
    def reads_2_mean_length_before_filtering(self):
        try:
            return self.report["summary"]["before_filtering"]["read2_mean_length"]
        except KeyError:
            raise ValueError("No r2 reads")

    @property
    def reads_2_mean_length_after_filtering(self):
        try:
            return self.report["summary"]["after_filtering"]["read2_mean_length"]
        except KeyError:
            raise ValueError("No r2 reads")


def generate_fastp_cmd(
    r1_path,
    out_r1_path,
    report_path,
    html_report_path=None,
    fastp_bin_path=FASTP_BIN,
    num_threads=3,
    allow_file_overwrite=False,
    r2_path=None,
    out_r2_path=None,
    min_base_qual_for_good_pos=15,
    percent_bad_bases_for_filtering_out=40,
    trim_front_n_bases=0,
    trim_tail_n_bases=0,
    cut_tail=False,
    cut_front=False,
    cut_window_size=4,
    cut_mean_quality=15,
    read_adapter_seq_r1=None,
    read_adapter_seq_r2=None,
    merge_reads=False,
    overlap_len_require=30,
    overlap_diff_limit=5,
    overlap_diff_percent_limit=20,
    length_required=40,
    merged_path=None,
    overlap_correction=True,
) -> list[str]:
    r1_path = Path(r1_path)
    out_r1_path = Path(out_r1_path)
    report_path = Path(report_path)
    fastp_bin_path = Path(fastp_bin_path)
    if r2_path is not None:
        r2_path = Path(r2_path)
        out_r2_path = Path(out_r2_path)
    if merged_path is not None:
        merged_path = Path(merged_path)
    if html_report_path is not None:
        html_report_path = Path(html_report_path)

    if merge_reads:
        out_1_param = "--out1"
        out_2_param = "--out2"
    else:
        out_1_param = "-o"
        out_2_param = "-O"

    cmd = [
        str(fastp_bin_path),
        "-i",
        str(r1_path),
        out_1_param,
        str(out_r1_path),
        "--length_required",
        str(40),
    ]

    if html_report_path:
        cmd.extend(["--html", str(html_report_path)])

    if r2_path:
        cmd.extend(["-I", str(r2_path), out_2_param, str(out_r2_path)])

    if merge_reads:
        if not r2_path or not merged_path:
            raise ValueError("To merge reads, r2 path and merged_path should be given")
        cmd.extend(
            [
                "--merge",
                "--overlap_len_require",
                str(overlap_len_require),
                "--overlap_diff_limit",
                str(overlap_diff_limit),
                "--overlap_diff_percent_limit",
                str(overlap_diff_percent_limit),
                "--merged_out",
                str(merged_path),
            ]
        )
        if overlap_correction:
            cmd.append("--correction")

    if min_base_qual_for_good_pos:
        cmd.extend(["--qualified_quality_phred", str(min_base_qual_for_good_pos)])
        cmd.extend(
            [
                "--unqualified_percent_limit",
                str(percent_bad_bases_for_filtering_out),
            ]
        )
    if trim_front_n_bases:
        cmd.extend(["--trim_front1", str(trim_front_n_bases)])
    if trim_tail_n_bases:
        cmd.extend(["--trim_tail1", str(trim_tail_n_bases)])

    if cut_tail:
        cmd.append("--cut_tail")
    if cut_front:
        cmd.append("--cut_front")
    if cut_tail or cut_front:
        cmd.extend(
            [
                "--cut_window_size",
                str(cut_window_size),
                "--cut_mean_quality",
                str(cut_mean_quality),
            ]
        )

    if not (read_adapter_seq_r1 and read_adapter_seq_r2) and r2_path:
        cmd.append("--detect_adapter_for_pe")

    if read_adapter_seq_r1:
        cmd.append(f"--adapter_sequence={read_adapter_seq_r1}")
    if read_adapter_seq_r2:
        cmd.append(f"--adapter_sequence_r2={read_adapter_seq_r2}")

    if not allow_file_overwrite:
        cmd.append("--dont_overwrite")

    cmd.extend(["--json", str(report_path)])
    cmd.extend(["--thread", str(num_threads)])

    return cmd


"""
TODO


--length_required 40

--dedup
--dup_calc_accuracy 4

"""
