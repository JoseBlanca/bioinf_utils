from functools import partial
from itertools import islice
from tempfile import NamedTemporaryFile
from subprocess import run

from rich import print

from bioinf_utils.seqs import (
    generate_random_seqs,
    create_qual_from_qual_list,
    write_fastq,
    reverse_complement,
    add_seqs,
    Seq,
)
from bioinf_utils.fastp import (
    generate_fastp_cmd,
    FastpReport,
    DEFAULT_READ_ILLUMINA_ADAPTER_SEQ_R1,
    DEFAULT_READ_ILLUMINA_ADAPTER_SEQ_R2,
)


def test_basic_run():
    seq_len = 150
    num_seqs = 100
    qual_generator = partial(create_qual_from_qual_list, quals=[40] * seq_len)
    seq_generator = generate_random_seqs(
        seq_length_distrib=seq_len, qual_generator=qual_generator
    )
    seqs = islice(seq_generator, num_seqs)
    seqs2 = islice(seq_generator, num_seqs)

    with (
        NamedTemporaryFile(suffix=".r1.fastq", mode="wt") as r1_fhand,
        NamedTemporaryFile(suffix=".r1.out.fastq") as r1_out_fhand,
        NamedTemporaryFile(suffix=".r2.fastq", mode="wt") as r2_fhand,
        NamedTemporaryFile(suffix=".r2.out.fastq") as r2_out_fhand,
        NamedTemporaryFile(suffix=".report.json") as report_fhand,
        NamedTemporaryFile(suffix=".report.json") as report2_fhand,
    ):
        write_fastq(seqs, r1_fhand)
        r1_fhand.flush()
        write_fastq(seqs2, r2_fhand)
        r2_fhand.flush()
        cmd = generate_fastp_cmd(
            r1_path=r1_fhand.name,
            out_r1_path=r1_out_fhand.name,
            report_path=report_fhand.name,
            allow_file_overwrite=True,
        )

        run(cmd, check=True, capture_output=True)
        report = FastpReport(report_fhand)

        cmd = generate_fastp_cmd(
            r1_path=r1_fhand.name,
            r2_path=r2_fhand.name,
            out_r1_path=r1_out_fhand.name,
            out_r2_path=r2_out_fhand.name,
            report_path=report2_fhand.name,
            allow_file_overwrite=True,
        )

        run(cmd, check=True, capture_output=True)
        report = FastpReport(report2_fhand)
        assert report.num_reads2_before_filtering == num_seqs
        assert report.num_reads2_after_filtering == num_seqs


def test_filter_out_bad_qual_reads():
    seq_len = 150
    num_seqs = 100
    qual_generator = partial(create_qual_from_qual_list, quals=[14] * seq_len)
    seq_generator = generate_random_seqs(
        seq_length_distrib=seq_len, qual_generator=qual_generator
    )
    seqs = islice(seq_generator, num_seqs)

    with (
        NamedTemporaryFile(suffix=".r1.fastq", mode="wt") as r1_fhand,
        NamedTemporaryFile(suffix=".r1.out.fastq") as r1_out_fhand,
        NamedTemporaryFile(suffix=".report.json") as report_fhand,
    ):
        write_fastq(seqs, r1_fhand)
        r1_fhand.flush()
        cmd = generate_fastp_cmd(
            r1_path=r1_fhand.name,
            out_r1_path=r1_out_fhand.name,
            report_path=report_fhand.name,
            allow_file_overwrite=True,
        )

        run(cmd, check=True, capture_output=True)
        report = FastpReport(report_fhand)
        assert report.num_reads1_before_filtering == num_seqs
        assert not report.num_reads1_after_filtering


def test_trim_n_nucleotides():
    seq_len = 150
    num_seqs = 100
    qual_generator = partial(create_qual_from_qual_list, quals=[40] * seq_len)
    seq_generator = generate_random_seqs(
        seq_length_distrib=seq_len, qual_generator=qual_generator
    )
    seqs = islice(seq_generator, num_seqs)

    with (
        NamedTemporaryFile(suffix=".r1.fastq", mode="wt") as r1_fhand,
        NamedTemporaryFile(suffix=".r1.out.fastq") as r1_out_fhand,
        NamedTemporaryFile(suffix=".report.json") as report_fhand,
    ):
        write_fastq(seqs, r1_fhand)
        r1_fhand.flush()
        cmd = generate_fastp_cmd(
            r1_path=r1_fhand.name,
            out_r1_path=r1_out_fhand.name,
            report_path=report_fhand.name,
            allow_file_overwrite=True,
            trim_front_n_bases=10,
            trim_tail_n_bases=5,
        )

        run(cmd, check=True, capture_output=True)
        report = FastpReport(report_fhand)
        assert report.reads_1_mean_length_after_filtering == 135


def test_trim_qual():
    seq_len = 150
    num_seqs = 100
    quals = ([40] * (seq_len - 20)) + ([10] * 20)
    qual_generator = partial(create_qual_from_qual_list, quals=quals)
    seq_generator = generate_random_seqs(
        seq_length_distrib=seq_len, qual_generator=qual_generator
    )
    seqs = islice(seq_generator, num_seqs)

    with (
        NamedTemporaryFile(suffix=".r1.fastq", mode="wt") as r1_fhand,
        NamedTemporaryFile(suffix=".r1.out.fastq") as r1_out_fhand,
        NamedTemporaryFile(suffix=".report.json") as report_fhand,
    ):
        write_fastq(seqs, r1_fhand)
        r1_fhand.flush()
        cmd = generate_fastp_cmd(
            r1_path=r1_fhand.name,
            out_r1_path=r1_out_fhand.name,
            report_path=report_fhand.name,
            allow_file_overwrite=True,
            cut_tail=True,
        )

        run(cmd, check=True, capture_output=True)
        report = FastpReport(report_fhand)
        assert report.reads_1_mean_length_after_filtering == 130


def test_illumina_adapter_trimming():
    adapter_r1 = DEFAULT_READ_ILLUMINA_ADAPTER_SEQ_R1
    adapter_r1 = Seq("adapter_r1", adapter_r1, [40] * len(adapter_r1))

    adapter_r2 = DEFAULT_READ_ILLUMINA_ADAPTER_SEQ_R2
    adapter_r2 = Seq("adapter_r2", adapter_r2, [40] * len(adapter_r2))
    adapter_len = len(adapter_r1.seq)
    assert adapter_len == len(adapter_r2.seq)
    seq_len = 150
    seq_len_no_adapt = seq_len - adapter_len
    num_seqs = 100
    qual_generator = partial(create_qual_from_qual_list, quals=[40] * seq_len_no_adapt)
    seq_generator = generate_random_seqs(
        seq_length_distrib=seq_len_no_adapt, qual_generator=qual_generator
    )
    seqs_r1 = islice(seq_generator, num_seqs)
    seqs_r1 = [add_seqs(seq, adapter_r1) for seq in seqs_r1]
    seqs_r2 = islice(seq_generator, num_seqs)
    seqs_r2 = [add_seqs(seq, adapter_r2) for seq in seqs_r2]

    with (
        NamedTemporaryFile(suffix=".r1.fastq", mode="wt") as r1_fhand,
        NamedTemporaryFile(suffix=".r1.out.fastq") as r1_out_fhand,
        NamedTemporaryFile(suffix=".r2.fastq", mode="wt") as r2_fhand,
        NamedTemporaryFile(suffix=".r2.out.fastq") as r2_out_fhand,
        NamedTemporaryFile(suffix=".report.json") as report_fhand,
        NamedTemporaryFile(suffix=".report.json") as report2_fhand,
    ):
        write_fastq(seqs_r1, r1_fhand)
        r1_fhand.flush()
        write_fastq(seqs_r2, r2_fhand)
        r2_fhand.flush()
        cmd = generate_fastp_cmd(
            r1_path=r1_fhand.name,
            out_r1_path=r1_out_fhand.name,
            report_path=report_fhand.name,
            allow_file_overwrite=True,
        )

        run(cmd, check=True, capture_output=True)
        report = FastpReport(report_fhand)

        cmd = generate_fastp_cmd(
            r1_path=r1_fhand.name,
            r2_path=r2_fhand.name,
            out_r1_path=r1_out_fhand.name,
            out_r2_path=r2_out_fhand.name,
            report_path=report2_fhand.name,
            allow_file_overwrite=True,
            read_adapter_seq_r1=DEFAULT_READ_ILLUMINA_ADAPTER_SEQ_R1,
            read_adapter_seq_r2=DEFAULT_READ_ILLUMINA_ADAPTER_SEQ_R2,
        )

        run(cmd, check=True, capture_output=True)
        report = FastpReport(report2_fhand)
        assert report.num_reads2_before_filtering == num_seqs
        assert report.num_reads2_after_filtering == num_seqs

        assert report.reads_1_mean_length_before_filtering == seq_len
        assert report.reads_2_mean_length_before_filtering == seq_len
        assert report.reads_1_mean_length_after_filtering == seq_len - adapter_len
        assert report.reads_2_mean_length_after_filtering == seq_len - adapter_len
