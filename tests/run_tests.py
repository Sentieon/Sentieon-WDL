#!/usr/bin/env python3

import json
import os
import pytest
import subprocess
import tempfile

base_cmd = "java -Dconfig.file={cromwell_config} -jar {cromwell} run -i {inputs_json} {wdl}"

base_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
wdl_files = {
    "germline_fastq": f"{base_dir}/pipelines/germline_DNAseq-DNAscope/germline_fastqToVcf.wdl",
    "germline_fastq_split": f"{base_dir}/pipelines/germline_DNAseq-DNAscope/germline_fastqToVcf_split.wdl",
    "germline_bam_split": f"{base_dir}/pipelines/germline_DNAseq-DNAscope/germline_BamToVcf_split.wdl",
}

input_jsons = {
    "germline_fastq": f"{base_dir}/pipelines/germline_DNAseq-DNAscope/germline_fastqToVcf.example_inputs.json"
}

QUICKSTART_URL = "https://s3.amazonaws.com/sentieon-release/other/sentieon_quickstart.tar.gz"
DNASCOPE_MODEL_URL = "https://s3.amazonaws.com/sentieon-release/other/SentieonDNAscopeModel1.1.model"

TMP_DIR = os.environ.get("SENTIEON_TMPDIR") or os.environ.get("TMPDIR") or "/tmp"

def _download_quickstart():
    cmd = f"curl -L {QUICKSTART_URL} | tar -C {TMP_DIR} -zxf -"
    if not os.path.isfile(f"{TMP_DIR}/sentieon_quickstart/1.fastq.gz"):
        subprocess.run(cmd, shell=True, check=True)

@pytest.fixture
def downloaded_quickstart():
    _download_quickstart()
    d = {
        "sentieon_germline.r1_fastq": [f"{TMP_DIR}/sentieon_quickstart/1.fastq.gz"],
        "sentieon_germline.r2_fastq": [f"{TMP_DIR}/sentieon_quickstart/2.fastq.gz"],
        "sentieon_germline.read_groups": ["@RG\\tID:sample-1\\tSM:sample\\tPL:ILLUMINA"],
        "sentieon_germline.ref_fasta": f"{TMP_DIR}/sentieon_quickstart/reference/ucsc.hg19_chr22.fasta",
        "sentieon_germline.dbsnp_vcf": f"{TMP_DIR}/sentieon_quickstart/reference/dbsnp_135.hg19_chr22.vcf",
        "sentieon_germline.dbsnp_vcf_tbi": f"{TMP_DIR}/sentieon_quickstart/reference/dbsnp_135.hg19_chr22.vcf.idx",
    }
    for suffix in ("fai", "amb", "ann", "amb", "bwt", "pac", "sa"):
        d[f"sentieon_germline.ref_{suffix}"] = f"{TMP_DIR}/sentieon_quickstart/reference/ucsc.hg19_chr22.fasta.{suffix}"
    vcfs = [
        "1000G_phase1.snps.high_confidence.hg19_chr22.sites.vcf",
        "dbsnp_135.hg19_chr22.vcf",
        "Mills_and_1000G_gold_standard.indels.hg19_chr22.sites.vcf",
    ]
    d["sentieon_germline.bqsr_vcfs"] = [f"{TMP_DIR}/sentieon_quickstart/reference/" + v for v in vcfs]
    d["sentieon_germline.bqsr_vcf_tbis"] = [f"{TMP_DIR}/sentieon_quickstart/reference/" + v + ".idx"for v in vcfs]
    return d

@pytest.fixture
def quickstart_bam():
    _download_quickstart()
    sorted_bam = f"{TMP_DIR}/sentieon_quickstart/sorted.bam"
    recal_table = f"{TMP_DIR}/sentieon_quickstart/recal.table"

    cmd = "sentieon bwa mem -R '@RG\\tID:sample\\tSM:sample' -K 10000000 -t 32 {ref} {fq1} {fq2} | sentieon util sort -t 32 -o {sorted_bam} -i - --sam2bam"
    bwa_cmd = cmd.format(
        ref = f"{TMP_DIR}/sentieon_quickstart/reference/ucsc.hg19_chr22.fasta",
        fq1 = f"{TMP_DIR}/sentieon_quickstart/1.fastq.gz",
        fq2 = f"{TMP_DIR}/sentieon_quickstart/2.fastq.gz",
        sorted_bam = sorted_bam,
    )
    cmd = "sentieon driver -i {sorted_bam} -r {ref} --algo QualCal {recal_table}"
    qcal_cmd = cmd.format(
        sorted_bam = sorted_bam,
        ref = f"{TMP_DIR}/sentieon_quickstart/reference/ucsc.hg19_chr22.fasta",
        recal_table = recal_table,
    )

    if not os.path.isfile(f"{TMP_DIR}/sentieon_quickstart/sorted.bam"):
        subprocess.run(bwa_cmd, shell=True, check=True)
        subprocess.run(qcal_cmd, shell=True, check=True)
    d = {
        "sentieon_germline.bam": [sorted_bam],
        "sentieon_germline.bam_index": [sorted_bam + ".bai"],
        "sentieon_germline.bqsr_table": recal_table,
    }
    return d

@pytest.fixture
def downloaded_model():
    TMP_DIR = os.environ.get("SENTIEON_TMPDIR") or os.environ.get("TMPDIR") or "/tmp"

    cmd = f"mkdir -p {TMP_DIR}/dnascope_models/illumina_wgs; curl -o {TMP_DIR}/dnascope_models/illumina_wgs/SentieonDNAscopeModel1.1.model -L {DNASCOPE_MODEL_URL}"
    if not os.path.isfile(f"{TMP_DIR}/dnascope_models/illumina_wgs/SentieonDNAscopeModel1.1.model"):
        subprocess.run(cmd, shell=True, check=True)

    d = {
        "sentieon_germline.calling_algo": "DNAscope",
        "sentieon_germline.dnascope_model": f"{TMP_DIR}/dnascope_models/illumina_wgs/SentieonDNAscopeModel1.1.model",
    }
    return d

def set_licsrvr_docker():
    d = {
        "sentieon_germline.sentieon_license_server": os.environ["SENTIEON_LICENSE"],
        "sentieon_germline.sentieon_docker": "sentieon/sentieon-wdl:latest",
    }
    return d

def gen_command(wdl, json_fn):
    d = {
        "wdl": wdl,
        "inputs_json": json_fn,
        "cromwell_config": os.environ["CROMWELL_CONFIG"],
        "cromwell": os.environ["CROMWELL_JAR"],
    }
    return base_cmd.format(**d)


def test_germline_fastq(downloaded_quickstart):
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as test_json:
        test_json.close()
        d = {**downloaded_quickstart, **set_licsrvr_docker()}
        json.dump(d, open(test_json.name, 'w'))

        cmd = gen_command(wdl_files["germline_fastq"], test_json.name)
        subprocess.run(cmd, shell=True, check=True)

        cmd = gen_command(wdl_files["germline_fastq_split"], test_json.name)
        subprocess.run(cmd, shell=True, check=True)
        os.unlink(test_json.name)

def test_germline_fastq_gvcf(downloaded_quickstart):
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as test_json:
        test_json.close()
        d = {**downloaded_quickstart, **set_licsrvr_docker()}
        d["sentieon_germline.output_gvcf"] = True
        json.dump(d, open(test_json.name, 'w'))

        cmd = gen_command(wdl_files["germline_fastq"], test_json.name)
        subprocess.run(cmd, shell=True, check=True)

        cmd = gen_command(wdl_files["germline_fastq_split"], test_json.name)
        subprocess.run(cmd, shell=True, check=True)
        os.unlink(test_json.name)

def test_germline_fastq_dnascope(downloaded_quickstart, downloaded_model):
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as test_json:
        test_json.close()
        d = {**downloaded_quickstart, **downloaded_model, **set_licsrvr_docker()}
        d["sentieon_germline.output_gvcf"] = True
        json.dump(d, open(test_json.name, 'w'))

        cmd = gen_command(wdl_files["germline_fastq"], test_json.name)
        subprocess.run(cmd, shell=True, check=True)

        cmd = gen_command(wdl_files["germline_fastq_split"], test_json.name)
        subprocess.run(cmd, shell=True, check=True)
        os.unlink(test_json.name)

def test_germline_bam(downloaded_quickstart, quickstart_bam, downloaded_model):
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as test_json, tempfile.NamedTemporaryFile(mode='w', suffix=".bed", delete=False) as test_bed:
        test_json.close()

        print("chr22	10000000	40000000", file=test_bed)
        test_bed.close()

        d = {
            **downloaded_quickstart,
            **set_licsrvr_docker(),
            **quickstart_bam,
            **downloaded_model,
        }
        d["sentieon_germline.bqsr_intervals"] = test_bed.name
        d["sentieon_germline.calling_intervals"] = test_bed.name

        for k in [
            "sentieon_germline.read_groups",
            "sentieon_germline.r1_fastq",
            "sentieon_germline.r2_fastq"
        ]:
            del d[k]
        json.dump(d, open(test_json.name, 'w'))

        cmd = gen_command(wdl_files["germline_bam_split"], test_json.name)
        subprocess.run(cmd, shell=True, check=True)
        os.unlink(test_json.name)
        os.unlink(test_bed.name)

def test_germline_bam_dnascope(downloaded_quickstart, quickstart_bam):
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as test_json:
        test_json.close()
        d = {**downloaded_quickstart, **set_licsrvr_docker(), **quickstart_bam}
        d["sentieon_germline.output_gvcf"] = True
        d["sentieon_germline.realign_input"] = True
        d["sentieon_germline.run_dedup_and_qc"] = True
        d["sentieon_germline.run_bqsr"] = True

        for k in [
            "sentieon_germline.read_groups",
            "sentieon_germline.r1_fastq",
            "sentieon_germline.r2_fastq"
        ]:
            del d[k]
        json.dump(d, open(test_json.name, 'w'))

        cmd = gen_command(wdl_files["germline_bam_split"], test_json.name)
        subprocess.run(cmd, shell=True, check=True)
        os.unlink(test_json.name)
