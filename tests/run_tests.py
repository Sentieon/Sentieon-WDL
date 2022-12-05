#!/usr/bin/env python3

import json
import os
import pytest
import subprocess
import tempfile

base_cmd = "java -Dconfig.file={cromwell_config} -jar {cromwell} run -i {inputs_json} {wdl}"

base_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
wdl_files = {
    "germline_fastq": f"{base_dir}/pipelines/germline_fastqToVcf.wdl",
    "germline_fastq_split": f"{base_dir}/pipelines/germline_fastqToVcf_split.wdl",
}

input_jsons = {
    "germline_fastq": f"{base_dir}/pipelines/germline_fastqToVcf.example_inputs.json"
}

QUICKSTART_URL = "https://s3.amazonaws.com/sentieon-release/other/sentieon_quickstart.tar.gz"
DNASCOPE_MODEL_URL = "https://s3.amazonaws.com/sentieon-release/other/SentieonDNAscopeModel1.1.model"

@pytest.fixture
def downloaded_quickstart():
    tmp_dir = os.environ.get("SENTIEON_TMPDIR") or os.environ.get("TMPDIR") or "/tmp"
    cmd = f"curl -L {QUICKSTART_URL} | tar -C {tmp_dir} -zxf -"
    if not os.path.isfile(f"{tmp_dir}/sentieon_quickstart/1.fastq.gz"):
        subprocess.run(cmd, shell=True, check=True)

    d = {
        "sentieon_germline.r1_fastq": [f"{tmp_dir}/sentieon_quickstart/1.fastq.gz"],
        "sentieon_germline.r2_fastq": [f"{tmp_dir}/sentieon_quickstart/2.fastq.gz"],
        "sentieon_germline.read_groups": ["@RG\\tID:sample-1\\tSM:sample\\tPL:ILLUMINA"],
        "sentieon_germline.ref_fasta": f"{tmp_dir}/sentieon_quickstart/reference/ucsc.hg19_chr22.fasta",
        "sentieon_germline.dbsnp_vcf": f"{tmp_dir}/sentieon_quickstart/reference/dbsnp_135.hg19_chr22.vcf",
        "sentieon_germline.dbsnp_vcf_tbi": f"{tmp_dir}/sentieon_quickstart/reference/dbsnp_135.hg19_chr22.vcf.idx",
    }
    for suffix in ("fai", "amb", "ann", "amb", "bwt", "pac", "sa"):
        d[f"sentieon_germline.ref_{suffix}"] = f"{tmp_dir}/sentieon_quickstart/reference/ucsc.hg19_chr22.fasta.{suffix}"
    vcfs = [
        "1000G_phase1.snps.high_confidence.hg19_chr22.sites.vcf",
        "dbsnp_135.hg19_chr22.vcf",
        "Mills_and_1000G_gold_standard.indels.hg19_chr22.sites.vcf",
    ]
    d["sentieon_germline.bqsr_vcfs"] = [f"{tmp_dir}/sentieon_quickstart/reference/" + v for v in vcfs]
    d["sentieon_germline.bqsr_vcf_tbis"] = [f"{tmp_dir}/sentieon_quickstart/reference/" + v + ".idx"for v in vcfs]
    return d

@pytest.fixture
def downloaded_model():
    tmp_dir = os.environ.get("SENTIEON_TMPDIR") or os.environ.get("TMPDIR") or "/tmp"

    cmd = f"mkdir -p {tmp_dir}/dnascope_models/illumina_wgs; curl -o {tmp_dir}/dnascope_models/illumina_wgs/SentieonDNAscopeModel1.1.model -L {DNASCOPE_MODEL_URL}"
    if not os.path.isfile(f"{tmp_dir}/dnascope_models/illumina_wgs/SentieonDNAscopeModel1.1.model"):
        subprocess.run(cmd, shell=True, check=True)

    d = {
        "sentieon_germline.calling_algo": "DNAscope",
        "sentieon_germline.dnascope_model": f"{tmp_dir}/dnascope_models/illumina_wgs/SentieonDNAscopeModel1.1.model",
    }
    return d

def test_germline_fastq(downloaded_quickstart):
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as test_json:
        test_json.close()
        d = downloaded_quickstart
        d["sentieon_germline.sentieon_license_server"] = os.environ["SENTIEON_LICENSE"]
        d["sentieon_germline.sentieon_docker"] = "wdl-docker"
        json.dump(d, open(test_json.name, 'w'))

        d = {
            "wdl": wdl_files["germline_fastq"],
            "inputs_json": test_json.name,
            "cromwell_config": os.environ["CROMWELL_CONFIG"],
            "cromwell": os.environ["CROMWELL_JAR"],
        }
        cmd = base_cmd.format(**d)
        subprocess.run(cmd, shell=True, check=True)

        d.update({"wdl": wdl_files["germline_fastq_split"]})
        cmd = base_cmd.format(**d)
        subprocess.run(cmd, shell=True, check=True)
        os.unlink(test_json.name)

def test_germline_fastq_gvcf(downloaded_quickstart):
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as test_json:
        test_json.close()
        d = downloaded_quickstart
        d["sentieon_germline.sentieon_license_server"] = os.environ["SENTIEON_LICENSE"]
        d["sentieon_germline.sentieon_docker"] = "wdl-docker"
        d["sentieon_germline.output_gvcf"] = True
        json.dump(d, open(test_json.name, 'w'))

        d = {
            "wdl": wdl_files["germline_fastq"],
            "inputs_json": test_json.name,
            "cromwell_config": os.environ["CROMWELL_CONFIG"],
            "cromwell": os.environ["CROMWELL_JAR"],
        }
        cmd = base_cmd.format(**d)
        subprocess.run(cmd, shell=True, check=True)

        d.update({"wdl": wdl_files["germline_fastq_split"]})
        cmd = base_cmd.format(**d)
        subprocess.run(cmd, shell=True, check=True)
        os.unlink(test_json.name)

def test_germline_fastq_dnascope(downloaded_quickstart, downloaded_model):
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as test_json:
        test_json.close()
        d = downloaded_quickstart
        d.update(downloaded_model)
        d["sentieon_germline.sentieon_license_server"] = os.environ["SENTIEON_LICENSE"]
        d["sentieon_germline.sentieon_docker"] = "wdl-docker"
        d["sentieon_germline.output_gvcf"] = True
        json.dump(d, open(test_json.name, 'w'))

        d = {
            "wdl": wdl_files["germline_fastq"],
            "inputs_json": test_json.name,
            "cromwell_config": os.environ["CROMWELL_CONFIG"],
            "cromwell": os.environ["CROMWELL_JAR"],
        }
        cmd = base_cmd.format(**d)
        subprocess.run(cmd, shell=True, check=True)

        d.update({"wdl": wdl_files["germline_fastq_split"]})
        cmd = base_cmd.format(**d)
        subprocess.run(cmd, shell=True, check=True)
        os.unlink(test_json.name)
