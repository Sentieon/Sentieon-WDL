version 1.0

task GermlineCalling {
  input {
    # Input files
    Array[File] aligned_reads
    Array[File] aligned_index
    File? bqsr_table

    # Reference genome files
    File ref_fasta
    File ref_fai

    # Known sites VCFs
    File? dbsnp_vcf
    File? dbsnp_vcf_tbi

    # Interval files - Both are optional but recommended
    File? calling_intervals

    # Workflow arguments
    String sample_name = "sample"
    String calling_algo = "Haplotyper"
    Boolean output_gvcf = false
    String calling_driver_xargs = ""
    String calling_algo_xargs = "--pcr_indel_model none"
    File? dnascope_model

    # Sentieon license configuration
    File? sentieon_license_file
    String? sentieon_license_server

    # Execution
    String n_threads = "32"
    String memory = "10 GB"
    Int preemptible_tries = 3
    Int disk_size = 50
    String sentieon_docker

    # Expressions
    String vcf_suffix = if output_gvcf then ".g.vcf.gz" else ".vcf.gz"
    String vcf_index = if output_gvcf then ".g.vcf.gz.tbi" else ".vcf.gz.tbi"

    # Disk size
    Int job_disk_size = disk_size
  }

  command <<<
    set -exvuo pipefail

    # License configuration
    license_file=~{default="" sentieon_license_file}
    license_server=~{default="" sentieon_license_server}
    if [[ -n "$license_file" ]]; then
      export SENTIEON_LICENSE=~{sentieon_license_file}
    elif [[ -n "$license_server" ]]; then
      export SENTIEON_LICENSE=~{sentieon_license_server}
    else
      echo "Error: Please supply either a license file or a license server address"
      exit 1
    fi

    # Sanity check the calling algo
    calling_algo=~{calling_algo}
    if ! [[ "$calling_algo" =~ ^(Haplotyper|DNAscope|Genotyper)$ ]]; then
      echo "Supported calling algos are Haplotyper, DNAscope, and Genotyper"
      exit 1
    fi

    # Ensure each input is adjecent to it's index
    aln_files=("~{sep='" "' aligned_reads}")
    aln_idxs=("~{sep='" "' aligned_index}")
    input_str=""
    for i in $(seq 1 ${#aln_files[@]}); do
      i=$((i - 1))
      aln="${aln_files[$i]}"
      aln_suffix=".bam"
      idx_suffix=".bam.bai"
      if [[ $aln == *.cram ]]; then
        aln_suffix=".cram"
        idx_suffix=".cram.crai"
      fi
      ln -s "${aln_files[$i]}" ./input_"$i"."$aln_suffix"
      ln -s "${aln_idxs[$i]}" ./input_"$i"."$idx_suffix"
      input_str="$input_str -i ./input_$i.$aln_suffix"
    done

    # Variant calling
    calling_intervals="~{default='' calling_intervals}"
    dbsnp_vcf="~{default='' dbsnp_vcf}"
    bqsr_table="~{default='' bqsr_table}"
    dnascope_model=""
    if [[ "~{calling_algo}" == "DNAscope" ]]; then
      dnascope_model="~{default="" dnascope_model}"
    fi
    sentieon driver -r ~{ref_fasta} $input_str \
      ${bqsr_table:+-q "$bqsr_table"} \
      ${calling_intervals:+--interval "$calling_intervals"} \
      ~{calling_driver_xargs} \
      --algo ~{calling_algo} \
      ~{true="--emit_mode gvcf" false="" output_gvcf} \
      ${dbsnp_vcf:+-d "$dbsnp_vcf"} \
      ${dnascope_model:+--model "$dnascope_model"} \
      ~{calling_algo_xargs} \
      "~{sample_name}.~{true="g.vcf.gz" false="vcf.gz" output_gvcf}"
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: sentieon_docker
    memory: memory
    cpu: n_threads
    disks: job_disk_size  # Disk usage should be ~1.5x the size of the input BAM
  }
  output {
    # VCF files
    File calls_vcf = "~{sample_name}~{vcf_suffix}"
    File calls_vcf_tbi = "~{sample_name}~{vcf_index}"
  }
}

task DNAModelApply {
  input {
    # Input files
    File vcf
    File vcf_tbi

    # Reference genome files
    File ref_fasta
    File ref_fai

    # Workflow arguments
    String sample_name = "sample"
    Boolean output_gvcf = false
    File? dnascope_model

    # Sentieon license configuration
    File? sentieon_license_file
    String? sentieon_license_server

    # Execution
    String n_threads = "32"
    String memory = "10 GB"
    Int preemptible_tries = 2
    Int disk_size = 50
    String sentieon_docker

    # Expressions
    String vcf_suffix = if output_gvcf then ".g.vcf.gz" else ".vcf.gz"
    String vcf_index = if output_gvcf then ".g.vcf.gz.tbi" else ".vcf.gz.tbi"

    # Disk size
    Int job_disk_size = disk_size
  }

  command <<<
    set -exvuo pipefail

    # License configuration
    license_file=~{default="" sentieon_license_file}
    license_server=~{default="" sentieon_license_server}
    if [[ -n "$license_file" ]]; then
      export SENTIEON_LICENSE=~{sentieon_license_file}
    elif [[ -n "$license_server" ]]; then
      export SENTIEON_LICENSE=~{sentieon_license_server}
    else
      echo "Error: Please supply either a license file or a license server address"
      exit 1
    fi

    # Run DNAModelApply
    sentieon driver -r ~{ref_fasta} \
      --algo DNAModelApply --model "~{default='' dnascope_model}" -v "~{vcf}" \
      "~{sample_name}.~{true="g.vcf.gz" false="vcf.gz" output_gvcf}"
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: sentieon_docker
    memory: memory
    cpu: n_threads
    disks: job_disk_size  # Disk usage should be ~3x the size of the input fastq
  }
  output {
    # VCF files
    File calls_vcf = "~{sample_name}~{vcf_suffix}"
    File calls_vcf_tbi = "~{sample_name}~{vcf_index}"
  }
}
