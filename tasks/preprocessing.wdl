version 1.0

task DedupAndQc {
  input {
    # Input files
    Array[File] bams
    Array[File] bams_idx

    # Reference genome files
    File ref_fasta
    File ref_fai

    # Workflow arguments
    String sample_name = "sample"
    Boolean output_cram = true
    String lc_xargs = ""
    String dedup_xargs = "--cram_write_options version=3.0,compressor=rans"
    String alnstat_adapter_seq = ""  # The adapter sequence for the AlignmentStat algo

    # Sentieon license configuration
    File? sentieon_license_file
    String? sentieon_license_server

    # Execution
    String n_threads = "16"
    String memory = "10 GB"
    Int preemptible_tries = 2
    Int disk_size = 50
    String sentieon_docker

    # Expressions
    String alignment_suffix = if output_cram then "cram" else "bam"
    String alignment_index = if output_cram then "crai" else "bai"

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

    # Ensure each BAM is adjecent to it's index
    bam_files=("~{sep='" "' bams}")
    bai_files=("~{sep='" "' bams_idx}")
    bam_input_str=""
    for i in $(seq 1 ${#bam_files[@]}); do
      i=$((i - 1))
      ln -s "${bam_files[$i]}" ./input_bam_"$i".bam
      ln -s "${bai_files[$i]}" ./input_bam_"$i".bam.bai
      bam_input_str="$bam_input_str -i ./input_bam_$i.bam"
    done

    # Dedup
    # Perform duplicate marking and QC
    sentieon driver $bam_input_str -r ~{ref_fasta} \
      --algo LocusCollector ~{lc_xargs} "~{sample_name}_score.txt.gz" \
      --algo MeanQualityByCycle "~{sample_name}_mq_metrics.txt" \
      --algo QualDistribution "~{sample_name}_qd_metrics.txt" \
      --algo GCBias --summary "~{sample_name}_gc_summary.txt" "~{sample_name}_gc_metrics.txt" \
      --algo AlignmentStat --adapter_seq '~{alnstat_adapter_seq}' "~{sample_name}_aln_metrics.txt" \
      --algo InsertSizeMetricAlgo "~{sample_name}_is_metrics.txt"

    sentieon driver $bam_input_str -r ~{ref_fasta} --algo Dedup \
      ~{dedup_xargs} --score_info "~{sample_name}_score.txt.gz" \
      --metrics "~{sample_name}_dedup_metrics.txt" \
      "~{sample_name}_aligned.~{true="cram" false="bam" output_cram}"

    # Plot the metrics output
    sentieon plot GCBias -o "~{sample_name}_gc-report.pdf" "~{sample_name}_gc_metrics.txt" &
    sentieon plot QualDistribution -o "~{sample_name}_qd-report.pdf" "~{sample_name}_qd_metrics.txt" &
    sentieon plot MeanQualityByCycle -o "~{sample_name}_mq-report.pdf" "~{sample_name}_mq_metrics.txt" &
    sentieon plot InsertSizeMetricAlgo -o "~{sample_name}_is-report.pdf" "~{sample_name}_is_metrics.txt" &
    wait
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: sentieon_docker
    memory: memory
    cpu: n_threads
    disks: job_disk_size  # Disk usage should be 2x the size of the input BAMs
  }
  output {
    # Alignment files
    File aligned_reads = "~{sample_name}_aligned.~{alignment_suffix}"
    File aligned_index = "~{sample_name}_aligned.~{alignment_suffix}.~{alignment_index}"

    # QC output metrics
    File dedup_metrics = "~{sample_name}_dedup_metrics.txt"
    File mq_metrics = "~{sample_name}_mq_metrics.txt"
    File qd_metrics = "~{sample_name}_qd_metrics.txt"
    File gc_summary = "~{sample_name}_gc_summary.txt"
    File gc_metrics = "~{sample_name}_gc_metrics.txt"
    File as_metrics = "~{sample_name}_aln_metrics.txt"
    File is_metrics = "~{sample_name}_is_metrics.txt"

    # QC output plots
    File mq_plot = "~{sample_name}_mq-report.pdf"
    File qd_plot = "~{sample_name}_qd-report.pdf"
    File gc_plot = "~{sample_name}_gc-report.pdf"
    File is_plot = "~{sample_name}_is-report.pdf"
  }
}

task ReadWriter {
  input {
    # Input files
    Array[File] bams
    Array[File] bams_idx

    # Reference genome files
    File ref_fasta
    File ref_fai

    # Workflow arguments
    String sample_name = "sample"
    Boolean output_cram = true
    String rw_xargs = "--cram_write_options version=3.0,compressor=rans"

    # Sentieon license configuration
    File? sentieon_license_file
    String? sentieon_license_server

    # Execution
    String n_threads = "16"
    String memory = "10 GB"
    Int preemptible_tries = 2
    Int disk_size = 50
    String sentieon_docker

    # Expressions
    String alignment_suffix = if output_cram then "cram" else "bam"
    String alignment_index = if output_cram then "crai" else "bai"

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

    # Ensure each BAM is adjecent to it's index
    bam_files=("~{sep='" "' bams}")
    bai_files=("~{sep='" "' bams_idx}")
    bam_input_str=""
    for i in $(seq 1 ${#bam_files[@]}); do
      i=$((i - 1))
      ln -s "${bam_files[$i]}" ./input_bam_"$i".bam
      ln -s "${bai_files[$i]}" ./input_bam_"$i".bam.bai
      bam_input_str="$bam_input_str -i ./input_bam_$i.bam"
    done

    # Merge the aligned BAM files
    sentieon driver $bam_input_str -r ~{ref_fasta} --algo ReadWriter \
      ~{rw_xargs} \
      "~{sample_name}_aligned.~{true="cram" false="bam" output_cram}"
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: sentieon_docker
    memory: memory
    cpu: n_threads
    disks: job_disk_size  # Disk usage should be ~3x the size of the input fastq
  }
  output {
    # Alignment files
    File aligned_reads = "~{sample_name}_aligned.~{alignment_suffix}"
    File aligned_index = "~{sample_name}_aligned.~{alignment_suffix}.~{alignment_index}"
  }
}

task QualCal {
  input {
    # Input files
    File aligned_reads
    File aligned_index

    # Reference genome files
    File ref_fasta
    File ref_fai

    # Known sites VCFs
    Array[File] bqsr_vcfs
    Array[File] bqsr_vcf_tbis

    # Interval files - Both are optional but recommended
    File? bqsr_intervals

    # Workflow arguments
    String sample_name = "sample"
    String qcal_xargs = ""

    # Sentieon license configuration
    File? sentieon_license_file
    String? sentieon_license_server

    # Execution
    String n_threads = "16"
    String memory = "10 GB"
    Int preemptible_tries = 2
    Int disk_size = 50
    String sentieon_docker

    # Expressions
    Boolean has_bqsr_vcfs = length(bqsr_vcfs) > 0

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

    # BQSR
    bqsr_intervals="~{default='' bqsr_intervals}"
    has_bqsr_vcfs=~{true="true" false="" has_bqsr_vcfs}
    sentieon driver -r ~{ref_fasta} -i "~{aligned_reads}" \
      ${bqsr_intervals:+--interval "$bqsr_intervals"} --algo QualCal \
      ${has_bqsr_vcfs:+-k }~{sep=" -k " bqsr_vcfs} \
      ~{qcal_xargs} "~{sample_name}_recal.table"
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: sentieon_docker
    memory: memory
    cpu: n_threads
    disks: job_disk_size  # Disk usage should be ~3x the size of the input fastq
  }
  output {
    File bqsr_table = "~{sample_name}_recal.table"
  }
}
