version 1.0

# Germline variant calling with Sentieon DNAseq or DNAscope

workflow sentieon_germline {
  input {
    # Input fastq files
    Array[File] r1_fastq
    Array[File] r2_fastq
    Array[String] read_groups

    # Reference genome files
    File ref_fasta
    File ref_fai
    File? ref_alt
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac

    # Known sites VCFs
    Array[File] bqsr_vcfs
    Array[File] bqsr_vcf_tbis
    File? dbsnp_vcf
    File? dbsnp_vcf_tbi

    # Interval files - Both are optional but recommended
    File? bqsr_intervals
    File? calling_intervals

    # Workflow arguments
    String sample_name = "sample"
    Boolean output_cram = true
    Boolean run_dedup_and_qc = true # Mark duplicates and perform QC in addition to alignment
    Boolean run_bqsr = true # Calibrate a BQSR model and use it during variant calling
    Boolean run_calling = true # Run variant calling
    String calling_algo = "Haplotyper"
    Boolean output_gvcf = false

    # Optional process arguments
    String bwa_xargs = ""
    String bwa_karg = "10000000"
    String sort_xargs = "--bam_compression 1"
    String lc_xargs = ""
    String dedup_xargs = "--cram_write_options version=3.0,compressor=rans"
    String rw_xargs = "--cram_write_options version=3.0,compressor=rans"
    String qcal_xargs = ""
    String calling_driver_xargs = ""
    String calling_algo_xargs = "--pcr_indel_model none"
    File? dnascope_model
    String alnstat_adapter_seq = ""  # The adapter sequence for the AlignmentStat algo

    # Sentieon license configuration
    File? sentieon_license_file
    String? sentieon_license_server

    # Execution
    String sentieon_docker
    Int preemptible_tries = 2

    # call-specific execution configuration
    String bwa_n_threads = "32"
    String bwa_memory = "60 GB"
    Int? bwa_disk_size

    String preprocess_n_threads = "16"
    String preprocess_memory = "10 GB"
    Int? preprocess_disk_size

    String calling_threads = "32"
    String calling_memory = "10 GB"
    Int? calling_disk_size
  }

  Array[Int] pair_range = range(length(r1_fastq))
  scatter(i in pair_range) {
    call SentieonBWA {
      input:
        r1_fastq = r1_fastq[i],
        r2_fastq = r2_fastq[i],
        read_group = read_groups[i],

        sample_name = sample_name,
        bwa_xargs = bwa_xargs,
        bwa_karg = bwa_karg,
        sort_xargs = sort_xargs,

        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_alt = ref_alt,
        ref_bwt = ref_bwt,
        ref_sa = ref_sa,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,

        sentieon_license_file = sentieon_license_file,
        sentieon_license_server = sentieon_license_server,

        n_threads = bwa_n_threads,
        memory = bwa_memory,
        disk_size = bwa_disk_size,
        preemptible_tries = preemptible_tries,
        sentieon_docker = sentieon_docker,
    }
  }

  Array[File] bwa_aligned_reads = flatten(SentieonBWA.aligned_reads)
  Array[File] bwa_aligned_index = flatten(SentieonBWA.aligned_index)

  if (run_dedup_and_qc) {
    call DedupAndQc {
      input:
        bams = bwa_aligned_reads,
        bams_idx = bwa_aligned_index,
        output_cram = output_cram,

        sample_name = sample_name,
        lc_xargs = lc_xargs,
        alnstat_adapter_seq = alnstat_adapter_seq,
        dedup_xargs = dedup_xargs,

        ref_fasta = ref_fasta,
        ref_fai = ref_fai,

        sentieon_license_file = sentieon_license_file,
        sentieon_license_server = sentieon_license_server,

        n_threads = preprocess_n_threads,
        memory = preprocess_memory,
        disk_size = preprocess_disk_size,
        preemptible_tries = preemptible_tries,
        sentieon_docker = sentieon_docker,
    }
  }
  if (!run_dedup_and_qc) {
    call ReadWriter {
      input:
        bams = bwa_aligned_reads,
        bams_idx = bwa_aligned_index,
        output_cram = output_cram,

        sample_name = sample_name,
        rw_xargs = rw_xargs,

        ref_fasta = ref_fasta,
        ref_fai = ref_fai,

        sentieon_license_file = sentieon_license_file,
        sentieon_license_server = sentieon_license_server,

        n_threads = preprocess_n_threads,
        memory = preprocess_memory,
        disk_size = preprocess_disk_size,
        preemptible_tries = preemptible_tries,
        sentieon_docker = sentieon_docker,
    }
  }

  File merged_aln = select_first([DedupAndQc.aligned_reads, ReadWriter.aligned_reads])
  File merged_aln_idx = select_first([DedupAndQc.aligned_index, ReadWriter.aligned_index])
  if (run_bqsr) {
    call QualCal {
      input:
        aligned_reads = merged_aln,
        aligned_index = merged_aln_idx,
        bqsr_intervals = bqsr_intervals,
        bqsr_vcfs = bqsr_vcfs,
        bqsr_vcf_tbis = bqsr_vcf_tbis,

        sample_name = sample_name,
        qcal_xargs = qcal_xargs,

        ref_fasta = ref_fasta,
        ref_fai = ref_fai,

        sentieon_license_file = sentieon_license_file,
        sentieon_license_server = sentieon_license_server,

        n_threads = preprocess_n_threads,
        memory = preprocess_memory,
        disk_size = preprocess_disk_size,
        preemptible_tries = preemptible_tries,
        sentieon_docker = sentieon_docker,
    }
  }

  if (run_calling) {
    call GermlineCalling {
      input:
        aligned_reads = merged_aln,
        aligned_index = merged_aln_idx,
        bqsr_table = QualCal.bqsr_table,
        calling_intervals = calling_intervals,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_tbi = dbsnp_vcf_tbi,

        sample_name = sample_name,
        calling_algo = calling_algo,
        output_gvcf = output_gvcf,
        calling_driver_xargs = calling_driver_xargs,
        calling_algo_xargs = calling_algo_xargs,
        dnascope_model = dnascope_model,

        ref_fasta = ref_fasta,
        ref_fai = ref_fai,

        sentieon_license_file = sentieon_license_file,
        sentieon_license_server = sentieon_license_server,

        n_threads = calling_threads,
        memory = calling_memory,
        disk_size = calling_disk_size,
        preemptible_tries = preemptible_tries,
        sentieon_docker = sentieon_docker,
    }

    if (defined(dnascope_model)) {
      # DNAModelApply is a no-op if calling_algo != 'DNAscope'
      call DNAModelApply {
        input:
          vcf = GermlineCalling.calls_vcf,
          vcf_tbi = GermlineCalling.calls_vcf_tbi,

          sample_name = sample_name,
          dnascope_model = dnascope_model,

          ref_fasta = ref_fasta,
          ref_fai = ref_fai,

          sentieon_license_file = sentieon_license_file,
          sentieon_license_server = sentieon_license_server,

          n_threads = calling_threads,
          memory = calling_memory,
          disk_size = calling_disk_size,
          preemptible_tries = preemptible_tries,
          sentieon_docker = sentieon_docker,
      }
    }
    File out_calls_vcf = select_first([DNAModelApply.calls_vcf, GermlineCalling.calls_vcf])
    File out_calls_vcf_tbi = select_first([DNAModelApply.calls_vcf_tbi, GermlineCalling.calls_vcf_tbi])
  }

  output {
    # Core alignment files
    File aligned_reads = merged_aln
    File aligned_index = merged_aln_idx

    # DNAseq outputs
    File? calls_vcf = out_calls_vcf
    File? calls_vcf_tbi = out_calls_vcf_tbi

    # QC output metrics
    File? dedup_metrics = DedupAndQc.dedup_metrics
    File? mq_metrics = DedupAndQc.mq_metrics
    File? qd_metrics = DedupAndQc.qd_metrics
    File? gc_summary = DedupAndQc.gc_summary
    File? gc_metrics = DedupAndQc.gc_metrics
    File? as_metrics = DedupAndQc.as_metrics
    File? is_metrics = DedupAndQc.is_metrics

    # QC output plots
    File? mq_plot = DedupAndQc.mq_plot
    File? qd_plot = DedupAndQc.qd_plot
    File? gc_plot = DedupAndQc.gc_plot
    File? is_plot = DedupAndQc.is_plot
  }
}

task SentieonBWA {
  input {
    # The input fastq files
    File r1_fastq
    File r2_fastq
    String read_group

    # Reference genome files
    File ref_fasta
    File ref_fai
    File? ref_alt
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac

    # Workflow arguments
    String sample_name = "sample"
    String bwa_xargs = ""
    String bwa_karg = "10000000"
    String sort_xargs = "--bam_compression 1"

    # Sentieon license configuration
    File? sentieon_license_file
    String? sentieon_license_server

    # Execution
    String n_threads = "32"
    String memory = "64 GB"
    Int preemptible_tries = 2
    Int? disk_size
    String sentieon_docker

    # Disk size
    Float r1_fq_size = size(r1_fastq, "GB")
    Float r2_fq_size = size(r2_fastq, "GB")
    Float fq_size = select_first([r1_fq_size, 1.0]) + select_first([r2_fq_size, 1.0])
    Int? fq_disk_size = ceil(select_first([fq_size, 2.0]) * 4)
    Int job_disk_size = select_first([disk_size, fq_disk_size, 8])
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

    # Get the NUMA configuration
    numa_nodes=$(lscpu | grep "NUMA node(s):" | sed 's/^NUMA node.* //')
    numa_cpulist=()
    for i in $(seq 1 "$numa_nodes"); do
        i=$((i - 1))
        numa_cpulist+=($(lscpu | grep "NUMA node$i CPU" | sed 's/^NUMA.* //'))
    done

    # Configuration
    nt=$(nproc)

    mem_kb=$(cat /proc/meminfo | grep "MemTotal" | awk '{print $2}')
    bwt_mem=$((mem_kb / 1024 / 1024 / numa_nodes - 6))
    export bwt_max_mem="$bwt_mem"G

    # Sanity check the sample name
    sample_name="~{sample_name}"
    pattern="	| "
    if [[ "$sample_name" =~ "$pattern" ]]; then
      echo "Sample name should not contain whitespace"
      exit 1
    fi

    # Alignment with BWA
    for j in $(seq 1 "$numa_nodes"); do
      j=$((j - 1))
      cpulist="${numa_cpulist[$j]}"

      # Alignment command
      perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)';
      taskset -c "$cpulist" sentieon bwa mem -R "~{read_group}" \
        ~{bwa_xargs} -K ~{bwa_karg} -t $nt -p "~{ref_fasta}" \
        <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
          sentieon fqidx extract -F "$j"/"$numa_nodes" -K ~{bwa_karg} \
          <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
            igzip -dc "~{r1_fastq}") \
          <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
            igzip -dc "~{r2_fastq}")) | \
        taskset -c "$cpulist" sentieon util sort -t $nt --sam2bam \
        -o "~{sample_name}_sorted_${i}_${j}.bam" -i - ~{sort_xargs} &
    done
    wait
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
    Array[File] aligned_reads = glob("*_sorted_*.bam")
    Array[File] aligned_index = glob("*_sorted_*.bam.bai")
  }
}

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

task GermlineCalling {
  input {
    # Input files
    File aligned_reads
    File aligned_index
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

    # Variant calling
    calling_intervals="~{default='' calling_intervals}"
    dbsnp_vcf="~{default='' dbsnp_vcf}"
    bqsr_table="~{default='' bqsr_table}"
    dnascope_model=""
    if [[ "~{calling_algo}" == "DNAscope" ]]; then
      dnascope_model="~{default="" dnascope_model}"
    fi
    sentieon driver -r ~{ref_fasta} -i "~{aligned_reads}" \
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
