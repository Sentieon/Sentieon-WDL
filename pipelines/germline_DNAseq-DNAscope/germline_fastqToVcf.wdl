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
    String n_threads = "32"
    String memory = "64 GB"
    Int preemptible_tries = 3
    Int? disk_size
    String sentieon_docker
  }
  call SentieonGermline {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      read_groups = read_groups,

      ref_fasta = ref_fasta,
      ref_fai = ref_fai,
      ref_alt = ref_alt,
      ref_bwt = ref_bwt,
      ref_sa = ref_sa,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_pac = ref_pac,

      bqsr_vcfs = bqsr_vcfs,
      bqsr_vcf_tbis = bqsr_vcf_tbis,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_tbi = dbsnp_vcf_tbi,

      bqsr_intervals = bqsr_intervals,
      calling_intervals = calling_intervals,

      sample_name = sample_name,
      output_cram = output_cram,
      run_dedup_and_qc = run_dedup_and_qc,
      run_bqsr = run_bqsr,
      run_calling = run_calling,
      calling_algo = calling_algo,
      output_gvcf = output_gvcf,

      bwa_xargs = bwa_xargs,
      bwa_karg = bwa_karg,
      sort_xargs = sort_xargs,
      lc_xargs = lc_xargs,
      dedup_xargs = dedup_xargs,
      rw_xargs = rw_xargs,
      qcal_xargs = qcal_xargs,
      calling_driver_xargs = calling_driver_xargs,
      calling_algo_xargs = calling_algo_xargs,
      dnascope_model = dnascope_model,
      alnstat_adapter_seq = alnstat_adapter_seq,

      sentieon_license_file = sentieon_license_file,
      sentieon_license_server = sentieon_license_server,

      n_threads = n_threads,
      memory = memory,
      preemptible_tries = preemptible_tries,
      disk_size = disk_size,
      sentieon_docker = sentieon_docker,
  }
  output {
    # Core alignment files
    File aligned_reads = SentieonGermline.aligned_reads
    File aligned_index = SentieonGermline.aligned_index

    # DNAseq outputs
    File? bqsr_table = SentieonGermline.bqsr_table
    File? calls_vcf = SentieonGermline.calls_vcf
    File? calls_vcf_tbi = SentieonGermline.calls_vcf_tbi

    # QC output metrics
    File? dedup_metrics = SentieonGermline.dedup_metrics
    File? mq_metrics = SentieonGermline.mq_metrics
    File? qd_metrics = SentieonGermline.qd_metrics
    File? gc_summary = SentieonGermline.gc_summary
    File? gc_metrics = SentieonGermline.gc_metrics
    File? as_metrics = SentieonGermline.as_metrics
    File? is_metrics = SentieonGermline.is_metrics

    # QC output plots
    File? mq_plot = SentieonGermline.mq_plot
    File? qd_plot = SentieonGermline.qd_plot
    File? gc_plot = SentieonGermline.gc_plot
    File? is_plot = SentieonGermline.is_plot
  }
}

task SentieonGermline {
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
    String n_threads = "32"
    String memory = "64 GB"
    Int preemptible_tries = 3
    Int? disk_size
    String sentieon_docker

    # Expressions
    Boolean has_bqsr_vcfs = length(bqsr_vcfs) > 0
    String alignment_suffix = if output_cram then "cram" else "bam"
    String alignment_index = if output_cram then "crai" else "bai"
    String vcf_suffix = if output_gvcf then ".g.vcf.gz" else ".vcf.gz"
    String vcf_index = if output_gvcf then ".g.vcf.gz.tbi" else ".vcf.gz.tbi"

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

    r1_fastq=("~{sep='" "' r1_fastq}")
    r2_fastq=("~{sep='" "' r2_fastq}")
    read_groups=("~{sep='" "' read_groups}")

    # Sanity check the input
    if [[ ${#r1_fastq[@]} -ne ${#read_groups[@]} ]]; then
      echo "The number of readgroups does not match the number of input fastq files"
      exit 1
    fi
    if [[ ${#r1_fastq[@]} -ne ${#r2_fastq[@]} ]]; then
      echo "The number of r1 fastq does not equal the number of r2 fastq"
      exit 1
    fi

    # Sanity check the calling algo
    calling_algo=~{calling_algo}
    if ! [[ "$calling_algo" =~ ^(Haplotyper|DNAscope|Genotyper)$ ]]; then
      echo "Supported calling algos are Haplotyper, DNAscope, and Genotyper"
      exit 1
    fi

    # Sanity check the sample name
    sample_name="~{sample_name}"
    pattern="	| "
    if [[ "$sample_name" =~ "$pattern" ]]; then
      echo "Sample name should not contain whitespace"
      exit 1
    fi

    # Alignment with BWA
    alignment_output=()
    bam_str=""
    for i in $(seq 1 ${#r1_fastq[@]}); do
      i=$((i - 1))
      r1="${r1_fastq[$i]}"
      r2="${r2_fastq[$i]}"
      rg="${read_groups[$i]}"

      for j in $(seq 1 "$numa_nodes"); do
        j=$((j - 1))
        cpulist="${numa_cpulist[$j]}"

        # Alignment command
        perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)';
        taskset -c "$cpulist" sentieon bwa mem -R "$rg" \
          ~{bwa_xargs} -K ~{bwa_karg} -t $nt -p "~{ref_fasta}" \
          <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
            sentieon fqidx extract -F "$j"/"$numa_nodes" -K ~{bwa_karg} \
            <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
              igzip -dc "$r1") \
            <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
              igzip -dc "$r2")) | \
          taskset -c "$cpulist" sentieon util sort -t $nt --sam2bam \
          -o "~{sample_name}_sorted_${i}_${j}.bam" -i - ~{sort_xargs} &
        alignment_output+=("~{sample_name}_sorted_${i}_${j}.bam")
        bam_str="$bam_str -i ~{sample_name}_sorted_${i}_${j}.bam"
      done
      wait

      (rm "$r1" "$r2" || (exit 0)) &  # Try removing the input files to save disk space
    done

    # Dedup
    run_dedup=~{true="true" false="false" run_dedup_and_qc}
    output_aln="~{sample_name}_aligned.~{true="cram" false="bam" output_cram}"
    if [[ "$run_dedup" == "true" ]]; then
      # Perform duplicate marking and QC
      sentieon driver $bam_str -r ~{ref_fasta} \
        --algo LocusCollector ~{lc_xargs} "~{sample_name}_score.txt.gz" \
        --algo MeanQualityByCycle "~{sample_name}_mq_metrics.txt" \
        --algo QualDistribution "~{sample_name}_qd_metrics.txt" \
        --algo GCBias --summary "~{sample_name}_gc_summary.txt" "~{sample_name}_gc_metrics.txt" \
        --algo AlignmentStat --adapter_seq '~{alnstat_adapter_seq}' "~{sample_name}_aln_metrics.txt" \
        --algo InsertSizeMetricAlgo "~{sample_name}_is_metrics.txt"

      sentieon driver $bam_str -r ~{ref_fasta} --algo Dedup \
        ~{dedup_xargs} --score_info "~{sample_name}_score.txt.gz" \
        --metrics "~{sample_name}_dedup_metrics.txt" \
        "$output_aln"

      # Plot the metrics output
      sentieon plot GCBias -o "~{sample_name}_gc-report.pdf" "~{sample_name}_gc_metrics.txt" &
      sentieon plot QualDistribution -o "~{sample_name}_qd-report.pdf" "~{sample_name}_qd_metrics.txt" &
      sentieon plot MeanQualityByCycle -o "~{sample_name}_mq-report.pdf" "~{sample_name}_mq_metrics.txt" &
      sentieon plot InsertSizeMetricAlgo -o "~{sample_name}_is-report.pdf" "~{sample_name}_is_metrics.txt" &
    else
      # Merge the aligned BAM files
      sentieon driver $bam_str -r ~{ref_fasta} --algo ReadWriter \
        ~{rw_xargs} "$output_aln"
    fi
    rm ${alignment_output[@]}  # Remove intermediate files to save space

    # BQSR
    run_bqsr=~{true="true" false="" run_bqsr}
    bqsr_out=${run_bqsr:+"~{sample_name}_recal.table"}
    bqsr_intervals="~{default='' bqsr_intervals}"
    has_bqsr_vcfs=~{true="true" false="" has_bqsr_vcfs}
    if [[ "$run_bqsr" == "true" ]]; then
      sentieon driver -r ~{ref_fasta} -i "$output_aln" \
        ${bqsr_intervals:+--interval "$bqsr_intervals"} --algo QualCal \
        ${has_bqsr_vcfs:+-k }~{sep=" -k " bqsr_vcfs} \
        ~{qcal_xargs} "$bqsr_out"
    fi

    # Variant calling
    run_calling=~{true="true" false="" run_calling}
    calling_intervals="~{default='' calling_intervals}"
    dbsnp_vcf="~{default='' dbsnp_vcf}"
    dnascope_model=""
    if [[ "~{calling_algo}" == "DNAscope" ]]; then
      dnascope_model="~{default="" dnascope_model}"
    fi
    output_vcf="~{sample_name}.~{true="g.vcf.gz" false="vcf.gz" output_gvcf}"
    if [[ "$run_calling" == "true" ]]; then
      sentieon driver -r ~{ref_fasta} -i "$output_aln" \
        ${bqsr_out:+-q "$bqsr_out"} \
        ${calling_intervals:+--interval "$calling_intervals"} \
        ~{calling_driver_xargs} \
        --algo ~{calling_algo} \
        ~{true="--emit_mode gvcf" false="" output_gvcf} \
        ${dbsnp_vcf:+-d "$dbsnp_vcf"} \
        ${dnascope_model:+--model "$dnascope_model"} \
        ~{calling_algo_xargs} \
        "$output_vcf"

      # Run DNAModelApply when DNAscope is run with a model file
      if [[ -n "$dnascope_model" ]]; then
        tmp_vcf="~{sample_name}.tmp.~{true="g.vcf.gz" false="vcf.gz" output_gvcf}"
        mv "$output_vcf" "$tmp_vcf"
        mv "$output_vcf".tbi "$tmp_vcf".tbi
        sentieon driver -r ~{ref_fasta} \
          --algo DNAModelApply --model "$dnascope_model" -v "$tmp_vcf" \
          "$output_vcf"
      fi
    fi
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
    File aligned_reads = "~{sample_name}_aligned.~{alignment_suffix}"
    File aligned_index = "~{sample_name}_aligned.~{alignment_suffix}.~{alignment_index}"

    # VCF files
    File? bqsr_table = "~{sample_name}_recal.table"
    File? calls_vcf = "~{sample_name}~{vcf_suffix}"
    File? calls_vcf_tbi = "~{sample_name}~{vcf_index}"

    # QC output metrics
    File? dedup_metrics = "~{sample_name}_dedup_metrics.txt"
    File? mq_metrics = "~{sample_name}_mq_metrics.txt"
    File? qd_metrics = "~{sample_name}_qd_metrics.txt"
    File? gc_summary = "~{sample_name}_gc_summary.txt"
    File? gc_metrics = "~{sample_name}_gc_metrics.txt"
    File? as_metrics = "~{sample_name}_aln_metrics.txt"
    File? is_metrics = "~{sample_name}_is_metrics.txt"

    # QC output plots
    File? mq_plot = "~{sample_name}_mq-report.pdf"
    File? qd_plot = "~{sample_name}_qd-report.pdf"
    File? gc_plot = "~{sample_name}_gc-report.pdf"
    File? is_plot = "~{sample_name}_is-report.pdf"
  }
}

