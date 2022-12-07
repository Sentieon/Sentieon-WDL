version 1.0

# Germline variant calling with Sentieon DNAseq or DNAscope

import "../../tasks/alignment.wdl" as Alignment
import "../../tasks/preprocessing.wdl" as Preprocessing
import "../../tasks/variant_calling.wdl" as Calling

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
    call Alignment.SentieonBWA {
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
    call Preprocessing.DedupAndQc {
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
    call Preprocessing.ReadWriter {
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
    call Preprocessing.QualCal {
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
    call Calling.GermlineCalling {
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
      call Calling.DNAModelApply {
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
