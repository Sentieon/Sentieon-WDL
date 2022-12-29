version 1.0

# Germline variant calling with Sentieon DNAseq or DNAscope

import "../../tasks/alignment.wdl" as Alignment
import "../../tasks/preprocessing.wdl" as Preprocessing
import "../../tasks/variant_calling.wdl" as Calling

workflow sentieon_germline {
  input {
    # Input BAM files
    Array[File] bam
    Array[File] bam_index
    File? bqsr_table

    # Reference genome files
    File ref_fasta
    File ref_fai
    File? ref_alt
    File? ref_bwt
    File? ref_sa
    File? ref_amb
    File? ref_ann
    File? ref_pac

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
    ## Realignment with BWA
    Boolean realign_input = false # Realign the input BAM file with BWA
    ## Dedup and merged BAM output
    Boolean run_dedup_and_qc = false # Mark duplicates and perform QC in addition to alignment
    Boolean output_cram = true
    ## Run an optional ReadWriter step
    Boolean run_readwriter = false
    ## BQSR
    Boolean run_bqsr = false # Calibrate a BQSR model and use it during variant calling
    ## Variant calling
    Boolean run_calling = true # Run variant calling
    String calling_algo = "Haplotyper"
    Boolean output_gvcf = false

    # Optional process arguments
    ## Realignment with BWA
    String fastq_xargs = ""
    String bwa_xargs = ""
    String bwa_karg = "10000000"
    String sort_xargs = "--bam_compression 1"
    ## Dedup and merged BAM output
    String lc_xargs = ""
    String dedup_xargs = "--cram_write_options version=3.0,compressor=rans"
    String rw_xargs = "--cram_write_options version=3.0,compressor=rans"
    String alnstat_adapter_seq = ""  # The adapter sequence for the AlignmentStat algo
    ## BQSR
    String qcal_xargs = ""
    ## Variant calling
    String calling_driver_xargs = ""
    String calling_algo_xargs = "--pcr_indel_model none"
    File? dnascope_model

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

  # Optional BAM realignment
  if (realign_input) {
    Array[Int] bam_range = range(length(bam))

    # bwa index files are required for realignment
    File bwa_ref_bwt = select_first([ref_bwt])
    File bwa_ref_sa = select_first([ref_sa])
    File bwa_ref_amb = select_first([ref_amb])
    File bwa_ref_ann = select_first([ref_ann])
    File bwa_ref_pac = select_first([ref_pac])

    scatter(i in bam_range) {
      call Alignment.RealignBam {
        input:
          bam = bam[i],

          sample_name = sample_name,
          bwa_xargs = bwa_xargs,
          bwa_karg = bwa_karg,
          sort_xargs = sort_xargs,

          ref_fasta = ref_fasta,
          ref_fai = ref_fai,
          ref_alt = ref_alt,
          ref_bwt = bwa_ref_bwt,
          ref_sa = bwa_ref_sa,
          ref_amb = bwa_ref_amb,
          ref_ann = bwa_ref_ann,
          ref_pac = bwa_ref_pac,

          sentieon_license_file = sentieon_license_file,
          sentieon_license_server = sentieon_license_server,

          n_threads = bwa_n_threads,
          memory = bwa_memory,
          disk_size = bwa_disk_size,
          preemptible_tries = preemptible_tries,
          sentieon_docker = sentieon_docker,
      }
    }
  }
  Array[File] input_aln = select_first([RealignBam.aligned_reads, bam])
  Array[File] input_aln_idx = select_first([RealignBam.aligned_index, bam_index])

  if (run_dedup_and_qc) {
    call Preprocessing.DedupAndQc {
      input:
        bams = input_aln,
        bams_idx = input_aln_idx,
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
  if (!run_dedup_and_qc && run_readwriter) {
    call Preprocessing.ReadWriter {
      input:
        bams = input_aln,
        bams_idx = input_aln_idx,
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

  if (run_dedup_and_qc || run_readwriter) {
     File merged_aln = select_first([DedupAndQc.aligned_reads, ReadWriter.aligned_reads])
     File merged_aln_idx = select_first([DedupAndQc.aligned_index, ReadWriter.aligned_index])
     Array[File] merged_aln_files = [merged_aln]
     Array[File] merged_aln_idxs = [merged_aln_idx]
  }
  Array[File] calling_alns = select_first([merged_aln_files, bam])
  Array[File] calling_idxs = select_first([merged_aln_idxs, bam_index])

  if (run_bqsr) {
    call Preprocessing.QualCal {
      input:
        aligned_reads = calling_alns,
        aligned_index = calling_idxs,
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
  if (run_bqsr || defined(bqsr_table)) {
    File calling_bqsr_table = select_first([bqsr_table, QualCal.bqsr_table])
  }

  if (run_calling) {
    call Calling.GermlineCalling {
      input:
        aligned_reads = calling_alns,
        aligned_index = calling_idxs,
        bqsr_table = calling_bqsr_table,
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
    File? aligned_reads = merged_aln
    File? aligned_index = merged_aln_idx

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

