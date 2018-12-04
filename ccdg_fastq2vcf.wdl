workflow sentieon_ccdg_fastq_vcf {
  # Inputs
  Array[String] fastq_r1s
  Array[String] fastq_r2s
  File? fof_fastq_r1s
  Array[File] file_fastq_r1s = if defined(fof_fastq_r1s) then read_lines(fof_fastq_r1s) else []
  File? fof_fastq_r2s
  Array[File] file_fastq_r2s = if defined(fof_fastq_r2s) then read_lines(fof_fastq_r2s) else []
  Array[String] read_groups
  String sample_name

  # Known sites
  File? dbsnp_vcf = "gs://sentieon-test/pipeline_test/reference/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
  File? dbsnp_index = "gs://sentieon-test/pipeline_test/reference/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
  Array[File] bqsr_vcfs = ["gs://sentieon-test/pipeline_test/reference/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz", "gs://sentieon-test/pipeline_test/reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz", "gs://sentieon-test/pipeline_test/reference/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"]
  Array[File] bqsr_tbis = ["gs://sentieon-test/pipeline_test/reference/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi", "gs://sentieon-test/pipeline_test/reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi", "gs://sentieon-test/pipeline_test/reference/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"]

  # Reference files
  File ref_fasta = "gs://sentieon-test/pipeline_test/reference/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
  File ref_fasta_index = "gs://sentieon-test/pipeline_test/reference/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
  File ref_dict = "gs://sentieon-test/pipeline_test/reference/hg38/GRCh38_full_analysis_set_plus_decoy_hla.dict"
  File? ref_alt = "gs://sentieon-test/pipeline_test/reference/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa.alt"
  File ref_bwt = "gs://sentieon-test/pipeline_test/reference/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt"
  File ref_sa = "gs://sentieon-test/pipeline_test/reference/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa.sa"
  File ref_amb = "gs://sentieon-test/pipeline_test/reference/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa.amb"
  File ref_ann = "gs://sentieon-test/pipeline_test/reference/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa.ann"
  File ref_pac = "gs://sentieon-test/pipeline_test/reference/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa.pac"

  # Workflow configurations
  ## Optional workflow stages
  String? output_bucket
  Boolean output_align_file = false # If false, output cram will be empty
  Boolean upload_gvcf = true
  Boolean output_gvcf = false
  Boolean output_vcf = false
  ## BQSR intervals
  String? bqsr_intervals = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
  ## Readwriter read_filter args, starting after the QualCalFilter table
  String readwriter_readfilter_args = ",prior=-1.0,indel=false,levels=10/20/30,min_qual=6"
  ## Variant calling algorithm
  String calling_algo = "Haplotyper"
  String? calling_intervals = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
  ## Extra driver parameters
  String lc_driver_args = "--traverse_param=200000/10000"
  String dedup_driver_args = "--traverse_param=200000/10000"
  String readwriter_driver_args = ""
  String bqsr_driver_args = ""
  String calling_driver_args = ""
  String genotyping_driver_args = ""
  ## Extra algo parameters
  String bwa_args = "-Y"
  String bwa_chunk_size = "10000000"
  String sort_args = "--block_size 512M --bam_compression 1"
  String lc_args = ""
  String dedup_args = ""
  String readwriter_args = "--cram_write_options version=3.0,compressor=rans"
  String bqsr_args = ""
  String calling_args = ""
  String genotyping_args = ""
  ## Alignment file formats
  Boolean alignment_cram = false
  Boolean dedup_cram = false
  Boolean readwriter_cram = true

  # Sentieon License configuration
  File? sentieon_license_file
  String? sentieon_license_server = "gcp.sentieon.com:9003"
  Boolean use_instance_metadata = true
  String? sentieon_auth_mech = "GOOGLE"
  String? sentieon_license_key

  # Execution configuration
  String disk = "local-disk 256 LOCAL"
  String threads = "64"
  String memory = "55 GB"
  Int preemptible_tries = 3
  String sentieon_version = "201808.01"
  String docker = "sentieon/sentieon-google-cloud:${sentieon_version}"
  String sentieon_release_dir = "/opt/sentieon/sentieon-genomics-${sentieon_version}"
  
  # One big shell script
  call SentieonFastqToVcf {
    input:
      # Inputs
      fastq_r1s = fastq_r1s,
      fastq_r2s = fastq_r2s,
      file_fastq_r1s = file_fastq_r1s,
      file_fastq_r2s = file_fastq_r2s,
      read_groups = read_groups,
      sample_name = sample_name,
      # Known sites
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_index = dbsnp_index,
      bqsr_vcfs = bqsr_vcfs,
      bqsr_tbis = bqsr_tbis,
      # Reference files
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_alt = ref_alt,
      ref_bwt = ref_bwt,
      ref_sa = ref_sa,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_pac = ref_pac,
      # Workflow configurations
      ## Optional workflow stages
      output_bucket = output_bucket,
      output_align_file = output_align_file,
      upload_gvcf = upload_gvcf,
      output_gvcf = output_gvcf,
      output_vcf = output_vcf,
      ## BQSR intervals
      bqsr_intervals = bqsr_intervals,
      ## BQSR read_filter args, starting with QualCalFilter
      readwriter_readfilter_args = readwriter_readfilter_args,
      ## Variant calling parameters
      calling_algo = calling_algo,
      calling_intervals = calling_intervals,
      ## Extra driver parameters
      lc_driver_args = lc_driver_args,
      dedup_driver_args = dedup_driver_args,
      readwriter_driver_args = readwriter_driver_args,
      bqsr_driver_args = bqsr_driver_args,
      calling_driver_args = calling_driver_args,
      genotyping_driver_args = genotyping_driver_args,
      ## Extra algo parameters
      bwa_args = bwa_args,
      bwa_chunk_size = bwa_chunk_size,
      sort_args = sort_args,
      lc_args = lc_args,
      dedup_args = dedup_args,
      readwriter_args = readwriter_args,
      bqsr_args = bqsr_args,
      calling_args = calling_args,
      genotyping_args = genotyping_args,
      ## Alignment file formats
      alignment_cram = alignment_cram,
      dedup_cram = dedup_cram,
      readwriter_cram = readwriter_cram,
      # Sentieon License configuration
      sentieon_license_server = sentieon_license_server,
      sentieon_license_file = sentieon_license_file,
      use_instance_metadata = use_instance_metadata,
      sentieon_auth_mech = sentieon_auth_mech,
      sentieon_license_key = sentieon_license_key,
      # Execution configuration
      disk = disk,
      threads = threads,
      memory = memory,
      preemptible_tries = preemptible_tries,
      docker = docker,
      sentieon_release_dir = sentieon_release_dir
  }
  output {
    File preprocessed_alignment = SentieonFastqToVcf.alignment
    File preprocessed_alignment_index = SentieonFastqToVcf.alignment_index
    File gvcf = SentieonFastqToVcf.gvcf
    File gvcf_index = SentieonFastqToVcf.gvcf_index
    File vcf = SentieonFastqToVcf.vcf
    File vcf_index = SentieonFastqToVcf.vcf_index
  }
}

task SentieonFastqToVcf {
  # Inputs
  Array[String] fastq_r1s
  Array[String] fastq_r2s
  Array[File] file_fastq_r1s
  Array[File] file_fastq_r2s
  Array[String] read_groups
  String sample_name

  # Known sites
  File? dbsnp_vcf
  File? dbsnp_index
  Array[File] bqsr_vcfs
  Array[File] bqsr_tbis

  # Reference files
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File? ref_alt
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac

  # Workflow configurations
  ## Optional workflow stages
  String? output_bucket
  Boolean output_align_file
  Boolean upload_gvcf
  Boolean output_gvcf
  Boolean output_vcf
  ## BQSR intervals
  String? bqsr_intervals
  ## Readwriter read_filter args, starting after the QualCalFilter table
  String readwriter_readfilter_args
  ## Variant calling algorithm
  String calling_algo
  String calling_intervals
  ## Extra driver parameters
  String lc_driver_args
  String dedup_driver_args
  String readwriter_driver_args
  String bqsr_driver_args
  String calling_driver_args
  String genotyping_driver_args
  ## Extra algo parameters
  String bwa_args
  String bwa_chunk_size
  String sort_args
  String lc_args
  String dedup_args
  String readwriter_args
  String bqsr_args
  String calling_args
  String genotyping_args
  ## Alignment file formats
  Boolean alignment_cram
  Boolean dedup_cram
  Boolean readwriter_cram

  # Sentieon License configuration
  File? sentieon_license_file
  String? sentieon_license_server
  Boolean use_instance_metadata
  String? sentieon_auth_mech
  String? sentieon_license_key

  # Execution configuration
  String disk
  String threads
  String memory
  Int preemptible_tries
  String docker
  String sentieon_release_dir
  
  # Some preprocessing
  Boolean make_gvcf = upload_gvcf || output_gvcf
  Boolean call_variants = make_gvcf || output_vcf
  Boolean emit_vcf = output_vcf && (!make_gvcf)
  Boolean run_genotyper = make_gvcf && output_vcf
  String readwriter_suffix = if readwriter_cram then "cram" else "bam"
  String readwriter_index_suffix = if readwriter_cram then "cram.crai" else "bam.bai"
  String dollar = "$"
  command <<<
    set -exo pipefail

    # Check that the configuration is valid.
    # Supported variant callers are Genotyper, Haplotyper and DNAscope
    if [[ "${calling_algo}" != "Genotyper" && "${calling_algo}" != "Haplotyper" && "${calling_algo}" != "DNAscope" ]]; then
      echo "${calling_algo} is not a supported variant caller. Please set calling_algo to 'Genotyper', 'Haplotyper' or 'DNAscope'" >&2
      exit 1
    fi
    # Must supply an output bucket to upload the gVCF
    if [[ -n '${true="y" false="" upload_gvcf}' && -z '${default="" output_bucket}' ]]; then
      echo "Must supply an output bucket to upload a gVCF" >&2
      exit 1
    fi
    # Number of readgroups must match the number of fastq files
    first_fastq=(${sep=" " fastq_r1s})
    second_fastq=(${sep=" " fastq_r2s})
    first_fastq+=(${sep=" " file_fastq_r1s})
    second_fastq+=(${sep=" " file_fastq_r2s})
    read_groups=('${sep="' '" read_groups}')
    if [[ ${dollar}{#first_fastq[@]} -ne ${dollar}{#read_groups[@]} ]]; then
      echo "The number of fastq files for r1 does not equal the number for r2"
      exit 1
    fi
    if [[ ${dollar}{#first_fastq[@]} -ne ${dollar}{#read_groups[@]} ]]; then
      echo "The number of fastq pairs does not equal the number of supplied readgroups"
      exit 1
    fi

    wait_list=()
    # License server setup
    license_file=${default="" sentieon_license_file}
    if [[ -n "$license_file" ]]; then
      # Using a license file
      export SENTIEON_LICENSE=${default="" sentieon_license_file}
    elif [[ -n '${true="yes" false="" use_instance_metadata}' ]]; then
      python /opt/sentieon/gen_credentials.py ~/credentials.json ${default="''" sentieon_license_key} &
      sleep 5
      export SENTIEON_LICENSE=${default="" sentieon_license_server}
      export SENTIEON_AUTH_MECH=${default="" sentieon_auth_mech}
      export SENTIEON_AUTH_DATA=~/credentials.json
      read -r SENTIEON_JOB_TAG < ~/credentials.json.project
      export SENTIEON_JOB_TAG
    else
      export SENTIEON_LICENSE=${default="" sentieon_license_server}
      export SENTIEON_AUTH_MECH=${default="" sentieon_auth_mech}
    fi

    # Optimizations
    mem_kb=$(cat /proc/meminfo | grep "MemTotal" | awk '{print $2}')
    export bwt_max_mem="$((mem_kb / 1024 / 1024 - 2))g"
    export MALLOC_CONF=lg_dirty_mult:-1

    # Alignment with BWA
    aligned_data=()
    for i in $(seq 1 ${dollar}{#first_fastq[@]}); do
      i=$((i - 1))
      LD_PRELOAD=${sentieon_release_dir}/lib/libjemalloc.so sentieon bwa mem -t ${threads} ${bwa_args} -K ${bwa_chunk_size} -R "${dollar}{read_groups[$i]}" ${ref_fasta} <(gsutil cp ${dollar}{first_fastq[$i]} -) <(gsutil cp ${dollar}{second_fastq[$i]} -) | \
        samblaster --addMateTags -a | \
        LD_PRELOAD=${sentieon_release_dir}/lib/libjemalloc.so sentieon util sort ${sort_args} -t ${threads} -i - --sam2bam -o ${sample_name}_sorted_${dollar}{i}.${true="cram" false="bam" alignment_cram}
      aligned_data+=(${sample_name}_sorted_${dollar}{i}.${true="cram" false="bam" alignment_cram})
    done

    export LD_PRELOAD=${sentieon_release_dir}/lib/libjemalloc.so

    # Dedup
    bam_input=""
    for f in ${dollar}{aligned_data[@]}; do
      bam_input="$bam_input -i $f"
    done
    sentieon driver ${lc_driver_args} -t ${threads} $bam_input --algo LocusCollector ${lc_args} ${sample_name}_score.txt
    sentieon driver ${dedup_driver_args} -t ${threads} $bam_input --algo Dedup ${dedup_args} --score_info ${sample_name}_score.txt --output_dup_read_name ${sample_name}_dup_qname.txt
    sentieon driver ${dedup_driver_args} -t ${threads} $bam_input --algo Dedup ${dedup_args} --dup_read_name ${sample_name}_dup_qname.txt ${sample_name}_deduped.${true="cram" false="bam" dedup_cram}
    for f in ${dollar}{aligned_data[@]}; do
      rm "$f" &
    done

    # BQSR
    # Group the known sites files with their indicies
    sites_vcfs=(${sep=" " bqsr_vcfs})
    sites_tbis=(${sep=" " bqsr_tbis})
    sites_str=""
    mkdir -p inputs
    for i in $(seq 1 ${dollar}{#sites_vcfs[@]}); do
      i=$((i - 1))
      ln -s ${dollar}{sites_vcfs[$i]} inputs/bqsr_sites_${dollar}{i}.vcf.gz
      ln -s ${dollar}{sites_tbis[$i]} inputs/bqsr_sites_${dollar}{i}.vcf.gz.tbi
      sites_str+=" -k inputs/bqsr_sites_${dollar}{i}.vcf.gz "
    done
    sentieon driver ${"--interval " + bqsr_intervals} -r ${ref_fasta} -t ${threads} -i ${sample_name}_deduped.${true="cram" false="bam" dedup_cram} ${bqsr_driver_args} --algo QualCal $sites_str ${bqsr_args} ${sample_name}_recal.table
    
    # ReadWriter
    sentieon driver -r ${ref_fasta} -t ${threads} -i ${sample_name}_deduped.${true="cram" false="bam" dedup_cram} --read_filter QualCalFilter,table=${sample_name}_recal.table${readwriter_readfilter_args} ${readwriter_driver_args} --algo ReadWriter ${readwriter_args} ${sample_name}_recal.${readwriter_suffix}
    rm ${sample_name}_deduped.${true="cram" false="bam" dedup_cram} &

    # Optionally uploaded the aligned data 
    if [[ -n '${default="" output_bucket}' ]]; then
      gsutil cp ${sample_name}_recal.${readwriter_suffix} "${default='' output_bucket}" &
      wait_list+=($!)
      gsutil cp ${sample_name}_recal.${readwriter_index_suffix} "${default='' output_bucket}" &
      wait_list+=($!)
    fi

    # Ensure all output files are present so Cromwell does not error if they are streamed
    touch ${sample_name}_${calling_algo}.g.vcf.gz ${sample_name}_${calling_algo}.g.vcf.gz.tbi ${sample_name}_${calling_algo}.vcf.gz ${sample_name}_${calling_algo}.vcf.gz.tbi
    if [[ -n '${true="y" false="" call_variants}' ]]; then
      # Call variants
      sentieon driver -r ${ref_fasta} -t ${threads} -i ${sample_name}_recal.${readwriter_suffix} ${calling_driver_args} ${"--interval " + calling_intervals} --algo ${calling_algo} ${"-d " + dbsnp_vcf} ${calling_args} ${true="--emit_mode GVCF" false="" make_gvcf} ${sample_name}_${calling_algo}${true=".g" false="" make_gvcf}.vcf.gz
      # Optionally upload the gVCF
      if [[ -n '${true="y" false="" upload_gvcf}' ]]; then
        gsutil cp ${sample_name}_${calling_algo}${true=".g" false="" make_gvcf}.vcf.gz "${default='' output_bucket}" &
        wait_list+=($!)
        gsutil cp ${sample_name}_${calling_algo}${true=".g" false="" make_gvcf}.vcf.gz.tbi "${default='' output_bucket}" &
        wait_list+=($!)
      fi
      if [[ -n '${true="y" false="" run_genotyper}' ]]; then
        # Genotype the GVCF
        sentieon driver -r ${ref_fasta} -t ${threads} ${genotyping_driver_args} --algo GVCFtyper ${genotyping_args} ${sample_name}_${calling_algo}.vcf.gz ${sample_name}_${calling_algo}.g.vcf.gz
      fi
    fi

    # Wait
    for pid in ${dollar}{wait_list[@]}; do
      wait $pid
    done

    # Finalize output files
    if [[ -n '${true="" false="yes" output_align_file}' ]]; then
      rm ${sample_name}_recal.${readwriter_suffix} ${sample_name}_recal.${readwriter_index_suffix}
      touch ${sample_name}_recal.${readwriter_suffix} ${sample_name}_recal.${readwriter_index_suffix}
    fi
    if [[ -n '${true="" false="yes" output_gvcf}' ]]; then
      rm ${sample_name}_${calling_algo}.g.vcf.gz ${sample_name}_${calling_algo}.g.vcf.gz.tbi
      touch ${sample_name}_${calling_algo}.g.vcf.gz ${sample_name}_${calling_algo}.g.vcf.gz.tbi
    fi
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: docker
    memory: memory
    cpu: threads
    disks: disk
  }
  output {
    File alignment = "${sample_name}_recal.${readwriter_suffix}"
    File alignment_index = "${sample_name}_recal.${readwriter_index_suffix}"
    File gvcf = "${sample_name}_${calling_algo}.g.vcf.gz"
    File gvcf_index = "${sample_name}_${calling_algo}.g.vcf.gz.tbi"
    File vcf = "${sample_name}_${calling_algo}.vcf.gz"
    File vcf_index = "${sample_name}_${calling_algo}.vcf.gz.tbi"
  }
}
