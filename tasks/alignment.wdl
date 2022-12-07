version 1.0

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
