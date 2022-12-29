# Sentieon DNAseq and DNAscope
WDL pipelines for germline variant calling

## Basic usage - Sentieon DNAseq
The pipeline is tested with [Cromwell](https://github.com/broadinstitute/cromwell), but should be compatible with other tools for WDL executors. Here is the basic usage for Cromwell:
```sh
java -Dconfig.file=<cromwell.config> -jar cromwell-83.jar run -i <inputs.json> <germline.wdl>
```
We provide example `inputs.json` files in the [germline_fastqToVcf.example_inputs.json](germline_fastqToVcf.example_inputs.json) and [germline_BamToVcf.example_inputs.json](germline_BamToVcf.example_inputs.json) files in this repository.

For more information on WDL, please see https://openwdl.org/.

### Using the `_split.wdl` pipeline
We provide two implementations of some pipelines in this repository, a standard implementation and a `_split.wdl` implementation. The `_split.wdl` implementation follows the more common approach of splitting the pipeline into multiple tasks, while the other WDL implements the full pipeline as a single task. Both implementations accept the same inputs and should produce identical output.

The standard pipeline is recommended for most use-cases, and will provide a faster overall turnaround time on large machines due to reduced data transfer and scheduling overhead. The `_split.wdl` may be more cost-effective in environments with frequent machine preemption.

### Output files
By default, the fastq-to-VCF pipeline will output duplicate-marked read alignments in the CRAM format, called variants in the VCF format, a Sentieon BQSR table, and sample QC metrics and plots. The BAM-to-VCF pipeline will output called variants in the VCF format, with default arguments.

As an alternative to the CRAM format, the pipeline may also output aligned reads in the BAM format by modifying the `inputs.json` file:
```json
"sentieon_germline.output_cram": false,
```

Other supported pipelines and arguments are described below.

## Alternative pipelines

The default arguments will run a Sentieon DNAseq pipeline for Fastq->VCF or BAM/CRAM->VCF processing. The `inputs.json` file can be modified to perform different data processing operations.

### Sentieon DNAscope

The `inputs.json` file can be modified to run Sentieon DNAscope instead of Sentieon DNAseq. Here are the most common arguments for Sentieon DNAscope:
```json
"sentieon_germline.run_bqsr": false,
"sentieon_germline.calling_algo": "DNAscope",
"sentieon_germline.dnascope_model": "<path/to/file.model>",
"sentieon_germline.calling_algo_xargs": "--pcr_indel_model none",
```

The `run_bqsr` argument should be set to `false`, as BQSR is not recommended with DNAscope.  
The Sentieon DNAscope pipeline requires a platform-specific model file that is passed through the `calling_algo` argument. You can find the updated DNAscope models for your platform in our [Sentieon-models repository](https://github.com/Sentieon/sentieon-models).  
The `calling_algo_xargs` argument can be used to appropriately set the `--pcr_indel_model` argument in DNAscope. This should be `--pcr_indel_model none` for PCR-free samples, and can be empty for samples sequenced with PCR library preparations.

You can find more information on Sentieon's DNAscope pipeline in Sentieon's [DNAscope appnote](https://support.sentieon.com/appnotes/dnascope_ml/).

### gVCF output

By default, both pipelines will output variants in the VCF format. The `inputs.json` file can be updated to output variants in the gVCF format:
```json
"sentieon_germline.output_gvcf": true,
```

gVCF output format is supported in both Sentieon DNAseq and DNAscope.

### Alignment and pre-processing only

By default, the fastqToVcf pipeline will implement the full Sentieon DNAseq pipeline from Fastq->VCF. The `inputs.json` file can be updated to skip variant calling:
```json
"sentieon_germline.run_calling": false,
```
The pipeline will output an aligned and duplicate-marked BAM/CRAM file along with a BQSR table for the sample.

It is also possible to skip duplicate marking and metrics collection in addition to variant calling so that only alignment with Sentieon BWA is performed:
```json
"sentieon_germline.run_calling": false,
"sentieon_germline.run_bqsr": false,
"sentieon_germline.run_dedup_and_qc": false,
```

### BAM realignment and pre-processing

By default, the BamToVcf pipeline will only perform variant calling from an aligned BAM/CRAM file. The following arguments can be set to re-align an input BAM file to the provided reference genome and perform duplicate marking and BQSR:
```json
"sentieon_germline.realign_input": false,
"sentieon_germline.run_dedup_and_qc": false,
"sentieon_germline.run_bqsr": false,
"sentieon_germline.run_calling": false,
```
