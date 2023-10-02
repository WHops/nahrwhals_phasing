# Nahrwhals Phasing 

Phase a VCF into 1000 genomes SNPs by NYCG.

## Table of Contents

1. [Description](#Description)
2. [Workflow](#Workflow)
3. [Post-processing Steps](#Post-processing-Steps)
4. [Usage](#Usage)

## Description

The Snakefile takes as input a series of phased SNP VCFs (one per chromosome), and the NW output folders. It checks all aligned reads, assigns them a phase, and determines statistics for the most likely phase. 

Final output, at the moment, is `actions_all.tsv`, a table with instructions. There is an `execute_actions` function, but it is not included in the Snakemake workflow so far. 

<img src="https://github.com/WHops/nahrwhals_phasing/blob/main/dag_one_region.png?raw=true">

## Workflow

Given the different post-processing steps, running NWs in its entirety requires some steps: 

### Initial Steps
1. Run the `ntk_interesting` pipeline which repeatedly calls Nahrwhals:
    ```shell
    ssh seneca
    cd /g/korbel/hoeps/projects/nahr/ntk_interesting_loci/ntk-scan-snakemake
    ## Edit configs ##
    conda activate nahrwhals-big
    ```

   This will produce an `all.tsv`.

2. Convert `all.tsv` to `all.vcf` using:
    ```shell
    R
    nahrwhals::write_vcf('/path/to/all.tsv', '/path/to/all.vcf')
    ```

### Phasing with External Reference
If you also want to phase with respect to an external phasing:

3. Run the `nahrwhals_phasing` pipeline:
    ```shell
    cd /g/korbel/hoeps/projects/nahr/phaselab/snake_approach
    ## Edit configs to match those in step 1 and link to all.vcf ##
    conda activate nahrwhals-big
    snakemake --cores n
    ```

   This process will produce an externally phased VCF, named `NW_flipped_50000.vcf`.


## Usage

Edit the `config/Snake.config.json` to indicate binaries for a lot of tools, and links to the reference fasta, NYGC VCF, and NAHRwhals output path. Files that have to be custom-made for this workflow are: 

- `region_file`:
    ```
    chr1-103543499-103821114
    chr1-108148088-108537579
    chr1-109662781-109726958
    ...
    ```

- `sample_file`:
    ```
    GM19129
    GM19434
    HG00171
    ...
    ```

Ensure you have a working R environment loaded. If on Seneca, use:
    ```shell
    conda deactivate # (a few times if necessary)
    conda activate nahrwhals-big
    ```
x
Run the workflow with:

    ```
    snakemake --cores n
    ```

## Post-processing Steps

### Merging BCF Files
After running the Snakefile, you need to post-process the results to merge the BCF files:

1. First, convert the `NW_flipped_50000.vcf` produced from the Snakemake part to a BCF file:
    ```shell
    bcftools view -O b NW_flipped_50000.vcf -o NW_flipped_50000.bcf
    ```
2. Concatenate the individual `allsample.vcfs` produced from the Snakemake using `bedtools` (if you have multiple):
    ```shell
    bedtools concat -i input1.vcf,input2.vcf -o allregions_allsamples.vcf.gz
    ```
3. Run the BCF merging script with the generated BCF files:
    ```shell
    ./merge_bcf.sh NW_flipped_50000.bcf allregions_allsamples.vcf.gz merged_output.vcf.gz
    ```
