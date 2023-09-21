# Nahrwhals Phasing 

Phase a vcf into 1000 genomes snps by NYCG.

## Table of Contents

1. [Description](#Description)
2. [Usage](#Usage)

## Description

The Snakefile takes as input a series of phased SNP vcfs (one per chromsoome), and the NW output folders. 
It then checks all aligned reads, assigns them a phase, and determines statistics for the most likely phase. 

Final output, at the moment, is actions_all.tsv, a table with instructions. There is an execute_actions somewhere but it's not included in the snakething so far. 

<img src="https://github.com/WHops/nahrwhals_phasing/blob/main/dag_one_region.png?raw=true">

## Workflow of a phased Nahrwhals run.

Given the different post-processing steps, running NWs in its entirety requires some steps. 

1) Run the ntk_interesting pipeline which repeatedly calls nahrwhals
   ```
   ssh seneca
   cd /g/korbel/hoeps/projects/nahr/ntk_interesting_loci/ntk-scan-snakemake
   ## Edit configs ##
   conda activate nahrwhals-big
   ## Ups, for the new version of NW we still need to make a new snakefile ##
   ```

    Now we get an all.tsv

2) Convert all.tsv to all.vcf using
   ```
   R
   nahrwhals::write_vcf('/path/to/all.tsv', '/path/to/all.vcf')
   ```

Now we are essentially done. However, if we also want to phase with respect to an outside phasing:

3) Run the nahrwhals_phasing pipeline
   ```
   cd /g/korbel/hoeps/projects/nahr/phaselab/snake_approach
   ## Edit configs to resemble your configs in step 1 and it links to all.vcf ##
   conda activate nahrwhals-big
   snakemake --cores n
   ```

Now you have a externally phased vcf.


## Usage

Edit the config/Snake.config.json to indicate binaries for a lot of tools, and links to the reference fasta, NYGC vcf and NAHRwhals output path. Files that have to be custom-made for this workflow are: 

- region_file
```
chr1-103543499-103821114
chr1-108148088-108537579
chr1-109662781-109726958
...
```

- sample_file
```
GM19129
GM19434
HG00171
...
```

Make sure that you have a working R environment loaded. For me personally, on seneca this is
```
conda deactivate (a few times...)
conda activate nahrwhals-big
```

Run the workflow with 

```
snakemake --cores n
```

