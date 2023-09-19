# Nahrwhals Phasing 

Phase a vcf into 1000 genomes snps by NYCG.

## Table of Contents

1. [Description](#Description)
2. [Usage](#Usage)

## Description

<img src="https://github.com/WHops/nahrwhals_phasing/blob/main/dag_one_region.png?raw=true">

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

Run the workflow with 

```
snakemake --cores n
```

