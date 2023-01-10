# ![SAVANA](/docs/SAVANA_logo_transparent.png)

SAVANA is a somatic structural variant (SV) caller for long-read data, currently in Beta testing. Briefly, it takes aligned tumour and normal BAM files, examines the reads for evidence of SVs, clusters adjacent potential SVs together, and finally calls consensus breakpoints and outputs them in BEDPE and VCF format.

SAVANA has been tested on ONT and PacBio HiFi reads aligned with minimap2 and winnowmap. It requires a Unix-based operating system and has been developed and tested on Linux.

## Contents
* [Installation](#installation)
* [Run SAVANA](#run-savana)
  + [Mandatory Arguments](#mandatory-arguments)
  + [Optional Arguments](#optional-arguments)
* [Troubleshooting](#troubleshooting)
* [License](#license)

## Installation

SAVANA requires Python 3.9 with the following dependencies:
- pysam
- pybedtools
- gitpython

All of which can be installed via conda OR pip:
### Install Dependencies with Conda 
To intall and manage dependencies with conda, create a new environment and install dependencies with the `environment.yml` file:
```
conda env create --name <env> --file environment.yml
```
### Install Dependencies with pip
Alternatively, you can install and manage dependencies with pip using the `requirements.txt` file:
```
pip install -r requirements.txt
```

## Run SAVANA

After cloning the repository, you can run SAVANA with the minumum set of arguments:
```
python savana/savana.py --tumour <tumour-bam> --normal <normal-bam> --outdir <outdir> --ref <ref-fasta>
```
### Mandatory Arguments
Argument|Description
---|---
tumour|Tumour BAM file (must have index in .bai format)
normal|Normal BAM file (must have index in .bai format)
outdir|Output directory (can exist but must be empty)
ref|Full path to reference genome that was used to align the `tumour` and `normal` BAM

### Optional Arugments
Argument|Description
---|---
ref_index|Full path to reference genome fasta index (ref path + ".fai" used by default)
contigs|Contigs/chromosomes to consider (default is all in fai file). Example in `example/chr.hg38.txt`
length|Minimum length SV to consider (default=30)
mapq|Minimum MAPQ of reads to consider (default=5)
buffer|Buffer to add when clustering adjacent potential breakpoints (default=10)
depth|Threshold for number of supporting reads (default=3)
threads|Number of threads to use (default is maximum available)
debug|Optional flag to output extra debugging info and files
validation|Optional VCF file to validate against

## Troubleshooting

SAVANA is currently in Beta development. Please raise a GitHub issue if you encounter issues installing or using it.

## License
**SAVANA is free for academic use only**. If you are not a member of a public funded academic and/or education and/or research institution you must obtain a commercial license from EMBL Enterprise Management GmbH (EMBLEM); please email EMBLEM (info@embl-em.de).
