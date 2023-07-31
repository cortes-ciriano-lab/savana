# ![SAVANA](/docs/SAVANA_logo_transparent.png)

SAVANA is a somatic structural variant (SV) caller for long-read data. It takes aligned tumour and normal BAM files, examines the reads for evidence of SVs, clusters adjacent potential SVs together, and finally calls consensus breakpoints, classifies somatic events, and outputs them in BEDPE and VCF format.

SAVANA has been tested on ONT and PacBio HiFi reads aligned with minimap2 and winnowmap. It requires a Unix-based operating system and has been developed and tested on Linux.

## Contents
* [Installation](#installation)
  + [Install SAVANA with Conda](#install-savana-with-conda)
  + [Install SAVANA from Source](#install-savana-from-source)
  + [Check SAVANA Installation](#check-savana-installation)
* [Run SAVANA](#run-savana)
  + [Mandatory Arguments](#mandatory-arguments)
  + [Optional Arguments](#optional-arguments)
  + [Optional Flags](#optional-flags)
  + [Output Files](#output-files)
* [Advanced Options](#advanced-options)
  + [Alternate Classification Methods](#alternate-classification-methods)
  + [Label Known Variants](#label-known-variants)
  + [Train Custom Model](#train-custom-model)
  + [Re-classify Variants](#re-classify-variants)
* [Troubleshooting](#troubleshooting)
* [License](#license)

## Installation

### Install SAVANA with Conda

The easiest and recommended way to install SAVANA is via conda:
```
conda install -c bioconda savana
```

This will install all dependencies and allow you to use SAVANA on the command-line.

### Install SAVANA from Source

_Alternately_, you can install SAVANA from source (note these steps are not required if you've installed SAVANA via conda)

To install from source, SAVANA requires Python 3.9 with the following dependencies:
- pysam
- pybedtools
- cyvcf2

All of which can be installed via conda __OR__ pip:
#### Install Dependencies with Conda
To intall and manage dependencies with conda, create a new environment and install dependencies (including Python 3.9.6) with the `environment.yml` file:
```
conda env create --name <env> --file environment.yml
```

#### Install Dependencies with pip
If preferred, you can install and manage dependencies with pip instead using the `requirements.txt` file
```
pip install -r requirements.txt
```

#### Clone and Install SAVANA
Once you've installed the required dependencies with conda or pip, you can install SAVANA by cloning this repository, navigating to the main folder, and installing with pip:
```
git clone git@github.com:cortes-ciriano-lab/savana.git
cd savana
python3 -m pip install . -vv
```

### Check SAVANA Installation

You can test that SAVANA was installed successfully by running `savana --help`, which should display the following text:
```
usage: savana [-h] [--version] {run,classify,evaluate,train} ...

SAVANA - somatic SV caller

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

subcommands:
  {run,classify,evaluate,train}
                        SAVANA sub-commands
    run                 run SAVANA on tumour and normal long-read BAMs to detect SVs
    classify            classify VCF using model
    evaluate            label SAVANA VCF with somatic/germline/missing given VCF(s) to compare against
    train               train model on folder of input VCFs
```

## Run SAVANA

After installing, SAVANA can be with a minumum set of arguments:
```
savana --tumour <tumour-bam> --normal <normal-bam> --outdir <outdir> --ref <ref-fasta>
```

### Mandatory Arguments
Argument|Description
--------|-----------
tumour|Tumour BAM file (must have index in .bai format)
normal|Normal BAM file (must have index in .bai format)
outdir|Output directory (can exist but must be empty)
ref|Full path to reference genome that was used to align the `tumour` and `normal` BAM

### Optional Arguments
Argument|Description
--------|-----------
ref_index|Full path to reference genome fasta index (ref path + ".fai" used by default)
contigs|Contigs/chromosomes to consider (default is all in fai file). Example in `example/contigs.chr.hg38.txt`
length|Minimum length SV to consider (default=30)
mapq|Minimum MAPQ of reads to consider (default=5)
buffer|Buffer to add when clustering adjacent (non-insertion) potential breakpoints (default=10)
insertion_buffer|Buffer to add when clustering adjacent insertion potential breakpoints (default=100)
depth|Minumum number of supporting reads from tumour OR normal to consider variant (default=3)
threads|Number of threads to use (default is maximum available)
sample|Name to prepend to output files (default=tumour BAM filename without extension)

### Optional Flags
Argument | Description
-------- | -----------
debug | Optional flag to output extra debugging info and files
| ont | Flag to indicate that the Oxford Nanopore (ONT) trained model should be used to classify variants (default) |
| ont_noisy | Flag to indicate that a model trained on ONT data with relatively more noise should be used |
| predict_germline | Flag to indicate that a model that also predicts germline events should be used (a note that this reduced the accuracy of the somatic calls)|

### Output Files

#### Raw SV Breakpoints VCF

`{sample}_sv_breakpoints.vcf` contains all (unfiltered) variants with each edge of a breakpoint having a line in the VCF (i.e. all breakpoints have two lines except for insertions). A note that the `SV_TYPE` field in the `INFO` column of SAVANA VCF output files only denotes `BND` and `INS` types. We have elected not to call `SV_TYPE` beyond these types as it is not definitively possible to do so without copy number information in VCF v4.2 (GRIDSS has a more in-depth explanation for their decision to do this here: https://github.com/PapenfussLab/gridss#why-are-all-calls-bnd).

SAVANA reports the breakend orientation using brackets in the ALT field as described in section 5.4 of [VCFv4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf). We also report this in a `BP_NOTATION` field which can be converted to different nomenclatures as follows:

| Nomenclature | Deletion-like | Duplication-like | Head-to-head Inversion | Tail-to-tail Inversion |
| ------------ | ------------- | ---------------- | ---------------------- | ---------------------- |
| BP_NOTATION  | +- | -+ | ++ | -- |
| Brackets (VCF)  | N[chr:pos[ / ]chr:pos]N | ]chr:pos]N / N[chr:pos[ | N]chr:pos] / N]chr:pos] | [chr:pos[N / [chr:pos[N |
| 5' to 3'  | 3to5 | 5to3 | 3to3 | 5to5 |

SAVANA also reports information about each structural variant in the INFO field of the VCF:

| Field | Description |
| ----- | ----------- |
| SVTYPE | Type of structural variant (INS or BND) |
| MATEID | ID of mate breakend |
| NORMAL_SUPPORT | Number of variant-supporting normal reads |
| TUMOUR_SUPPORT | Number of variant-supporting tumour reads |
| SVLEN | Length of the SV (always >= 0) |
| BP_NOTATION | Notation of breakpoint from table above (_not_ flipped for mate breakpoint) |
| ORIGINATING_CLUSTER | Originating cluster id supporting variant (for debugging) |
| END_CLUSTER | End cluster id supporting variant (for debugging) |
| {ORIGIN\|END}_STARTS_STD_DEV | Cluster value for the standard deviation of the supporting breakpoints' starts |
| {ORIGIN\|END}_MAPQ_MEAN | Cluster value for the mean mapping quality (MAPQ) of the supporting reads |
| {ORIGIN\|END}_EVENT_SIZE_STD_DEV | Cluster value for the standard deviation of the supporting breakpoints' lengths |
| {ORIGIN\|END}_EVENT_SIZE_MEDIAN | Cluster value for the median of the supporting breakpoints' lengths |
| {ORIGIN\|END}_EVENT_SIZE_MEAN | Cluster value for the mean of the supporting breakpoints' lengths |
| {ORIGIN\|END}_TUMOUR_DP | Total depth/coverage (number of reads) in the tumour at SV location (one per breakpoint edge) |
| {ORIGIN\|END}_NORMAL_DP | Total depth/coverage (number of reads) in the normal at SV location (one per breakpoint edge) |

#### Classified Breakpoints VCF

By default, SAVANA classifies somatic variants using a random-forest classifier, trained on a range of somatic Oxford Nanopore data labelled with true somatic variants (as determined by supporting Illumina data). They can be found in the `{sample}.classified.sv_breakpoints.somatic.vcf`. We have found this yields the best results when running on somatic Oxford Nanopore data.

#### Raw Breakpoints BEDPE

`{sample}_sv_breakpoints.bedpe` contains the paired end breakpoints of all unfiltered variants along with their variant ID, length, cluster IDs (for debugging purposes), number of supporting reads from the tumour and normal (listed as "TUMOUR_x/NORMAL_y" - absence indicates 0), and breakpoint orientation (as listed in the table above).

#### Read-support TSV

`{sample}_sv_breakpoints_read_support.tsv` contains one line per structural variant with the variant ID in the first column, the comma-separated ids of the tumour-supporting reads in the second, and normal-supporting reads in the third.

## Alternate Classification Methods

By default, SAVANA uses a model, trained on a range of somatic data. However you may also use alternate classification methods.

### Classify by Parameters File

Given a custom parameters file, you can create your own filters via a JSON file. An example of which can be found in `example/classification-parameters.json`. See below:
```
{
        "somatic": {
                "MAX_NORMAL_SUPPORT": 0,
                "MIN_TUMOUR_SUPPORT": 10,
                "MAX_ORIGIN_STARTS_STD_DEV": 10,
                "MAX_ORIGIN_EVENT_SIZE_STD_DEV": 5
        },
        "germline": {
                "MIN_NORMAL_SUPPORT": 5,
                "MIN_TUMOUR_SUPPORT": 3
        }
}
```
Briefly, you can set limits on the minimum and maximum allowable values for different fields, listed below:
| Field | Description | Type |
| ----- | ----------- | ---- |
| TUMOUR_SUPPORT | Number of variant-supporting tumour reads | Int |
| NORMAL_SUPPORT | Number of variant-supporting normal reads | Int |
| TUMOUR_DP_0 | Total depth/coverage (number of reads) in the tumour at first breakpoint edge | Int |
| TUMOUR_DP_1 | Total depth/coverage (number of reads) in the tumour at second breakpoint edge | Int |
| NORMAL_DP_0 | Total depth/coverage (number of reads) in the normal at first breakpoint edge | Int |
| NORMAL_DP_1 | Total depth/coverage (number of reads) in the normal at second breakpoint edge |  Int |
| TUMOUR_SUPPORT_RATIO | Ratio of tumour-supporting reads to DP (averaged over both breakpoints) | Float |
| NORMAL_SUPPORT_RATIO | Ratio of normal-supporting reads to DP (averaged over both breakpoints) | Float |
| {ORIGIN\|END}_STARTS_STD_DEV | Cluster value for the standard deviation of the supporting breakpoints' starts | Float |
| {ORIGIN\|END}_MAPQ_MEAN | Cluster value for the mean mapping quality (MAPQ) of the supporting reads | Float |
| {ORIGIN\|END}_EVENT_SIZE_STD_DEV | Cluster value for the standard deviation of the supporting breakpoints' lengths | Float |
| {ORIGIN\|END}_EVENT_SIZE_MEDIAN | Cluster value for the median of the supporting breakpoints' lengths | Float |
| {ORIGIN\|END}_EVENT_SIZE_MEAN | Cluster value for the mean of the supporting breakpoints' lengths | Float |

### Classify by Legacy Methods

Alternately, you can use the `--legacy` flag to use filtering and classification methods used in the Beta version of SAVANA. This will output `strict` and `lenient` somatic VCF files which are informed by a decision-tree classifier (strict) and manually plotting data to determine cutoffs (lenient).

## Advanced Options

### Label Known Variants

If you have known somatic and germline variants in a VCF that you'd like to annotate in the output of SAVANA, you can provide them via the `--somatic` and `--germline` command-line arguments (N.B. you cannot provide germline variants without providing somatic ones). This will output a `{sample}.evaluation.sv_breakpoints.vcf` which contains the classified variants (if a model was used) or all variants (if custom filters or legacy methods were used) with a `LABEL` added to the `INFO` field in the VCF which indicates whether a variant was found in the `SOMATIC` or `GERMLINE` files or was `NOT_IN_COMPARISON`.

By default, a buffer of 100bp is used to consider two variants as overlapping. This can be modified via the `--overlap_buffer` command-line argument. Statistics about the overlapping variants are automatically written to a `{sample}.evaluation.stats` file, the name of which can be overwritten using the `--stats` argument. Tie-breakers (when two variants are both within the overlap window to a known variant) are by default broken by distance, with the closest variant being used. Optionally, you can also tie-break based on `SUPPORT` - e.g.) if two variants are within the overlap window to a known somatic variant, use the variant with the highest `TUMOUR_SUPPORT` (and vice versa for a germline variant and `NORMAL_SUPPORT`).

If you decide you want to label variants after SAVANA has already been run, you can do so via the sub-command `savana evaluate` like so:
```
savana evaluate --input ${savana_vcf} --somatic ${known_somatic_variants_vcf} --output ${savana_labelled_vcf}
```

See the table below for a full list of arguments:
| Argument | Description |
| -------- | ----------- |
| input | VCF file to evaluate |
| somatic | Somatic VCF file to evaluate against |
| germline | Germline VCF file to evaluate against (optional) |
| overlap_buffer | Buffer for considering an overlap (default=100) |
| output | Output VCF with LABEL added to INFO |
| stats | Output file for statistics on comparison if desired |
| by_support | Flag for comparison method: tie-break by read support |
| by_distance | Flag for comparison method: tie-break by min. distance (default) |

### Train Custom Model

SAVANA also provides functionality to train your own random forest model using raw labelled VCFs from the previous step.

An example of how to train a model:
```
savana train --vcfs ${folder_of_labelled_vcfs} --outdir ${output_directory_for_model}
```

The output directory contains a `.pkl` file which can be used to classify variants via the `--model` argument to SAVANA. If you'd like to classify an existing output file, you can do so via the `savana classify` sub-command (see next section, [Re-classify](#re-classify-variants)).

Additional output files include a `model_arguments.txt` file with the savana command used to train the model; a `model_stats.txt` file with the precision, recall, f1-score and number of variants in each class (0 is false);  a `confusion_matrix.png` of the TP/FP/TN/FN breakdown in the test set; and `test_set_incorrect.tsv` and `test_set_correct.tsv` files which list the variants from the test set (20% of the input by default - 80% is used for training) that were incorrectly and correctly categorized, along with their information.

See the table below for a full list of options and arguments to the `savana train` sub-command:
| Argument | Description |
| -------- | ----------- |
| vcfs | Folder of labelled VCF files to read in |
| recursive | Set flag to search recursively through input folder for input VCFs (default only one-level deep) |
| load_matrix | Pickle file of pre-processed VCFs (faster) |
| save_matrix | Optional output pickle file of processed VCFs (to be used in `load_matrix` argument for faster loading) |
| downsample | Fraction to downsample the majority class by - with 0 removing no data and .99 removing 99% of it (default=0.1) |
| germline_class | Train the model to predict germline and somatic variants (GERMLINE label must be present) |
| hyper | Perform a randomised search on hyper parameters and use best |
| outdir | Output directory (can exist but must be empty)


### Re-classify Variants

If you'd like to re-classify variants using an alternate method after SAVANA has already been run, you can do so via the `savana classify` sub-command. An example of reclassifying an existing savana output VCF using a custom model would be:
```
savana classify --vcf {raw_sv_breakpoints_vcf} --model {custom_trained_model_pkl} --output {output_classified_vcf}
```

You can also use the `savana classify` sub-command to re-classify using custom parameters (see [Classify by Parameters File](#classify-by-parameters-file)), or legacy methods ([Classify by Legacy Methods](#classify-by-legacy-methods)). See the table below for a full list of arguments:
| Argument | Description |
| -------- | ----------- |
| input | VCF file to classify |
| ont | Flag to indicate that the Oxford Nanopore (ONT) trained model should be used to classify variants (default) |
| ont_noisy | Flag to indicate that a model trained on ONT data with relatively more noise should be used |
| predict_germline | Flag to indicate that a model that also predicts germline events should be used (a note that this reduced the accuracy of the somatic calls)|
| model | Pickle file of machine-learning model |
| custom_params | JSON file of custom filtering parameters |
| legacy | Use legacy lenient/strict filtering |

## Troubleshooting

Please raise a GitHub issue if you encounter issues installing or using SAVANA.

## License
**SAVANA is free for academic use only**. If you are not a member of a public funded academic and/or education and/or research institution you must obtain a commercial license from EMBL Enterprise Management GmbH (EMBLEM); please email EMBLEM (info@embl-em.de).
