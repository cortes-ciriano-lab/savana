# ![SAVANA](/docs/SAVANA_logo_transparent.png)

SAVANA is a somatic structural variant (SV) and copy number aberrations (SCNA) caller for long-read data. It takes aligned tumour and normal BAM files, examines the reads for evidence of SVs, clusters adjacent potential SVs together, and finally calls consensus breakpoints, classifies somatic events, and outputs them in BEDPE and VCF format. It also identifies copy number abberations utilising somatic SV breakpoints and circular binary segmentation. SAVANA then estimates tumour purity using B-allele frequency values of heterozygous SNPs at regions with loss of heterozygosity, and performs absolute copy number fitting to determine tumour ploidy and allele-specific absolute copy number.

SAVANA has been tested on ONT and PacBio HiFi reads aligned with minimap2 and winnowmap. It requires a Unix-based operating system and has been developed and tested on Linux.

For further information, benchmarking and for citation, please refer to our [SAVANA preprint](https://doi.org/10.1101/2024.07.25.604944).

## Contents
* [Installation](#installation)
  + [Install SAVANA with Conda](#install-savana-with-conda)
  + [Alternately, Install SAVANA from Source](#alternately-install-savana-from-source)
  + [Check SAVANA Installation](#check-savana-installation)
* [Run SAVANA](#run-savana)
  + [Mandatory Arguments](#mandatory-arguments)
  + [Optional Arguments](#optional-arguments)
  + [Tumour-only Mode](#tumour-only-mode)
* [Output Files](#output-files)
  + [Output Files SV Algorithm](#output-files-sv-algorithm)
  + [Output Files CNA Algorithm](#output-files-cna-algorithm)
* [Phasing Information](#phasing-information)
  + [Generating Phased VCF](#generating-phased-vcf)
  + [Generating Phased BAMs](#generating-phased-bams)
* [Advanced Options](#advanced-options)
  + [Alternate Classification Methods](#alternate-classification-methods)
  + [Label Known Variants](#label-known-variants)
  + [Train Custom Model](#train-custom-model)
  + [Re-classify Variants](#re-classify-variants)
* [Troubleshooting](#troubleshooting)
* [License](#license)
* [Contacts](#contacts)

## Installation

### Install SAVANA with Conda

The easiest and recommended way to install SAVANA is via conda:
```
conda install -c bioconda savana
```

This will install all dependencies and allow you to use SAVANA on the command-line.

### Alternately, install SAVANA from Source

_Alternately_, you can install SAVANA from source (note these steps are not required if you've installed SAVANA via conda)

First, clone this repository:
```
git clone git@github.com:cortes-ciriano-lab/savana.git
```

To install from source, SAVANA requires Python 3.9 with dependencies as listed in the `requirements.txt` file

All of which can be installed via conda __OR__ pip:
#### Install Dependencies with Conda
To intall and manage dependencies with conda, create a new environment and install dependencies (including Python 3.9.6) with the `environment.yml` file in the top-level of the repository:
```
conda env create --name <env> --file environment.yml
```

#### Install Dependencies with pip
If preferred, you can install and manage dependencies with pip instead using the `requirements.txt` file
```
pip install -r requirements.txt
```

#### Install SAVANA
Once you've installed the required dependencies with conda or pip, you can install SAVANA by navigating to the cloned repo and running:
```
python3 -m pip install . -vv
```

### Check SAVANA Installation

You can test that SAVANA was installed successfully by running `savana --help`, which should display the following text:
```
usage: savana [-h] [--version] {run,classify,evaluate,train,cna,to} ...

SAVANA - somatic SV caller

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

subcommands:
  {run,classify,evaluate,train,cna,to}
                        SAVANA sub-commands
    run                 identify and cluster breakpoints - output raw variants without classification
    classify            classify VCF using model
    evaluate            label SAVANA VCF with somatic/germline/missing given VCF(s) to compare against
    train               train model on folder of input VCFs
    cna                 run copy number
    to                  identify somatic SVs without normal sample
```

## Run SAVANA

After installing, SAVANA can be run on long-read data with a minumum set of arguments:
```
savana --tumour <tumour-file> --normal <normal-file> --outdir <outdir> --ref <ref-fasta>
```

This will call somatic SVs. We recommend running with the --contigs argument (`--contigs example/contigs.chr.hg38.txt`) to only examine chromosomes of interest. We recommend phasing your input BAM files for optimal results (see [Phasing Information](#phasing-information)).

To compute copy number aberrations, you must provide a SNP VCF for the germline sample. Then, to call both SVs and CNAs you can run savana with:
```
savana --tumour <tumour-file> --normal <normal-file> --outdir <outdir> --ref <ref-fasta> --snp_vcf <vcf-file>
```
Alternatively, you can use the 1000 genome vcf files (provided within SAVANA) using `--g1000_vcf 1000g_hg38` for hg38 or `--g1000_vcf 1000g_hg19` for hg19 aligned bam files. However, we strongly recommend using a SNP VCF from the matched germline sample for best performance.

If you would like to provide a blacklist for CNA calling you can do so via:
```
 --blacklist <blacklist-bed-file>
```

If you have already generated the heterozygous SNP allele counts using the above command, you can skip this step by providing the <allele_counts_hetSNPs.bed> using the `--allele_counts_het_snps` parameter instead of a VCF (`--snp_vcf/--g1000_vcf`).
```
--allele_counts_het_snps <allele_counts_hetSNPs.bed>
```

### Quickstart

To test on a small BAM, download subset COLO829 cell-line files:
```
wget ftp://ftp.ebi.ac.uk/pub/databases/icortes-public/COLO829_subset/*.bam*
```

Then, to call somatic SVs, run:
```
savana --tumour ONT_COLO829_T_truthset_50k.phased.bam --normal ONT_COLO829_N_truthset_50k.phased.bam --outdir COLO829_subset_test --threads 8 --ref <hg38-fasta>
```
> To run the above on a SLURM cluster we requested an interactive job with 8 cpus (`-c 8`) and 16 GB of memory (`--mem 16`). The total time to call variants (as measured by Linux `time`) was 48.88 seconds. You may need to adjust requirements for your computing environment


### Tumour-only Mode

We strongly recommend running SAVANA conventionally using tumour and matched normal bam files for best performance. However, we have developed a SAVANA tumour-only `savana to` mode which can be run in the absence of a matched normal sample if required:

```
savana to --tumour <tumour-file> --outdir <outdir> --ref <ref-fasta> --g1000_vcf <vcf-file>
```

In the absence of a normal BAM, we still recommend phasing your tumour BAM file using population SNPs as this will produce optimal SV calls (see [Phasing Information](#phasing-information) for more details on how to do this).

If you use this mode, filtering the resulting SVs using external population and panel of normal (PoN) resources is **highly** recommended. Specifically, we recommend removing SVs that overlap with any SVs in [gnomadSV](https://gnomad.broadinstitute.org/data#v4-structural-variants) that have >=10% population allele-frequency (AF in `gnomad.v4.1.sv.sites.bed.gz`) as well as SVs present in the Hartwig Medical Foundation SV PoN (available at [hmf_dna_pipeline_resources](https://storage.googleapis.com/hmf-public/HMFtools-Resources/pipeline/oncoanalyser/2.0/38/hmf_pipeline_resources.38_v2.0.0--3.tar.gz) - `dna/sv/sv_pon.38.bedpe.gz` in the untarred directory).


### SAVANA Arguments

#### Mandatory Arguments
Argument|Description
--------|-----------
--tumour|Tumour BAM/CRAM file (must have index in .bai/.crai format)
--normal|Normal BAM/CRAM file (must have index in .bai/.crai format)
--outdir|Output directory (can exist but must be empty)
--ref|Full path to reference genome that was used to align the `tumour` and `normal` BAM

#### Optional Arguments
Argument|Description
--------|-----------
*Basic Arguments*
--snp_vcf| Path to SNP vcf file to extract heterozygous SNPs for allele counting
--g1000_vcf | Use 1000g biallelic vcf file for allele counting instead of SNP vcf from matched normal. Specify which genome version to use. choices={"1000g_hg38", "1000g_hg19", "1000g_t2t"}
--ont | Run on Nanopore data (default)
--pb | Use PacBio filters to classify variants ([Classify for PacBio](#classify-for-pacbio))
--sample| Name to prepend to output files (default=tumour BAM filename without extension)
--contigs| Contigs/chromosomes to consider (default is all in fai file). Example in `example/contigs.chr.hg38.txt`. Should be in order.
--length| Minimum length SV to consider (default=30)
--keep_inv_artefact| Do not remove breakpoints with foldback-inversion artefact pattern (default is to remove)
--mapq| Minimum MAPQ of reads to consider (default=0)
--min_support| Minimum supporting reads for a variant (default=3)
--min_af| Minimum allele-fraction (AF) for a variant (default=0.01)
--cna_resuce| Copy number abberation output file for this sample (used to rescue variants)
--cna_rescue_distance| Maximum distance from a copy number abberation for a variant to be rescued by it
--threads| Number of threads to use (default is maximum available)
--cna_threads| Number of threads to use for CNA calling (default is maximum available)
--ref_index| Full path to reference genome fasta index (ref path + ".fai" used by default)
--single_bnd| Report single breakend variants in addition to standard types (False by default)
--single_bnd_min_length| Minimum length of single breakend to consider (default=100)
--single_bnd_max_mapq| Convert supplementary alignments below this threshold to single breakend (default=5, must not exceed --mapq argument)
--confidence| If using a model to classify variants (default), you can use Mondrian Conformal Prediction with a chosen confidence (0.01-0.99)
*SV Algorithm Arguments*
--buffer| Buffer to add when clustering adjacent (non-insertion) potential breakpoints, excepting insertions (default=10)
--insertion_buffer| Buffer to add when clustering adjacent insertion potential breakpoints (default=100)
--end_buffer | Buffer to add when clustering the alternate edge of potential breakpoints, excepting insertions (default=100)
--coverage_binsize | Length used for coverage bins (default=5)
--min_reads_per_cluster | During initial clustering, discard clusters with fewer than n (default=3) reads supporting a putative event of any size/orientation
--chunksize | Chunksize to use when splitting genome for parallel analysis (default=1000000)
*CNA Algorithm Arguments*
--tmpdir| Temp directory for allele counting temp files (defaults to outdir)
--allele_counts_het_snps| If allele counting has already been performed provide the path for the allele counts of heterozygous SNPs to skip this step
--allele_mapq|    Mapping quality threshold for reads to be included in the allele counting (default = 5)
--allele_min_reads|    Minimum number of reads required per het SNP site for allele counting (default = 10)
--ac_window| Window size for allele counting to parallelise (default = 1200000; this should be >=500000)
--cn_binsize|  Bin window size in kbp (default=10)
--blacklist| Path to the blacklist file
--breakpoints| Path to SAVANA VCF file to incorporate savana breakpoints into copy number analysis
--chromosomes| Chromosomes to analyse. To run on all chromosomes leave unspecified (default). To run on a subset of chromosomes only specify the chromosome numbers separated by spaces. For x and y chromosomes, use 23 and 24, respectively.  E.g. use "--chromosomes 1 4 23 24" to run chromosomes 1, 4, X and Y
--readcount_mapq| Mapping quality threshold for reads to be included in the read counting (default = 5)
--no_blacklist| Don't use a blacklist
--bl_threshold| Percentage overlap between bin and blacklist threshold to tolerate for read counting (default = 5). Please specify percentage threshold as integer, e.g. "-t 5". Set "-t 0" if no overlap with blacklist is to be tolerated
--no_basesfilter| Do not filter bases
--bases_threshold| Percentage of known bases per bin required for read counting (default = 75). Please specify percentage threshold as integer, e.g. "-bt 95"
--smoothing_level|   Size of neighbourhood for smoothing (default = 10)
--trim| Trimming percentage to be used (0.025)
--min_segment_size|    Minimum size for a segement to be considered a segment (default = 5)
--shuffles|    Number of permutations (shuffles) to be performed during CBS (default = 1000)
--p_seg| p-value used to test segmentation statistic for a given interval during CBS using (shuffles) number of permutations (default = 0.05)
--p_val| p-value used to test validity of candidate segments from CBS using (shuffles) number of permutations (default = 0.01)
--quantile| Quantile of changepoint (absolute median differences across all segments) used to estimate threshold for segment merging (default = 0.2; set to 0 to avoid segment merging)
--min_ploidy| Minimum ploidy to be considered for copy number fitting (default = 1.5)
--max_ploidy| Maximum ploidy to be considered for copy number fitting (default = 5)
--ploidy_step| Ploidy step size for grid search space used during for copy number fitting (default = 0.01)
--min_cellularity| Minimum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided this is estimated during copy number fitting. Alternatively a purity value can be provided if the purity of the sample is already known
--max_cellularity| Maximum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided this is estimated during copy number fitting. Alternatively a purity value can be provided if the purity of the sample is already known
--cellularity_step| Cellularity step size for grid search space used during for copy number fitting (default = 0.01)
--cellularity_buffer| Cellularity buffer to define purity grid search space during copy number fitting (default = 0.1)
--overrule_cellularity| Set to sample`s purity if known. This value will overrule the cellularity estimated using hetSNP allele counts (not used by default)
--distance_function| Distance function to be used for copy number fitting. choices=[RMSD, MAD] (default = RMSD)
--distance_filter_scale_factor| Distance filter scale factor to only include solutions with distances < scale factor * min(distance)
--distance_precision|   Number of digits to round distance functions to (default = 3)
--max_proportion_zero| Maximum proportion of fitted copy numbers to be tolerated in the zero or negative copy number state (default = 0.1)
--min_proportion_close_to_whole_number| Minimum proportion of fitted copy numbers sufficiently close to whole number to be tolerated for a given fit (default = 0.5)
--max_distance_from_whole_number| Distance from whole number for fitted value to be considered sufficiently close to nearest copy number integer (default = 0.25)
--main_cn_step_change| Max main copy number step change across genome to be considered for a given solution
--min_block_size|   Minimum size (number of SNPs) for a genomic block to be considered for purity estimation (default = 10)
--min_block_length|   Minimum length (bps) for a genomic block to be considered for purity estimation (default = 500000)
*Additional Arguments*
--legacy | Use legacy filters (strict/lenient) to classify variants
--custom_model | Path to custom model pkl file
--custom_params | Path to custom paramaters JSON file for filtering
--somatic_output | Output VCF path for a file containing only PASS somatic variants
--germline_output | Output VCF path for a file containing only PASS germline variants (`--predict_germline` must be specified)
--somatic | VCF file containing somatic variants to evaluate PASS somatic variants against
--germline | VCF file containing germline variants to evaluate PASS germline variants against
--overlap_buffer | If comparing against a --somatic or --germline VCF file, buffer for considering variants to overlap
--by_support | When comparing to --somatic or --germline VCF, tie-break by read-support
--by_distance | When comparing to --somatic or --germline VCF, tie-break by distance (default)
--stats | Output filename for statistics on comparison to somatic/germline VCF

## Output Files
### Output Files SV Algorithm

#### Somatic Breakpoints

By default, SAVANA classifies somatic variants using a random-forest classifier, trained on a range of somatic Oxford Nanopore data labelled with true somatic variants (as determined by supporting Illumina data). We have found that this yields the best results when running on somatic Oxford Nanopore data. The SVs classified as somatic can be found in the `{sample}.classified.somatic.vcf` as well as a `{sample}.classified.somatic.bedpe`.

##### Somatic VCF

In SAVANA, each edge of a breakpoint has a line in the VCF (i.e. all breakpoints have two lines except for insertions). The `SV_TYPE` field in the `INFO` column of SAVANA VCF output files only denotes `BND` and `INS` types. We have elected not to call `SV_TYPE` beyond these types as it is not definitively possible to do so without copy number information in VCF v4.2 (GRIDSS has a more in-depth explanation for their decision to do this here: https://github.com/PapenfussLab/gridss#why-are-all-calls-bnd).

SAVANA reports the breakend orientation using brackets in the ALT field as described in section 5.4 of [VCFv4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf). We also report this in a `BP_NOTATION` field which can be converted to different nomenclatures as follows:

| Nomenclature | Deletion-like | Duplication-like | Head-to-head Inversion | Tail-to-tail Inversion |
| ------------ | ------------- | ---------------- | ---------------------- | ---------------------- |
| BP_NOTATION  | +- | -+ | ++ | -- |
| Brackets (VCF)  | N[chr:pos[ / ]chr:pos]N | ]chr:pos]N / N[chr:pos[ | N]chr:pos] / N]chr:pos] | [chr:pos[N / [chr:pos[N |
| 5' to 3'  | 3to5 | 5to3 | 3to3 | 5to5 |

SAVANA reports information about each structural variant in the INFO field of the VCF:
<details>
<summary>VCF INFO Fields</summary>

| Field | Description |
| ----- | ----------- |
| SVTYPE | Type of structural variant (INS or BND) |
| CLASS | Variant class (predicted or classified as SOMATIC or NOISE)
| MATEID | ID of mate breakend |
| TUMOUR_READ_SUPPORT | Number of variant-supporting tumour reads |
| TUMOUR_ALN_SUPPORT | Number of variant-supporting tumour alignments |
| NORMAL_READ_SUPPORT | Number of variant-supporting normal reads |
| NORMAL_ALN_SUPPORT | Number of variant-supporting normal alignments |
| TUMOUR_AF | Tumour allele-fraction: ratio of tumour-supporting reads to DP (averaged over both edges) |
| NORMAL_AF | Normal allele-fraction: ratio of normal-supporting reads to DP (averaged over both edges) |
| SVLEN | Length of the SV (always >= 0) |
| BP_NOTATION | Notation of breakpoint from table above (_not_ flipped for mate breakpoint) |
| SOURCE | Source of evidence for a breakpoint - CIGAR (INS, DEL, SOFTCLIP), SUPPLEMENTARY or mixture |
| CLUSTERED_READS_NORMAL | Total number of normal reads clustered at this location supporting any SV type |
| CLUSTERED_READS_TUMOUR | Total number of tumour reads clustered at this location supporting any SV type |
| TUMOUR_DP_{BEFORE\|AT\|AFTER} | Local tumour depth in bin before/at/after the breakpoint(s) of an SV |
| NORMAL_DP_{BEFORE\|AT\|AFTER} | Local normal depth in bin before/at/after the breakpoint(s) of an SV |
| TUMOUR_ALT_HP | Counts of SV-supporting reads belonging to each haplotype in the tumour sample (1/2/NA) |
| NORMAL_ALT_HP | Counts of SV-supporting reads belonging to each haplotype in the normal sample (1/2/NA) |
| TUMOUR_PS | List of unique phase sets from the tumour supporting reads |
| NORMAL_PS | List of unique phase sets from the normal supporting reads |
| TUMOUR_TOTAL_HP_AT | Counts of all reads at SV location belonging to each haplotype in the tumour sample (1/2/NA) |
| TUMOUR_TOTAL_HP_AT | Counts of all reads at SV location belonging to each haplotype in the normal sample (1/2/NA) |
| {ORIGIN\|END}_STARTS_STD_DEV | Cluster value for the standard deviation of the supporting breakpoints' starts |
| {ORIGIN\|END}_MAPQ_MEAN | Cluster value for the mean mapping quality (MAPQ) of the supporting reads |
| {ORIGIN\|END}_EVENT_SIZE_STD_DEV | Cluster value for the standard deviation of the supporting breakpoints' lengths |
| {ORIGIN\|END}_EVENT_SIZE_MEDIAN | Cluster value for the median of the supporting breakpoints' lengths |
| {ORIGIN\|END}_EVENT_SIZE_MEAN | Cluster value for the mean of the supporting breakpoints' lengths |
| {ORIGIN\|END}_TUMOUR_DP | Total depth/coverage (number of reads) in the tumour at SV location (one per breakpoint edge) |
| {ORIGIN\|END}_NORMAL_DP | Total depth/coverage (number of reads) in the normal at SV location (one per breakpoint edge) |
</details>

##### Somatic BEDPE

`{sample}.classified.somatic.bedpe` contains the paired end breakpoints of all SVs predicted as somatic along with their variant ID, length, and, number of supporting reads from the tumour and normal (listed as "TUMOUR_x/NORMAL_y" - absence indicates 0), and breakpoint orientation (as listed in the table above).

#### Raw Breakpoints

If you would like to access _all_ SV breakpoints, including those with normal support or predicted as noise by the model, they are in the `{sample}_sv_breakpoints.vcf` and `{sample}_sv_breakpoints.bedpe` files. The INFO columns and structure of the raw VCF is the same as in the [somatic VCF](#somatic-vcf) and the BEDPE as the [somatic BEDPE](#somatic-bedpe).

#### Read-support TSV

`{sample}_sv_breakpoints_read_support.tsv` contains one line per structural variant with the variant ID in the first column, the comma-separated ids of the tumour-supporting reads in the second, and normal-supporting reads in the third.

#### Inserted sequences FASTA

To enable the examination of inserted sequence, `{sample}.inserted_sequences.fa` contains supporting sequences for insertion variants. Each sequence identifier contains the variant ID, index of the inserted sequence (INSSEQ_0000n, INSSEQ_0000n+1, etc.), variant type (whether insertion or single-breakend), index of the insert (out of all inserts supporting the variant, i.e. 1/17), and length of the sequence. This identifier line is followed by the inserted bases.

### Output Files CNA Algorithm

#### Raw read counts TSV
`{sample}_{cn_binsize}_raw_read_counts.tsv` contains all raw and unfiltered read counts for each bin across the reference genome for the tumour and matched normal sample. In addition, SAVANA also outputs other intermediate files during copy number processing, including the filtered and matched normal normalised log2 transformed read counts (`{sample}_{cn_binsize}_read_counts_mnorm_log2r.tsv`).

#### Segmented log2r relative copy number TSV
`{sample}_{cn_binsize}_read_counts_mnorm_log2r_segmented.tsv` contains the final relative copy number (log2r) data post CBS segmentation. This includes the log2r relative copy number for each bin across the reference genome, as well as the segment IDs and according segmented log2r relative copy number values.

#### Fitted purity and ploidy TSV
`{sample}_{cn_binsize}_fitted_purity_ploidy.tsv` contains the final copy number fit (i.e. purity and ploidy values, as well as the distance function used during fitting) for a given sample. Note that SAVANA also outputs all viable solutions together with their distance functions and ranking prior to the final fit being selected (`{sample}_{cn_binsize}_ranked_solutions.tsv`).

#### Segmented absolute copy number
The final and main SAVANA CNA output file is `{sample}_{cn_binsize}_segmented_absolute_copy_number.tsv`, which contains the fitted total and minor absolute copy number values for each copy number segment (collapsed). Note that this output file (together with the classified somatic SAVANA SV calls) can be used to generate the Copy Number ReCon Plots, as outlined and described [here](https://github.com/cortes-ciriano-lab/ReConPlot).

## Phasing Information

### Generating Phased VCF

We recommend using [LongPhase](https://github.com/twolinin/longphase) or [WhatsHap](https://whatshap.readthedocs.io/en/stable/index.html) to generate phased VCF from matched normal samples or a set of population SNPs. As an example, WhatsHap can be run using the following command:

```
whatshap phase  --ignore-read-groups -o <phased.vcf.gz> --reference=<ref-fasta> <germline_snps.vcf> <normal-file>
```

Germline SNPs (`<germline_snps.vcf>`) can for example be obtained using [Clair3](https://github.com/HKU-BAL/Clair3) on the matched normal long-read BAM (or Strelka if using phased SNPs from a matched normal Illumina sample).

### Generating Phased BAMs

We recommend using [LongPhase](https://github.com/twolinin/longphase) or [WhatsHap](https://whatshap.readthedocs.io/en/stable/index.html) to tag sequencing reads by haplotype. This will generate phased BAMs for both the tumour and normal (if available) using the `<phased.vcf.gz>` generated using the steps above.
```
whatshap haplotag --ignore-read-groups -o <phased_tumour.bam> --reference <ref-fasta> <phased.vcf.gz> <tumour_bam> && samtools index <phased_tumour.bam>
whatshap haplotag --ignore-read-groups -o <phased_normal.bam> --reference <ref-fasta> <phased.vcf.gz> <normal_bam> && samtools index <phased_normal.bam>
```

## Advanced Options
### Alternate Classification Methods

By default, SAVANA uses a model, trained on a range of ONT somatic data. However you may also use alternate classification methods. You may also [train your own model](#train-custom-model).

#### Classify for PacBio

Currently, there is no model available in SAVANA which was trained on PacBio data. If the `--pb` flag is used, a set of filters (shown in the table below), will be used. The minimum allele-fraction (AF) and support can be modified with the `--min_support` and `--min_af` flags. By default `--min_support` is set to 5, but we recommend testing different values here - in our PacBio samples, increasing this value to 10 yielded the best results.

| Field | Description | PacBio Somatic Filter |
| ----- | ----------- | --------------------- |
| TUMOUR_SUPPORT | Number of variant-supporting tumour reads; modify with `--min_support` | >=7 |
| TUMOUR_AF | Tumour allele-fraction: ratio of tumour-supporting reads to DP (averaged over both edges); modify with `--min_af` | >=0.15 |
| NORMAL_SUPPORT | Number of variant-supporting tumour reads | ==0 |
| {ORIGIN\|END}_STARTS_STD_DEV | Cluster value for the standard deviation of the supporting breakpoints' starts | <=50.0 |
| {ORIGIN\|END}_MAPQ_MEAN | Cluster value for the mean mapping quality (MAPQ) of the supporting reads | >=40.0 |
| {ORIGIN\|END}_EVENT_SIZE_STD_DEV | Cluster value for the standard deviation of the supporting breakpoints' lengths | <=60.0 |
| CLUSTERED_READS_NORMAL | Number of co-clustered normal reads of any variant type | <=3 |

#### Classify by Parameters File

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
| TUMOUR_READ_SUPPORT | Number of variant-supporting tumour reads | Int |
| NORMAL_READ_SUPPORT | Number of variant-supporting normal reads | Int |
| TUMOUR_AF | Tumour allele-fraction: ratio of tumour-supporting reads to DP (averaged over both edges) | Float |
| NORMAL_AF | Normal allele-fraction: ratio of normal-supporting reads to DP (averaged over both edges) | Float |
| TUMOUR_DP_0 | Total depth/coverage (number of reads) in the tumour at first breakpoint edge | Int |
| TUMOUR_DP_1 | Total depth/coverage (number of reads) in the tumour at second breakpoint edge | Int |
| NORMAL_DP_0 | Total depth/coverage (number of reads) in the normal at first breakpoint edge | Int |
| NORMAL_DP_1 | Total depth/coverage (number of reads) in the normal at second breakpoint edge |  Int |
| {ORIGIN\|END}_STARTS_STD_DEV | Cluster value for the standard deviation of the supporting breakpoints' starts | Float |
| {ORIGIN\|END}_MAPQ_MEAN | Cluster value for the mean mapping quality (MAPQ) of the supporting reads | Float |
| {ORIGIN\|END}_EVENT_SIZE_STD_DEV | Cluster value for the standard deviation of the supporting breakpoints' lengths | Float |
| {ORIGIN\|END}_EVENT_SIZE_MEDIAN | Cluster value for the median of the supporting breakpoints' lengths | Float |
| {ORIGIN\|END}_EVENT_SIZE_MEAN | Cluster value for the mean of the supporting breakpoints' lengths | Float |

#### Classify by Legacy Methods

Alternately, you can use the `--legacy` flag to use filtering and classification methods used in the Beta version of SAVANA. This will output `strict` and `lenient` somatic VCF files which are informed by a decision-tree classifier (strict) and manually plotting data to determine cutoffs (lenient).

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
| --input | VCF file to evaluate |
| --somatic_output | Output VCF path for a file containing only PASS somatic variants |
| --germline_output | Output VCF path for a file containing only PASS germline variants (`--predict_germline` must be specified) |
| --somatic | VCF file containing somatic variants to evaluate PASS somatic variants against |
| --germline | VCF file containing germline variants to evaluate PASS germline variants against |
| --overlap_buffer | If comparing against a --somatic or --germline VCF file, buffer for considering variants to overlap (default=100) |
| --by_support | When comparing to --somatic or --germline VCF, tie-break by read-support |
| --by_distance | When comparing to --somatic or --germline VCF, tie-break by distance (default) |
| --stats | Output filename for statistics on comparison |

### Re-classify Variants

If you'd like to re-classify variants using an alternate method after SAVANA has already been run, you can do so via the `savana classify` sub-command. An example of reclassifying an existing savana output VCF using a custom model would be:
```
savana classify --vcf {raw_sv_breakpoints_vcf} --custom_model {custom_trained_model_pkl} --output {output_classified_vcf}
```

You can also use the `savana classify` sub-command to re-classify using custom parameters (see [Classify by Parameters File](#classify-by-parameters-file)), or legacy methods ([Classify by Legacy Methods](#classify-by-legacy-methods)). See the table below for a full list of arguments:
| Argument | Description |
| -------- | ----------- |
| --input | VCF file to classify |
| --ont | Flag to indicate that the Oxford Nanopore (ONT) trained model should be used to classify variants (default) |
| --pb | Use PacBio filters to classify variants ([see description of filters](classify-for-pacbio)) |
| --predict_germline | Flag to indicate that a model that also predicts germline events should be used (a note that this reduced the accuracy of the somatic calls) |
| --custom_model | Path to custom model pkl file |
| --custom_params | Path to custom paramaters JSON file for filtering |
| --legacy | Flag to use legacy lenient/strict filtering |

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

## Note on SV Types

The `SV_TYPE` field in the `INFO` column of SAVANA only denotes `BND` and `INS` types. We have elected not call `SV_TYPE` beyond these types as it is not definitively possible to do so without copy number information in VCF v4.2 (GRIDSS has a more in-depth explanation for their decision to do this here: https://github.com/PapenfussLab/gridss#why-are-all-calls-bnd).

For now, SAVANA reports the breakend orientation using brackets in the ALT field as described in section 5.4 of [VCFv4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf). We also report this in a `BP_NOTATION` field which can be converted to different nomenclatures as follows:

|Nomenclature  | Deletion-like | Duplication-like  | Head-to-head Inversion  | Tail-to-tail Inversion  |
| ------------------- | ----------------- | --------------------- |  ------------------------------- |  -------------------------- |
| BP_NOTATION  | +- | -+ | ++ | -- |
| Brackets (VCF)  | N[chr:pos[ / ]chr:pos]N | ]chr:pos]N / N[chr:pos[ | N]chr:pos] / N]chr:pos] | [chr:pos[N / [chr:pos[N |
| 5' to 3'  | 3to5 | 5to3 | 3to3 | 5to5|

## Troubleshooting

Please raise a GitHub issue if you encounter issues installing or using SAVANA.

## License

Apache 2.0 License

Copyright (c) 2024 - European Molecular Biology Laboratory (EMBL).
All rights reserved.

## Contacts

Hillary Elrick: helrick@ebi.ac.uk

Carolin Sauer: csauer@ebi.ac.uk

Isidro Cortes-Ciriano: icortes@ebi.ac.uk
