"""
SAVANA strucural variant caller for long-read data - main program
Created: 21/09/2021
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import sys
import os
import argparse

from time import time
from multiprocessing import cpu_count

import pysam

import savana.run as run
import savana.evaluate as evaluate
import savana.train as train
import savana.classify as classify
import savana.allele_counter as allele_counter
import savana.bin_generator as bin_generator
import savana.smooth as smooth
import savana.read_counter as read_counter
import savana.segment as segment
import savana.fit_absolute as fit_absolute
import savana.helper as helper

logo = """
███████  █████  ██    ██  █████  ███    ██  █████
██      ██   ██ ██    ██ ██   ██ ████   ██ ██   ██
███████ ███████ ██    ██ ███████ ██ ██  ██ ███████
     ██ ██   ██  ██  ██  ██   ██ ██  ██ ██ ██   ██
███████ ██   ██   ████   ██   ██ ██   ████ ██   ██
"""

def savana_run(args):
    """ identify raw breakpoints """
    args.tumour_only = False
    if not args.sample:
        # set sample name to default if req.
        args.sample = os.path.splitext(os.path.basename(args.tumour))[0]
    print(f'Running as sample {args.sample}')
    outdir = helper.check_outdir(args.outdir, args.overwrite, illegal='.vcf')
    # set number of threads to cpu count if none set
    if not args.threads:
        args.threads = cpu_count()
    # check if files are bam or cram (must have indices)
    if args.tumour.endswith('bam') and args.normal.endswith('bam'):
        args.is_cram = False
        aln_files = {
            'tumour': pysam.AlignmentFile(args.tumour, "rb"),
            'normal': pysam.AlignmentFile(args.normal, "rb")
        }
    elif args.tumour.endswith('cram') and args.normal.endswith('cram'):
        args.is_cram = True
        aln_files = {
            'tumour': pysam.AlignmentFile(args.tumour, "rc"),
            'normal': pysam.AlignmentFile(args.normal, "rc")
        }
    else:
        sys.exit('Unrecognized file extension. Tumour and normal files must be BAM/CRAM')

    # confirm ref and ref fasta index exist
    if not os.path.exists(args.ref):
        sys.exit(f'Provided reference: "{args.ref}" does not exist. Please provide full path')
    elif args.ref_index and not os.path.exists(args.ref_index):
        sys.exit(f'Provided reference fasta index: "{args.ref_index}" does not exist. Please provide full path')
    elif not os.path.exists(f'{args.ref}.fai'):
        sys.exit(f'Default reference fasta index: "{args.ref}.fai" does not exist. Please provide full path')
    else:
        args.ref_index = f'{args.ref}.fai' if not args.ref_index else args.ref_index
        print(f'Found {args.ref_index} to use as reference fasta index')
    # initialize timing
    checkpoints = [time()]
    time_str = []
    # run SAVANA processes
    checkpoints, time_str = run.spawn_processes(args, aln_files, checkpoints, time_str, outdir)
    # finish timing
    helper.time_function("Total time to call raw variants", checkpoints, time_str, final=True)

def savana_classify(args):
    """ main function for savana classify """
    # initialize timing
    checkpoints = [time()]
    time_str = []
    # perform logic checks and set defaults
    if not args.ont and not args.pb and not args.custom_model and not args.custom_params and not args.legacy:
        print(f'Using ONT model to classify variants')
        args.ont = True
    if args.predict_germline:
        if args.custom_model:
            sys.exit(f'The `predict_germline` flag cannot be used in conjunction with a custom model')
        elif args.custom_params:
            sys.exit(f'The `predict_germline` flag cannot be used in conjunction with custom parameters. Please include germline criteria in JSON input')
        elif args.legacy:
            sys.exit(f'The `predict_germline` flag cannot be used in conjunction with legacy filtering')

    if args.legacy:
        classify.classify_legacy(args, checkpoints, time_str)
    elif args.custom_params:
        classify.classify_by_params(args, checkpoints, time_str)
    elif args.custom_model:
        # using a model - check whether has been un-tarred
        import tarfile
        if tarfile.is_tarfile(args.custom_model):
            # see if untarred file exists without tar.gz extension
            untar_model = '.'.join(args.custom_model.split('.')[:-2])
            if os.path.isfile(untar_model):
                print(f'Using untarred model path at {untar_model}')
                args.custom_model = untar_model
            else:
                print(f'Will untar {args.custom_model} to {untar_model}')
                tar = tarfile.open(args.custom_model, "r:gz")
                model_dir = os.path.dirname(args.custom_model)
                tar.extractall(model_dir)
                tar.close()
                args.custom_model = untar_model
        elif not os.path.isfile(args.custom_model):
            print(f'Unable to locate model at {args.custom_model} - please check installation')
            return
        classify.classify_by_model(args, checkpoints, time_str)
    elif args.pb:
        classify.classify_pacbio(args, checkpoints, time_str)
    elif args.ont:
        # using a model - perform logic to determine which one
        if args.tumour_only:
            model_base = 'ont-somatic-to'
        elif args.predict_germline:
            model_base = 'ont-germline'
        else:
            model_base = 'ont-somatic'
        # check whether the model has been un-tarred
        models_dir = os.path.join(os.path.dirname(__file__), 'models')
        model_path = os.path.join(models_dir, model_base)
        model_pkl = model_path+'.pkl'
        model_tar = model_path+'.tar.gz'
        if os.path.isfile(model_pkl):
            print(f'Using untarred {model_pkl}')
            args.custom_model = model_pkl
        elif os.path.isfile(model_tar):
            import tarfile
            print(f'First time using model - will untar {model_tar}')
            tar = tarfile.open(model_tar, "r:gz")
            tar.extractall(models_dir)
            tar.close()
            args.custom_model = model_pkl
        else:
            print(f'Unable to locate model at {model_path}.* - please check installation')
            return
        classify.classify_by_model(args, checkpoints, time_str)

    if args.cna_rescue:
        classify.rescue_cna(args, checkpoints, time_str)

    # finish timing
    helper.time_function("Total time to classify variants", checkpoints, time_str, final=True)

def savana_evaluate(args):
    """ main function for savana evaluate """
    # check input VCFs
    vcf_string = ''
    if not os.path.exists(args.input):
        sys.exit(f'Provided input vcf: "{args.input}" does not exist. Please provide full path')
    if not os.path.exists(args.somatic):
        sys.exit(f'Provided somatic VCF: "{args.somatic}" does not exist. Please provide full path')
    else:
        vcf_string += f'somatic vcf: "{args.somatic}"'
    if args.germline:
        if not os.path.exists(args.germline):
            sys.exit(f'Provided germline VCF: "{args.germline}" does not exist. Please provide full path')
        else:
            vcf_string += f' and germline vcf: "{args.germline}"'

    # initialize timing
    checkpoints = [time()]
    time_str = []
    # perform validation
    evaluate.evaluate_vcf(args, checkpoints, time_str)
    # finish timing
    helper.time_function("Total time to evaluate variants", checkpoints, time_str, final=True)

def savana_train(args):
    """ main function for savana train """
    outdir = helper.check_outdir(args.outdir, args.overwrite, illegal='.pkl')
    data_matrix = None
    if args.vcfs:
        # read in and create matrix from VCF files
        data_matrix = train.create_dataframe(args)
    elif args.load_matrix:
        # load data matrix from pickle file
        data_matrix = train.load_matrix(args)
    features, target = train.prepare_data(data_matrix, germline_class=args.germline_class, tumour_only=args.tumour_only)
    print('Done preparing data.')
    classifier = train.cross_conformal_classifier(features, target, outdir, args.test_split, args.threads)
    train.save_model(args, classifier, outdir)

def savana_cna(args, as_workflow=False):
    """ main function for copy number analysis """
    if not as_workflow:
        # not being run under savana umbrella - do these checks/setup
        if not args.sample:
            # set sample name to default
            args.sample = os.path.splitext(os.path.basename(args.tumour))[0]
        outdir = helper.check_outdir(args.outdir, args.overwrite)
    else:
        outdir = args.outdir
    # define tmpdir
    tmpdir = helper.check_tmpdir(args.tmpdir, outdir, args.overwrite)
    # initialize timing
    checkpoints = [time()]
    time_str = []
    # cna threads
    if not args.cna_threads:
        args.cna_threads = cpu_count()
    # if args.cna_threads:
    #     args.threads = args.cna_threads
    # first do allele counting
    if not args.allele_counts_het_snps:
        if args.snp_vcf:
            args.g1000_vcf = None
            allele_counts_bed_path = allele_counter.perform_allele_counting(outdir, args.sample, args.chromosomes,  args.ref, args.snp_vcf, args.g1000_vcf, args.tumour, args.ac_window, args.allele_mapq, args.allele_min_reads, tmpdir, args.cna_threads)
            helper.time_function("Counted alleles heterozygous SNPs", checkpoints, time_str)
        if args.g1000_vcf:
            args.snp_vcf = None
            allele_counts_bed_path = allele_counter.perform_allele_counting(outdir, args.sample, args.chromosomes,  args.ref, args.snp_vcf, args.g1000_vcf, args.tumour, args.ac_window, args.allele_mapq, args.allele_min_reads, tmpdir, args.cna_threads)
            helper.time_function("Counted alleles of 1000g SNPs", checkpoints, time_str)
    else:
        allele_counts_bed_path = args.allele_counts_het_snps
    # now generate bins
    bin_annotations_path = bin_generator.generate_bins(outdir, args.sample, args.ref, args.chromosomes, args.cn_binsize, args.blacklist, args.breakpoints, args.cna_threads)
    helper.time_function("Binned reference genome", checkpoints, time_str)
    # perform the read counting
    read_counts_path = read_counter.count_reads(outdir, args.tumour, args.normal, args.sample, bin_annotations_path, args.readcount_mapq, args.blacklisting, args.bl_threshold, args.bases_filter, args.bases_threshold, args.cna_threads)
    helper.time_function("Performed read counting", checkpoints, time_str)
    # smooth the copy number data
    smoothened_cn_path = smooth.smooth_copy_number(outdir, read_counts_path, args.smoothing_level, args.trim)
    helper.time_function("Performed smoothing", checkpoints, time_str)
    # segment the copy number data
    log2r_cn_path = segment.segment_copy_number(outdir, smoothened_cn_path, args.min_segment_size, args.shuffles, args.p_seg, args.p_val, args.quantile, args.cna_threads)
    helper.time_function("Performed CBS", checkpoints, time_str)
    # fit absolute copy number
    fit_absolute.fit_absolute_cn(outdir, log2r_cn_path, allele_counts_bed_path, args.sample,
        args.min_ploidy, args.max_ploidy, args.ploidy_step, args.min_cellularity, args.max_cellularity, args.cellularity_step, args.cellularity_buffer, args.overrule_cellularity,
        args.distance_function, args.distance_filter_scale_factor, args.distance_precision,
        args.max_proportion_zero, args.min_proportion_close_to_whole_number, args.max_distance_from_whole_number, args.main_cn_step_change,
        args.min_block_size, args.min_block_length, args.cna_threads)
    helper.time_function("Fit absolute copy number", checkpoints, time_str)
    # cleanup tmpdir
    helper.clean_tmpdir(args.tmpdir, outdir)
    helper.time_function("Total time to perform copy number calling", checkpoints, time_str, final=True)

def savana_tumour_only(args):
    """ identify breakpoints and classify them with model - without matched normal """
    args.tumour_only = True
    # germline not possible - manually declare to prevent nonsensical argument to TO parser
    args.predict_germline, args.germline_output, args.germline = False, False, False

    if not args.sample:
        # set sample name to default if req.
        args.sample = os.path.splitext(os.path.basename(args.tumour))[0]
    print(f'Running as sample {args.sample}')
    outdir = helper.check_outdir(args.outdir, args.overwrite, illegal='.vcf')
    # set number of threads to cpu count if none set
    if not args.threads:
        args.threads = cpu_count()
    # check if files are bam or cram (must have indices)
    if args.tumour.endswith('bam'):
        args.is_cram = False
        aln_files = {
            'tumour': pysam.AlignmentFile(args.tumour, "rb")
        }
    elif args.tumour.endswith('cram'):
        args.is_cram = True
        aln_files = {
            'tumour': pysam.AlignmentFile(args.tumour, "rc")
        }
    else:
        sys.exit('Unrecognized file extension. Tumour files must be BAM/CRAM')

    # confirm ref and ref fasta index exist
    if not os.path.exists(args.ref):
        sys.exit(f'Provided reference: "{args.ref}" does not exist. Please provide full path')
    elif args.ref_index and not os.path.exists(args.ref_index):
        sys.exit(f'Provided reference fasta index: "{args.ref_index}" does not exist. Please provide full path')
    elif not os.path.exists(f'{args.ref}.fai'):
        sys.exit(f'Default reference fasta index: "{args.ref}.fai" does not exist. Please provide full path')
    else:
        args.ref_index = f'{args.ref}.fai' if not args.ref_index else args.ref_index
        print(f'Found {args.ref_index} to use as reference fasta index')
    # initialize timing
    checkpoints = [time()]
    time_str = []
    # run SAVANA processes
    checkpoints, time_str = run.spawn_processes(args, aln_files, checkpoints, time_str, outdir)
    # finish timing
    helper.time_function("Total time to call raw variants", checkpoints, time_str, final=True)

    # set inputs and outputs
    args.vcf = os.path.join(args.outdir,f'{args.sample}.sv_breakpoints.vcf') # previous step's output
    args.output = os.path.join(args.outdir,f'{args.sample}.classified.vcf')
    if not args.somatic_output:
        args.somatic_output = os.path.join(args.outdir,f'{args.sample}.classified.somatic.vcf')

    savana_classify(args)

    # perform evaluation against provided VCFs, if any
    if args.somatic:
        # evaluate against somatic "truthset" VCF
        args.input = args.somatic_output if args.somatic_output else os.path.join(args.outdir, f'{args.sample}.sv_breakpoints.vcf')
        args.output = os.path.join(args.outdir, f'{args.sample}.somatic.labelled.vcf')
        args.stats = os.path.join(args.outdir,f'{args.sample}.somatic.evaluation.stats')
        savana_evaluate(args)

    args.breakpoints = args.somatic_output

    if args.snp_vcf or args.g1000_vcf or args.allele_counts_het_snps:
        args.normal = None
        savana_cna(args, True)

    return


def savana_main(args):
    """ default workflow for savana: savana_run, savana_classify, savana_evaluate, savana_cna """

    # call raw breakpoints
    savana_run(args)

    # set inputs and outputs
    args.vcf = os.path.join(args.outdir,f'{args.sample}.sv_breakpoints.vcf') # previous step's output
    args.output = os.path.join(args.outdir,f'{args.sample}.classified.vcf')
    if not args.somatic_output:
        args.somatic_output = os.path.join(args.outdir,f'{args.sample}.classified.somatic.vcf')
    if args.germline and not args.germline_output:
        args.germline_output = os.path.join(args.outdir,f'{args.sample}.classified.germline.vcf')

    # perform classification by selected method
    savana_classify(args)

    # perform evaluation against provided VCFs, if any
    if args.somatic:
        # evaluate against somatic "truthset" VCF
        args.input = args.somatic_output if args.somatic_output else os.path.join(args.outdir, f'{args.sample}.sv_breakpoints.vcf')
        args.output = os.path.join(args.outdir, f'{args.sample}.somatic.labelled.vcf')
        args.stats = os.path.join(args.outdir,f'{args.sample}.somatic.evaluation.stats')
        savana_evaluate(args)
    if args.germline:
        # evaluate against germline "truthset" VCF
        args.input = args.germline_output if args.germline else os.path.join(args.outdir, f'{args.sample}.sv_breakpoints.vcf')
        args.output = os.path.join(args.outdir, f'{args.sample}.germline.labelled.vcf')
        args.stats = os.path.join(args.outdir,f'{args.sample}.germline.evaluation.stats')
        savana_evaluate(args)

    args.breakpoints = args.somatic_output

    if args.snp_vcf or args.g1000_vcf or args.allele_counts_het_snps:
        savana_cna(args, True)

    return

def parse_args(args):
    """ parse arguments - separated into subcommands """
    global_parser = argparse.ArgumentParser(prog="savana", description="SAVANA - somatic SV caller")
    global_parser.add_argument('--version', action='version', version=f'SAVANA {helper.__version__}')
    subparsers = global_parser.add_subparsers(title="subcommands", help='SAVANA sub-commands', dest='command')
    subparsers.required = False

    # savana run
    run_parser = subparsers.add_parser("run", help="identify and cluster breakpoints - output raw variants without classification")
    run_parser.add_argument('-t','--tumour', nargs='?', type=str, required=True, help='Tumour BAM file (must have index)')
    run_parser.add_argument('-n', '--normal', nargs='?', type=str, required=True, help='Normal BAM file (must have index)')
    run_parser.add_argument('-r', '--ref', nargs='?', type=str, required=True, help='Full path to reference genome')
    run_parser.add_argument('--ref_index', nargs='?', type=str, required=False, help='Full path to reference genome fasta index (ref path + ".fai" by default)')
    run_parser.add_argument('--contigs', nargs='?', type=str, help="Contigs/chromosomes to consider. See example at example/contigs.chr.hg38.txt (optional, default=All)")
    run_parser.add_argument('--length', nargs='?', type=int, default=30, help='Minimum length SV to consider (default=30)')
    run_parser.add_argument('--mapq', nargs='?', type=int, default=5, help='Minimum MAPQ to consider a read mapped (default=5)')
    run_parser.add_argument('--buffer', nargs='?', type=int, default=10, help='Buffer when clustering adjacent potential breakpoints, excepting insertions (default=10)')
    run_parser.add_argument('--insertion_buffer', nargs='?', type=int, default=250, help='Buffer when clustering adjacent potential insertion breakpoints (default=250)')
    run_parser.add_argument('--end_buffer', nargs='?', type=int, default=50, help='Buffer when clustering alternate edge of potential breakpoints, excepting insertions (default=50)')
    run_parser.add_argument('--threads', nargs='?', type=int, const=0, help='Number of threads to use (default=max)')
    run_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
    run_parser.add_argument('--overwrite', action='store_true', help='Use this flag to write to output directory even if files are present')
    run_parser.add_argument('-s','--sample', nargs='?', type=str, help="Name to prepend to output files (default=tumour BAM filename without extension)")
    run_parser.add_argument('--keep_inv_artefact', action='store_true', help="Do not remove breakpoints with foldback-inversion artefact pattern")
    run_parser.add_argument('--single_bnd', action='store_true', help='Report single breakend variants in addition to standard types (default=False)')
    run_parser.add_argument('--single_bnd_min_length', nargs='?', type=int, default=1000, help='Minimum length of single breakend to consider (default=100)')
    run_parser.add_argument('--single_bnd_max_mapq', nargs='?', type=int, default=5, help='Convert supplementary alignments below this threshold to single breakend (default=5, must not exceed --mapq argument)')
    run_parser.add_argument('--debug', action='store_true', help='output extra debugging info and files')
    run_parser.add_argument('--chunksize', nargs='?', type=int, default=100000000, help='Chunksize to use when splitting genome for parallel analysis - used to optimise memory (default=1000000)')
    run_parser.add_argument('--coverage_binsize', nargs='?', type=int, default=5, help='Length used for coverage bins (default=5)')
    run_parser.add_argument('--min_support', nargs='?', type=int, default=3, required=False, help='Minimum supporting reads for a PASS variant (default=3)')
    run_parser.add_argument('--min_af', nargs='?', type=helper.float_range(0.0, 1.0), default=0.01, required=False, help='Minimum allele-fraction for a PASS variant (default=0.01)')
    run_parser.set_defaults(func=savana_run)

    # savana classify
    classify_parser = subparsers.add_parser("classify", help="classify VCF using model")
    classify_parser.add_argument('--vcf', nargs='?', type=str, required=True, help='VCF file to classify')
    classify_parser.add_argument('--min_support', nargs='?', type=int, default=3, required=False, help='Minimum supporting reads for a PASS variant')
    classify_parser.add_argument('--min_af', nargs='?', type=helper.float_range(0.0, 1.0), default=0.01, required=False, help='Minimum allele-fraction for a PASS variant')
    classify_parser.add_argument('--cna_rescue', nargs='?', type=str, required=False, help='Copy number abberation output file for this sample (used to rescue variants)')
    classify_parser.add_argument('--cna_rescue_distance', nargs='?', type=int, default=50, required=False, help='Maximum distance from a copy number abberation for a variant to be rescued')
    group = classify_parser.add_mutually_exclusive_group()
    group.add_argument('--ont', action='store_true', help='Use the Oxford Nanopore (ONT) trained model to classify variants (default)')
    group.add_argument('--pb', action='store_true', help='Use PacBio thresholds to classify variants')
    # whether to use a germline-trained model
    classify_parser.add_argument('--predict_germline', action='store_true', help='Also predict germline events (reduces accuracy of somatic calls for model)')
    group.add_argument('--custom_model', nargs='?', type=str, required=False, help='Pickle file of custom machine-learning model')
    group.add_argument('--custom_params', nargs='?', type=str, required=False, help='JSON file of custom filtering parameters')
    group.add_argument('--legacy', action='store_true', help='Legacy lenient/strict filtering')
    classify_parser.add_argument('-to','--tumour_only', action='store_true', help='Classifying tumour-only data')
    classify_parser.add_argument('--output', nargs='?', type=str, required=True, help='Output VCF with PASS columns and CLASS added to INFO')
    classify_parser.add_argument('--somatic_output', nargs='?', type=str, required=False, help='Output VCF containing only PASS somatic variants')
    classify_parser.add_argument('--germline_output', nargs='?', type=str, required=False, help='Output VCF containing only PASS germline variants')
    classify_parser.add_argument('--confidence', nargs='?', type=helper.float_range(0.0, 1.0), default=None, help='Confidence level for mondrian conformal prediction - suggested range (0.70-0.99) (not used by default)')
    classify_parser.add_argument('--threads', nargs='?', type=int, default=16, const=0, help='Number of threads to use')
    classify_parser.set_defaults(func=savana_classify)

    # savana evaluate
    evaluate_parser = subparsers.add_parser("evaluate", help="label SAVANA VCF with somatic/germline/missing given VCF(s) to compare against")
    evaluate_parser.add_argument('--input', nargs='?', type=str, required=True, help='VCF file to evaluate')
    evaluate_parser.add_argument('--somatic', nargs='?', type=str, required=True, help='Somatic VCF file to evaluate against')
    evaluate_parser.add_argument('--germline', nargs='?', type=str, required=False, help='Germline VCF file to evaluate against (optional)')
    evaluate_parser.add_argument('--overlap_buffer', nargs='?', type=int, default=100, help='Buffer for considering an overlap (default=100)')
    evaluate_parser.add_argument('--output', nargs='?', type=str, required=True, help='Output VCF with LABEL added to INFO')
    evaluate_parser.add_argument('--stats', nargs='?', type=str, required=False, help='Output file for statistics on comparison if desired')
    evaluate_parser.add_argument('--curate', action='store_true', default=False, help='Attempt to reduce false labels for training (allow label to be used twice)')
    evaluate_parser.add_argument('--qual_filter', nargs='?', type=int, default=None, help='Impose a quality filter on comparator variants (default=None)')
    group = evaluate_parser.add_mutually_exclusive_group()
    group.add_argument('--by_support', action='store_true', help='Comparison method: tie-break by read support')
    group.add_argument('--by_distance', action='store_true', default=True, help='Comparison method: tie-break by min. distance (default)')
    evaluate_parser.set_defaults(func=savana_evaluate)

    # savana train
    train_parser = subparsers.add_parser("train", help="train model on folder of input VCFs")
    group = train_parser.add_mutually_exclusive_group()
    group.add_argument('--vcfs', nargs='?', type=str, required=False, help='Folder of labelled VCF files to read in')
    train_parser.add_argument('--recursive', action='store_true', help='Search recursively through input folder for input VCFs (default only one-level deep)')
    group.add_argument('--load_matrix', nargs='?', type=str, required=False, help='Pre-loaded pickle file of VCFs')
    train_parser.add_argument('--save_matrix', nargs='?', type=str, required=False, help='Output pickle file for data matrix of VCFs')
    train_parser.add_argument('--test_split', nargs='?', type=float, default=0.2, help='Fraction of data to use for test (default=0.2)')
    train_parser.add_argument('--germline_class', action='store_true', help='Train the model to predict germline and somatic variants (GERMLINE label must be present)')
    train_parser.add_argument('-to','--tumour_only', action='store_true', help='Training a model on tumour-only data')
    train_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
    train_parser.add_argument('--overwrite', action='store_true', help='Use this flag to write to output directory even if files are present')
    train_parser.add_argument('--threads', nargs='?', type=int, default=16, const=0, help='Number of threads to use')
    train_parser.set_defaults(func=savana_train)

    # savana cna
    cna_parser = subparsers.add_parser("cna", help="run copy number")
    cna_parser.add_argument('-t','--tumour', nargs='?', type=str, required=True, help='Tumour BAM file (must have index)')
    cna_parser.add_argument('-n', '--normal', nargs='?', type=str, required=False, help='Normal BAM file (must have index)')
    cna_parser.add_argument('--ref', nargs='?', type=str, required=True, help='Full path to reference genome')
    cna_parser.add_argument('--sample', nargs='?', type=str, help="Name to prepend to output files (default=tumour BAM filename without extension)")
    cna_parser.add_argument('--cna_threads', nargs='?', type=int, const=0, help='Number of threads to use for CNA (default=max)')
    cna_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
    cna_parser.add_argument('--overwrite', action='store_true', help='Use this flag to write to output directory even if files are present')
    cna_parser.add_argument('--tmpdir', nargs='?', required=False, default='tmp', help='Temp directory for allele counting temp files (defaults to outdir)')
    allele_group = cna_parser.add_mutually_exclusive_group()
    allele_group.add_argument('-v', '--snp_vcf', type=str, help='Path to matched normal SNP vcf file to extract heterozygous SNPs for allele counting.', required=False)
    allele_group.add_argument('-vg', '--g1000_vcf', type=str, choices={"1000g_hg38", "1000g_hg19", "1000g_t2t"}, help='Use 1000g biallelic vcf file for allele counting instead of SNP vcf from matched normal. Specify which genome version to use.', required=False)
    allele_group.add_argument('-ac', '--allele_counts_het_snps', type=str, help='Path to allele counts of heterozygous SNPs', required=False)
    cna_parser.add_argument('-q', '--allele_mapq', type=int,  default=5, help='Mapping quality threshold for reads to be included in the allele counting (default = 5)', required=False)
    cna_parser.add_argument('-mr', '--allele_min_reads', type=int,  default=10, help='Minimum number of reads required per het SNP site for allele counting (default = 10)', required=False)
    cna_parser.add_argument('-acw','--ac_window', type=int, default=1200000, help='Window size for allele counting to parallelise and use for purity estimation (default = 1200000; this should be >=500000)', required=False)
    cna_parser.add_argument('-w', '--cn_binsize', type=int, default=10, help='Bin window size in kbp', required=False)
    cna_parser.add_argument('-b', '--blacklist', type=str, help='Path to the blacklist file', required=False)
    cna_parser.add_argument('-bp', '--breakpoints', type=str, help='Path to SAVANA VCF file to incorporate savana breakpoints into copy number analysis', required=False)
    cna_parser.add_argument('-c', '--chromosomes', nargs='+', default='all', help='Chromosomes to analyse. To run on all chromosomes, leave unspecified (default). To run on a subset of chromosomes only, specify the chromosome numbers separated by spaces. For x and y chromosomes, use 23 and 24, respectively.  E.g. use "-c 1 4 23 24" to run chromosomes 1, 4, X and Y', required=False)
    cna_parser.add_argument('-rq', '--readcount_mapq', type=int,  default=5, help='Mapping quality threshold for reads to be included in the read counting (default = 5)', required=False)
    cna_parser.add_argument('--no_blacklist', dest='blacklisting', action='store_false')
    cna_parser.set_defaults(blacklisting=True)
    cna_parser.add_argument('-blt', '--bl_threshold', type=int,  default='5', help='Percentage overlap between bin and blacklist threshold to tolerate for read counting (default = 5). Please specify percentage threshold as integer, e.g. "-t 5". Set "-t 0" if no overlap with blacklist is to be tolerated', required=False)
    cna_parser.add_argument('--no_basesfilter', dest='bases_filter', action='store_false')
    cna_parser.set_defaults(bases_filter=True)
    cna_parser.add_argument('-bt', '--bases_threshold', type=int,  default='75', help='Percentage of known bases per bin required for read counting (default = 75). Please specify percentage threshold as integer, e.g. "-bt 95" ', required=False)
    cna_parser.add_argument('-sl', '--smoothing_level', type=int, default='10', help='Size of neighbourhood for smoothing.', required=False)
    cna_parser.add_argument('-tr', '--trim', type=float, default='0.025', help='Trimming percentage to be used.', required=False)
    cna_parser.add_argument('-ms', '--min_segment_size', type=int,  default=5, help='Minimum size for a segement to be considered a segment (default = 5).', required=False)
    cna_parser.add_argument('-sf', '--shuffles', type=int,  default=1000, help='Number of permutations (shuffles) to be performed during CBS (default = 1000).', required=False)
    cna_parser.add_argument('-ps', '--p_seg', type=float,  default=0.05, help='p-value used to test segmentation statistic for a given interval during CBS using (shuffles) number of permutations (default = 0.05).', required=False)
    cna_parser.add_argument('-pv', '--p_val', type=float,  default=0.01, help='p-value used to test validity of candidate segments from CBS using (shuffles) number of permutations (default = 0.01).', required=False)
    cna_parser.add_argument('-qt', '--quantile', type=float,  default=0.2, help='Quantile of changepoint (absolute median differences across all segments) used to estimate threshold for segment merging (default = 0.2; set to 0 to avoid segment merging).', required=False)
    cna_parser.add_argument('--min_ploidy', type=float, default=1.5, help='Minimum ploidy to be considered for copy number fitting.', required=False)
    cna_parser.add_argument('--max_ploidy', type=float, default=5, help='Maximum ploidy to be considered for copy number fitting.', required=False)
    cna_parser.add_argument('--ploidy_step', type=float, default=0.01, help='Ploidy step size for grid search space used during for copy number fitting.', required=False)
    cna_parser.add_argument('--min_cellularity', type=float, default=0.2, help='Minimum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided, this is estimated during copy number fitting. Alternatively, a purity value can be provided if the purity of the sample is already known.', required=False)
    cna_parser.add_argument('--max_cellularity', type=float, default=1, help='Maximum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided, this is estimated during copy number fitting. Alternatively, a purity value can be provided if the purity of the sample is already known.', required=False)
    cna_parser.add_argument('--cellularity_step', type=float, default=0.01, help='Cellularity step size for grid search space used during for copy number fitting.', required=False)
    cna_parser.add_argument('--cellularity_buffer', type=float, default=0.1, help='Cellularity buffer to define purity grid search space during copy number fitting (default = 0.1).', required=False)
    cna_parser.add_argument('--overrule_cellularity', type=float, default=None, help='Set to sample`s purity if known. This value will overrule the cellularity estimated using hetSNP allele counts (not used by default).', required=False)   
    cna_parser.add_argument('--distance_function', type=str, default='RMSD', help='Distance function to be used for copy number fitting.', choices=['RMSD', 'MAD'], required=False)
    cna_parser.add_argument('--distance_filter_scale_factor', type=float, default=1.25, help='Distance filter scale factor to only include solutions with distances < scale factor * min(distance).', required=False)
    cna_parser.add_argument('--distance_precision', type=int, default=3, help='Number of digits to round distance functions to', required=False)
    cna_parser.add_argument('--max_proportion_zero', type=float, default=0.1, help='Maximum proportion of fitted copy numbers to be tolerated in the zero or negative copy number state.', required=False)
    cna_parser.add_argument('--min_proportion_close_to_whole_number', type=float, default=0.5, help='Minimum proportion of fitted copy numbers sufficiently close to whole number to be tolerated for a given fit.', required=False)
    cna_parser.add_argument('--max_distance_from_whole_number', type=float, default=0.25, help='Distance from whole number for fitted value to be considered sufficiently close to nearest copy number integer.', required=False)
    cna_parser.add_argument('--main_cn_step_change', type=int, default=1, help='Max main copy number step change across genome to be considered for a given solution.', required=False)
    cna_parser.add_argument('--min_block_size', type=int, default=10, help='Minimum size (number of SNPs) for a genomic block to be considered for purity estimation.', required=False)
    cna_parser.add_argument('--min_block_length', type=int, default=500000, help='Minimum length (bps) for a genomic block to be considered for purity estimation.', required=False)
    cna_parser.set_defaults(func=savana_cna)

    # tumour only parser
    to_parser = subparsers.add_parser("to", help="identify somatic SVs without normal sample")
    to_parser.add_argument('-t','--tumour', nargs='?', type=str, required=True, help='Tumour BAM file (must have index)')
    to_parser.add_argument('-r', '--ref', nargs='?', type=str, required=True, help='Full path to reference genome')
    to_parser.add_argument('--ref_index', nargs='?', type=str, required=False, help='Full path to reference genome fasta index (ref path + ".fai" by default)')
    to_parser.add_argument('--contigs', nargs='?', type=str, help="Contigs/chromosomes to consider. See example at example/contigs.chr.hg38.txt (optional, default=All)")
    to_parser.add_argument('--length', nargs='?', type=int, default=30, help='Minimum length SV to consider (default=30)')
    to_parser.add_argument('--mapq', nargs='?', type=int, default=5, help='Minimum MAPQ to consider a read mapped (default=5)')
    to_parser.add_argument('--buffer', nargs='?', type=int, default=10, help='Buffer when clustering adjacent potential breakpoints, excepting insertions (default=10)')
    to_parser.add_argument('--insertion_buffer', nargs='?', type=int, default=250, help='Buffer when clustering adjacent potential insertion breakpoints (default=250)')
    to_parser.add_argument('--end_buffer', nargs='?', type=int, default=50, help='Buffer when clustering alternate edge of potential breakpoints, excepting insertions (default=50)')
    to_parser.add_argument('--threads', nargs='?', type=int, const=0, help='Number of threads to use (default=max)')
    to_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
    to_parser.add_argument('--overwrite', action='store_true', help='Use this flag to write to output directory even if files are present')
    to_parser.add_argument('-s','--sample', nargs='?', type=str, help="Name to prepend to output files (default=tumour BAM filename without extension)")
    to_parser.add_argument('--keep_inv_artefact', action='store_true', help="Do not remove breakpoints with foldback-inversion artefact pattern")
    to_parser.add_argument('--single_bnd', action='store_true', help='Report single breakend variants in addition to standard types (default=False)')
    to_parser.add_argument('--single_bnd_min_length', nargs='?', type=int, default=1000, help='Minimum length of single breakend to consider (default=100)')
    to_parser.add_argument('--single_bnd_max_mapq', nargs='?', type=int, default=5, help='Convert supplementary alignments below this threshold to single breakend (default=5, must not exceed --mapq argument)')
    to_parser.add_argument('--debug', action='store_true', help='Output extra debugging info and files')
    to_parser.add_argument('--chunksize', nargs='?', type=int, default=100000000, help='Chunksize to use when splitting genome for parallel analysis - used to optimise memory (default=1000000)')
    to_parser.add_argument('--coverage_binsize', nargs='?', type=int, default=5, help='Length used for coverage bins (default=5)')
    to_parser.add_argument('--min_support', nargs='?', type=int, default=3, required=False, help='Minimum supporting reads for a PASS variant (default=3)')
    to_parser.add_argument('--min_af', nargs='?', type=helper.float_range(0.0, 1.0), default=0.01, required=False, help='Minimum allele-fraction for a PASS variant (default=0.01)')
    # classify args
    to_parser.add_argument('--cna_rescue', nargs='?', type=str, required=False, help='Copy number abberation output file for this sample (used to rescue variants)')
    to_parser.add_argument('--cna_rescue_distance', nargs='?', type=int, default=50, required=False, help='Maximum distance from a copy number abberation for a variant to be rescued')
    classify_group = to_parser.add_mutually_exclusive_group()
    classify_group.add_argument('--ont', action='store_true', help='Use the Oxford Nanopore (ONT) trained model to classify variants (default)')
    classify_group.add_argument('--pb', action='store_true', help='Use PacBio thresholds to classify variants')
    classify_group.add_argument('--custom_model', nargs='?', type=str, required=False, help='Pickle file of custom machine-learning model')
    classify_group.add_argument('--custom_params', nargs='?', type=str, required=False, help='JSON file of custom filtering parameters')
    classify_group.add_argument('--legacy', action='store_true', help='Use legacy lenient/strict filtering')
    to_parser.add_argument('--somatic_output', nargs='?', type=str, required=False, help='Output VCF with only PASS somatic variants')
    to_parser.add_argument('--confidence', nargs='?', type=helper.float_range(0.0, 1.0), default=None, help='Confidence level for mondrian conformal prediction - suggested range (0.70-0.99) (not used by default)')
    # evaluate args
    to_parser.add_argument('--somatic', nargs='?', type=str, required=False, help='Somatic VCF file to evaluate against')
    to_parser.add_argument('--overlap_buffer', nargs='?', type=int, default=100, required=False, help='Buffer for considering an overlap (default=100)')
    to_parser.add_argument('--curate', action='store_true', default=False, help='Attempt to reduce false labels for training (allow label to be used twice)')
    to_parser.add_argument('--qual_filter', nargs='?', type=int, default=None, help='Impose a quality filter on comparator variants (default=None)')
    evaluate_group = to_parser.add_mutually_exclusive_group()
    evaluate_group.add_argument('--by_support', action='store_true', help='Comparison method: tie-break by read support')
    evaluate_group.add_argument('--by_distance', action='store_true', default=True, help='Comparison method: tie-break by min. distance (default)')
    # cna arguments for tumour only mode
    to_parser.add_argument('--cna_threads', nargs='?', type=int, const=0, help='Number of threads to use for CNA (default=max)')
    to_parser.add_argument('--tmpdir', nargs='?', required=False, default='tmp', help='Temp directory for allele counting temp files (defaults to outdir)')
    allele_group_to = to_parser.add_mutually_exclusive_group()
    allele_group_to.add_argument('-v', '--snp_vcf', type=str, help='Path to SNP vcf file to extract heterozygous SNPs for allele counting.', required=False)
    allele_group_to.add_argument('-vg', '--g1000_vcf', type=str, choices={"1000g_hg38", "1000g_hg19", "1000g_t2t"}, help='Use 1000g biallelic vcf file for allele counting instead of SNP vcf from matched normal. Specify which genome version to use.', required=False)
    allele_group_to.add_argument('-ac', '--allele_counts_het_snps', type=str, help='Path to allele counts of heterozygous SNPs', required=False)
    to_parser.add_argument('-q', '--allele_mapq', type=int,  default=5, help='Mapping quality threshold for reads to be included in the allele counting (default = 5)', required=False)
    to_parser.add_argument('-mr', '--allele_min_reads', type=int,  default=10, help='Minimum number of reads required per het SNP site for allele counting (default = 10)', required=False)
    to_parser.add_argument('-acw','--ac_window', type=int, default=1200000, help='Window size for allele counting to parallelise (default = 1200000; this should be >=500000)', required=False)
    to_parser.add_argument('-w', '--cn_binsize', type=int, default=10, help='Bin window size in kbp', required=False)
    to_parser.add_argument('-b', '--blacklist', type=str, help='Path to the blacklist file', required=False)
    to_parser.add_argument('-bp', '--breakpoints', type=str, help='Path to SAVANA VCF file to incorporate savana breakpoints into copy number analysis', required=False)
    to_parser.add_argument('-c', '--chromosomes', nargs='+', default='all', help='Chromosomes to analyse. To run on all chromosomes, leave unspecified (default). To run on a subset of chromosomes only, specify the chromosome numbers separated by spaces. For x and y chromosomes, use 23 and 24, respectively.  E.g. use "-c 1 4 23 24" to run chromosomes 1, 4, X and Y', required=False)
    to_parser.add_argument('-rq', '--readcount_mapq', type=int,  default=5, help='Mapping quality threshold for reads to be included in the read counting (default = 5)', required=False)
    to_parser.add_argument('--no_blacklist', dest='blacklisting', action='store_false')
    to_parser.set_defaults(blacklisting=True)
    to_parser.add_argument('-blt', '--bl_threshold', type=int,  default='5', help='Percentage overlap between bin and blacklist threshold to tolerate for read counting (default = 5). Please specify percentage threshold as integer, e.g. "-t 5". Set "-t 0" if no overlap with blacklist is to be tolerated', required=False)
    to_parser.add_argument('--no_basesfilter', dest='bases_filter', action='store_false')
    to_parser.set_defaults(bases_filter=True)
    to_parser.add_argument('-bt', '--bases_threshold', type=int,  default='75', help='Percentage of known bases per bin required for read counting (default = 75). Please specify percentage threshold as integer, e.g. "-bt 95" ', required=False)
    to_parser.add_argument('-sl', '--smoothing_level', type=int, default='10', help='Size of neighbourhood for smoothing.', required=False)
    to_parser.add_argument('-tr', '--trim', type=float, default='0.025', help='Trimming percentage to be used.', required=False)
    to_parser.add_argument('-ms', '--min_segment_size', type=int,  default=5, help='Minimum size for a segement to be considered a segment (default = 5).', required=False)
    to_parser.add_argument('-sf', '--shuffles', type=int,  default=1000, help='Number of permutations (shuffles) to be performed during CBS (default = 1000).', required=False)
    to_parser.add_argument('-ps', '--p_seg', type=float,  default=0.05, help='p-value used to test segmentation statistic for a given interval during CBS using (shuffles) number of permutations (default = 0.05).', required=False)
    to_parser.add_argument('-pv', '--p_val', type=float,  default=0.01, help='p-value used to test validity of candidate segments from CBS using (shuffles) number of permutations (default = 0.01).', required=False)
    to_parser.add_argument('-qt', '--quantile', type=float,  default=0.2, help='Quantile of changepoint (absolute median differences across all segments) used to estimate threshold for segment merging (default = 0.2; set to 0 to avoid segment merging).', required=False)
    to_parser.add_argument('--min_ploidy', type=float, default=1.5, help='Minimum ploidy to be considered for copy number fitting.', required=False)
    to_parser.add_argument('--max_ploidy', type=float, default=5, help='Maximum ploidy to be considered for copy number fitting.', required=False)
    to_parser.add_argument('--ploidy_step', type=float, default=0.01, help='Ploidy step size for grid search space used during for copy number fitting.', required=False)
    to_parser.add_argument('--min_cellularity', type=float, default=0.2, help='Minimum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided, this is estimated during copy number fitting. Alternatively, a purity value can be provided if the purity of the sample is already known.', required=False)
    to_parser.add_argument('--max_cellularity', type=float, default=1, help='Maximum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided, this is estimated during copy number fitting. Alternatively, a purity value can be provided if the purity of the sample is already known.', required=False)
    to_parser.add_argument('--cellularity_step', type=float, default=0.01, help='Cellularity step size for grid search space used during for copy number fitting.', required=False)
    to_parser.add_argument('--cellularity_buffer', type=float, default=0.1, help='Cellularity buffer to define purity grid search space during copy number fitting (default = 0.1).', required=False)
    to_parser.add_argument('--overrule_cellularity', type=float, default=None, help='Set to sample`s purity if known. This value will overrule the cellularity estimated using hetSNP allele counts (not used by default).', required=False)   
    to_parser.add_argument('--distance_function', type=str, default='RMSD', help='Distance function to be used for copy number fitting.', choices=['RMSD', 'MAD'], required=False)
    to_parser.add_argument('--distance_filter_scale_factor', type=float, default=1.25, help='Distance filter scale factor to only include solutions with distances < scale factor * min(distance).', required=False)
    to_parser.add_argument('--distance_precision', type=int, default=3, help='Number of digits to round distance functions to', required=False)
    to_parser.add_argument('--max_proportion_zero', type=float, default=0.1, help='Maximum proportion of fitted copy numbers to be tolerated in the zero or negative copy number state.', required=False)
    to_parser.add_argument('--min_proportion_close_to_whole_number', type=float, default=0.5, help='Minimum proportion of fitted copy numbers sufficiently close to whole number to be tolerated for a given fit.', required=False)
    to_parser.add_argument('--max_distance_from_whole_number', type=float, default=0.25, help='Distance from whole number for fitted value to be considered sufficiently close to nearest copy number integer.', required=False)
    to_parser.add_argument('--main_cn_step_change', type=int, default=1, help='Max main copy number step change across genome to be considered for a given solution.', required=False)
    to_parser.add_argument('--min_block_size', type=int, default=10, help='Minimum size (number of SNPs) for a genomic block to be considered for purity estimation.', required=False)
    to_parser.add_argument('--min_block_length', type=int, default=500000, help='Minimum length (bps) for a genomic block to be considered for purity estimation.', required=False)
    to_parser.set_defaults(func=savana_tumour_only)

    try:
        global_parser.exit_on_error = False
        subparser = global_parser.parse_args().command if not args else global_parser.parse_args(args).command
    except argparse.ArgumentError as _:
        # unable to parse args, set args to None
        subparser = None

    if not subparser:
        # arguments for default, main savana process: run, classify, evaluate, cna
        global_parser.add_argument('-t', '--tumour', nargs='?', type=str, required=True, help='Tumour BAM file (must have index)')
        global_parser.add_argument('-n','--normal', nargs='?', type=str, required=False, help='Normal BAM file (must have index)')
        global_parser.add_argument('--ref', nargs='?', type=str, required=True, help='Full path to reference genome')
        global_parser.add_argument('--ref_index', nargs='?', type=str, required=False, help='Full path to reference genome fasta index (ref path + ".fai" by default)')
        global_parser.add_argument('--contigs', nargs='?', type=str, help="Contigs/chromosomes to consider (optional, default=All)")
        global_parser.add_argument('--length', nargs='?', type=int, default=30, help='Minimum length SV to consider (default=30)')
        global_parser.add_argument('--mapq', nargs='?', type=int, default=5, help='Minimum MAPQ to consider a read mapped (default=5)')
        global_parser.add_argument('--buffer', nargs='?', type=int, default=10, help='Buffer when clustering adjacent potential breakpoints, excepting insertions (default=10)')
        global_parser.add_argument('--insertion_buffer', nargs='?', type=int, default=250, help='Buffer when clustering adjacent potential insertion breakpoints (default=250)')
        global_parser.add_argument('--end_buffer', nargs='?', type=int, default=50, help='Buffer when clustering alternate edge of potential breakpoints, excepting insertions (default=50)')
        global_parser.add_argument('--threads', nargs='?', type=int, const=0, help='Number of threads to use (default=max)')
        global_parser.add_argument('--outdir', nargs='?', required=True, help='Output directory (can exist but must be empty)')
        global_parser.add_argument('--overwrite', action='store_true', help='Use this flag to write to output directory even if files are present')
        global_parser.add_argument('--sample', nargs='?', type=str, help='Name to prepend to output files (default=tumour BAM filename without extension)')
        global_parser.add_argument('--keep_inv_artefact', action='store_true', help="Do not remove breakpoints with foldback-inversion artefact pattern")
        global_parser.add_argument('--single_bnd', action='store_true', help='Report single breakend variants in addition to standard types (default=False)')
        global_parser.add_argument('--single_bnd_min_length', nargs='?', type=int, default=1000, help='Minimum length of single breakend to consider (default=100)')
        global_parser.add_argument('--single_bnd_max_mapq', nargs='?', type=int, default=5, help='Convert supplementary alignments below this threshold to single breakend (default=5, must not exceed --mapq argument)')
        global_parser.add_argument('--debug', action='store_true', help='Output extra debugging info and files')
        global_parser.add_argument('--chunksize', nargs='?', type=int, default=100000000, help='Chunksize to use when splitting genome for parallel analysis - used to optimise memory (default=1000000)')
        global_parser.add_argument('--coverage_binsize', nargs='?', type=int, default=5, help='Length used for coverage bins (default=5)')
        # classify args
        global_parser.add_argument('--min_support', nargs='?', type=int, default=3, required=False, help='Minimum supporting reads for a PASS variant (default=3)')
        global_parser.add_argument('--min_af', nargs='?', type=helper.float_range(0.0, 1.0), default=0.01, required=False, help='Minimum allele-fraction for a PASS variant (default=0.01)')
        global_parser.add_argument('--cna_rescue', nargs='?', type=str, required=False, help='Copy number abberation output file for this sample (used to rescue variants)')
        global_parser.add_argument('--cna_rescue_distance', nargs='?', type=int, default=50, required=False, help='Maximum distance from a copy number abberation for a variant to be rescued')
        classify_group = global_parser.add_mutually_exclusive_group()
        classify_group.add_argument('--ont', action='store_true', help='Use the Oxford Nanopore (ONT) trained model to classify variants (default)')
        classify_group.add_argument('--pb', action='store_true', help='Use PacBio thresholds to classify variants')
        # whether to use a germline-trained model
        global_parser.add_argument('--predict_germline', action='store_true', help='Also predict germline events')
        classify_group.add_argument('--custom_model', nargs='?', type=str, required=False, help='Pickle file of custom machine-learning model')
        classify_group.add_argument('--custom_params', nargs='?', type=str, required=False, help='JSON file of custom filtering parameters')
        classify_group.add_argument('--legacy', action='store_true', help='Use legacy lenient/strict filtering')
        global_parser.add_argument('--somatic_output', nargs='?', type=str, required=False, help='Output VCF with only PASS somatic variants')
        global_parser.add_argument('--germline_output', nargs='?', type=str, required=False, help='Output VCF with only PASS germline variants')
        global_parser.add_argument('--confidence', nargs='?', type=helper.float_range(0.0, 1.0), default=None, help='Confidence level for mondrian conformal prediction - suggested range (0.70-0.99) (not used by default)')
        # evaluate args
        global_parser.add_argument('--somatic', nargs='?', type=str, required=False, help='Somatic VCF file to evaluate against')
        global_parser.add_argument('--germline', nargs='?', type=str, required=False, help='Germline VCF file to evaluate against (optional)')
        global_parser.add_argument('--overlap_buffer', nargs='?', type=int, default=100, required=False, help='Buffer for considering an overlap (default=100)')
        global_parser.add_argument('--curate', action='store_true', default=False, help='Attempt to reduce false labels for training (allow label to be used twice)')
        global_parser.add_argument('--qual_filter', nargs='?', type=int, default=None, help='Impose a quality filter on comparator variants (default=None)')
        evaluate_group = global_parser.add_mutually_exclusive_group()
        evaluate_group.add_argument('--by_support', action='store_true', help='Comparison method: tie-break by read support')
        evaluate_group.add_argument('--by_distance', action='store_true', default=True, help='Comparison method: tie-break by min. distance (default)')
        # cna args
        global_parser.add_argument('--cna_threads', nargs='?', type=int, const=0, help='Number of threads to use for CNA (default=max)')
        global_parser.add_argument('--tmpdir', nargs='?', required=False, default='tmp', help='Temp directory for allele counting temp files (defaults to outdir)')
        allele_group = global_parser.add_mutually_exclusive_group()
        allele_group.add_argument('-v', '--snp_vcf', type=str, help='Path to SNP vcf file to extract heterozygous SNPs for allele counting.', required=False)
        allele_group.add_argument('-vg','--g1000_vcf', type=str, choices={"1000g_hg38", "1000g_hg19", "1000g_t2t"}, help='Use 1000g biallelic vcf file for allele counting instead of SNP vcf from matched normal. Specify which genome version to use.', required=False)
        allele_group.add_argument('-ac', '--allele_counts_het_snps', type=str, help='Path to allele counts of heterozygous SNPs', required=False)
        global_parser.add_argument('-q', '--allele_mapq', type=int,  default=5, help='Mapping quality threshold for reads to be included in the allele counting (default = 5)', required=False)
        global_parser.add_argument('-mr', '--allele_min_reads', type=int,  default=10, help='Minimum number of reads required per het SNP site for allele counting (default = 10)', required=False)
        global_parser.add_argument('-acw','--ac_window', type=int, default=1200000, help='Window size for allele counting to parallelise (default = 1200000; this should be >=500000)', required=False)
        global_parser.add_argument('-w', '--cn_binsize', type=int, default=10, help='Bin window size in kbp', required=False)
        global_parser.add_argument('-b', '--blacklist', type=str, help='Path to the blacklist file', required=False)
        global_parser.add_argument('-bp', '--breakpoints', type=str, help='Path to SAVANA VCF file to incorporate savana breakpoints into copy number analysis', required=False)
        global_parser.add_argument('-c', '--chromosomes', nargs='+', default='all', help='Chromosomes to analyse. To run on all chromosomes, leave unspecified (default). To run on a subset of chromosomes only, specify the chromosome numbers separated by spaces. For x and y chromosomes, use 23 and 24, respectively.  E.g. use "-c 1 4 23 24" to run chromosomes 1, 4, X and Y', required=False)
        global_parser.add_argument('-rq', '--readcount_mapq', type=int,  default=5, help='Mapping quality threshold for reads to be included in the read counting (default = 5)', required=False)
        global_parser.add_argument('--no_blacklist', dest='blacklisting', action='store_false')
        global_parser.set_defaults(blacklisting=True)
        global_parser.add_argument('-blt', '--bl_threshold', type=int,  default='5', help='Percentage overlap between bin and blacklist threshold to tolerate for read counting (default = 0, i.e. no overlap tolerated). Please specify percentage threshold as integer, e.g. "-t 5" ', required=False)
        global_parser.add_argument('--no_basesfilter', dest='bases_filter', action='store_false')
        global_parser.set_defaults(bases_filter=True)
        global_parser.add_argument('-bt', '--bases_threshold', type=int,  default='75', help='Percentage of known bases per bin required for read counting (default = 0, i.e. no filtering). Please specify percentage threshold as integer, e.g. "-bt 95" ', required=False)
        global_parser.add_argument('-sl', '--smoothing_level', type=int, default='10', help='Size of neighbourhood for smoothing.', required=False)
        global_parser.add_argument('-tr', '--trim', type=float, default='0.025', help='Trimming percentage to be used.', required=False)
        global_parser.add_argument('-ms', '--min_segment_size', type=int,  default=5, help='Minimum size for a segement to be considered a segment (default = 5).', required=False)
        global_parser.add_argument('-sf', '--shuffles', type=int,  default=1000, help='Number of permutations (shuffles) to be performed during CBS (default = 1000).', required=False)
        global_parser.add_argument('-ps', '--p_seg', type=float,  default=0.05, help='p-value used to test segmentation statistic for a given interval during CBS using (shuffles) number of permutations (default = 0.05).', required=False)
        global_parser.add_argument('-pv', '--p_val', type=float,  default=0.01, help='p-value used to test validity of candidate segments from CBS using (shuffles) number of permutations (default = 0.01).', required=False)
        global_parser.add_argument('-qt', '--quantile', type=float,  default=0.2, help='Quantile of changepoint (absolute median differences across all segments) used to estimate threshold for segment merging (default = 0.2; set to 0 to avoid segment merging).', required=False)
        global_parser.add_argument('--min_ploidy', type=float, default=1.5, help='Minimum ploidy to be considered for copy number fitting.', required=False)
        global_parser.add_argument('--max_ploidy', type=float, default=5, help='Maximum ploidy to be considered for copy number fitting.', required=False)
        global_parser.add_argument('--ploidy_step', type=float, default=0.01, help='Ploidy step size for grid search space used during for copy number fitting.', required=False)
        global_parser.add_argument('--min_cellularity', type=float, default=0.2, help='Minimum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided, this is estimated during copy number fitting. Alternatively, a purity value can be provided if the purity of the sample is already known.', required=False)
        global_parser.add_argument('--max_cellularity', type=float, default=1, help='Maximum cellularity to be considered for copy number fitting. If hetSNPs allele counts are provided, this is estimated during copy number fitting. Alternatively, a purity value can be provided if the purity of the sample is already known.', required=False)
        global_parser.add_argument('--cellularity_step', type=float, default=0.01, help='Cellularity step size for grid search space used during for copy number fitting.', required=False)
        global_parser.add_argument('--cellularity_buffer', type=float, default=0.1, help='Cellularity buffer to define purity grid search space during copy number fitting (default = 0.1).', required=False)
        global_parser.add_argument('--overrule_cellularity', type=float, default=None, help='Set to sample`s purity if known. This value will overrule the cellularity estimated using hetSNP allele counts (not used by default).', required=False)   
        global_parser.add_argument('--distance_function', type=str, default='RMSD', help='Distance function to be used for copy number fitting.', choices=['RMSD', 'MAD'], required=False)
        global_parser.add_argument('--distance_filter_scale_factor', type=float, default=1.25, help='Distance filter scale factor to only include solutions with distances < scale factor * min(distance).', required=False)
        global_parser.add_argument('--distance_precision', type=int, default=3, help='Number of digits to round distance functions to', required=False)
        global_parser.add_argument('--max_proportion_zero', type=float, default=0.1, help='Maximum proportion of fitted copy numbers to be tolerated in the zero or negative copy number state.', required=False)
        global_parser.add_argument('--min_proportion_close_to_whole_number', type=float, default=0.5, help='Minimum proportion of fitted copy numbers sufficiently close to whole number to be tolerated for a given fit.', required=False)
        global_parser.add_argument('--max_distance_from_whole_number', type=float, default=0.25, help='Distance from whole number for fitted value to be considered sufficiently close to nearest copy number integer.', required=False)
        global_parser.add_argument('--main_cn_step_change', type=int, default=1, help='Max main copy number step change across genome to be considered for a given solution.', required=False)
        global_parser.add_argument('--min_block_size', type=int, default=10, help='Minimum size (number of SNPs) for a genomic block to be considered for purity estimation.', required=False)
        global_parser.add_argument('--min_block_length', type=int, default=500000, help='Minimum length (bps) for a genomic block to be considered for purity estimation.', required=False)
        global_parser.set_defaults(func=savana_main)
        parsed_args = global_parser.parse_args() if not args else global_parser.parse_args(args)
    else:
        global_parser.exit_on_error = True
        parsed_args = global_parser.parse_args() if not args else global_parser.parse_args(args)

    return parsed_args

def main(args=None):
    """ main function for SAVANA - collects command line arguments and executes algorithm """
    print(logo)
    print(f'Version {helper.__version__}')
    src_location = __file__
    print(f'Source: {src_location}\n')
    args = parse_args(args)
    args.func(args)

if __name__ == "__main__":
    main()
