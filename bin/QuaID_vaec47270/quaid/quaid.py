# Standard library imports
import logging
import os
# import glob
import argparse as ap
# 3rd party imports
from tqdm import tqdm
import pandas as pd
# QuaID imports
from utils.file_loaders import *
from utils.quaid_func import *

def main():
    parser = ap.ArgumentParser()
    parser.add_argument('vdb_path', help='path to the aggregated VDB output data')
    parser.add_argument('metadata_path', help='path to the metadata file')
    parser.add_argument('vcf_folder_path', help='path to the folder containing variant calling information')
    parser.add_argument('coverage_folder', help='path to the folder containing coverage information')
    parser.add_argument('--output-folder', help='path to the output folder', 
                        default='./output')
    parser.add_argument('--concat-output', help='create a single summary table instead of separate files per date',
                        action='store_true')
    parser.add_argument('--occurence-cutoff', help='minimum number of genomes required per lineage (lineages below this count are excluded)',
                        type=int, default=0)
    parser.add_argument('--inclusion-cutoff', help='minimum mutation prevalence fraction in a lineage for the inclusion into the quasi-unique set',
                        type=float, default=0.5)
    parser.add_argument('--exclusion-cutoff', help='minimum mutation prevalence fraction in a non-target lineage used for exclusion from the quasi-unique set',
                        type=float, default=0.5)
    parser.add_argument('--time-window', help='number of prior weeks of sequence data to consider (or -1 to use YTD)',
                        type=int, default=4)
    parser.add_argument('--level', help="PANGO lineage level for quasi-unique mutations (e.g. B.1.1.7 is level 4)",
                        type=int, default=4)
    parser.add_argument('--save_qu_to_disk', help='path to a file where to save quasi-unique mutations DataFrame, if None the file is not saved to disk',
                        default=None)
    parser.add_argument('--save-vcf', help='save individual mutation calls with lineage annotation',
                        action='store_true')
    parser.add_argument('--verbose', help='print progress information to stdout (useful for progress monitoring)',
                        action='store_true')

    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(level=(logging.DEBUG if args.verbose else logging.ERROR), 
                        format='%(asctime)s [%(name)s - %(levelname)s]: %(message)s', 
                        datefmt='%m/%d/%Y %I:%M:%S %p')
    logger = logging.getLogger('QuaID')

    # Load mutation information, metadata, and variant calling results
    vdb_df = load_vdb_df(args.vdb_path)
    metadata = load_metadata(args.metadata_path)
    vcf_df = load_vcf(args.vcf_folder_path)
    coverage_data = load_coverage_data(args.coverage_folder)
    alias_data = load_alias_data()

    # Preprocess data
    logger.info('Pre-processing count data')
    genome_counts = get_counts(metadata)
    WHO_variants = get_all_voc(args.level)
    logger.info('Merging count and variant calling data')
    count_df = build_count_df(genome_counts, vdb_df)

    dates = sorted(vcf_df.Date.unique())
    results = []
    logger.info('Analyzing samples')
    for t, _date in enumerate(tqdm(dates, colour='#CCCCCC', disable=(not args.verbose))):
        recent_nt_df = get_recent_nt_df(count_df, _date, time_window=args.time_window)
        quasiunique_df = get_quasiunique_df(recent_nt_df,
                                            WHO_variants,
                                            alias_data,
                                            occurence_cutoff=args.occurence_cutoff,
                                            incl_cutoff=args.inclusion_cutoff,
                                            excl_cutoff=args.exclusion_cutoff,
                                            level=args.level,
                                            save_to_disk=args.save_qu_to_disk,
                                            verbose=args.verbose,
                                            date=_date)
        result_table = get_result_table(vcf_df, quasiunique_df, coverage_data, WHO_variants, _date, 
                                            save=args.save_vcf, 
                                            verbose=args.verbose)
        
        if args.concat_output:
            results.append(result_table)
        else:
            date_str = pd.to_datetime(_date).strftime('%Y-%m-%d')
            try:                
                os.mkdir(f"{args.output_folder}/{date_str}")
            except FileExistsError:
                pass
            name_str = f"Summary-{pd.to_datetime(_date).strftime('%Y-%m-%d')}-o{args.occurence_cutoff}-i{args.inclusion_cutoff:.2}-e{args.exclusion_cutoff:.2}-t{args.time_window}-l{args.level}-QU-noX-mask.csv"
            result_table.to_csv(f"{args.output_folder}/{date_str}/{name_str}", index=False)

    if args.concat_output:
        logger.info(f'Writing summary output to {args.output_folder}')
        result_table = pd.concat(results, ignore_index=True)
        date_str = f"{pd.to_datetime(dates[0]).strftime('%Y-%m-%d')}-{pd.to_datetime(dates[-1]).strftime('%Y-%m-%d')}"
        name_str = f"Summary-{date_str}-o{args.occurence_cutoff}-i{args.inclusion_cutoff:.2}-e{args.exclusion_cutoff:.2}-t{args.time_window}-l{args.level}-QU-noX-mask.csv"
        try: 
            os.mkdir(f"{args.output_folder}")
        except FileExistsError:
            pass
        result_table.to_csv(f"{args.output_folder}/{name_str}", index=False)

    exit(0)


if __name__ == '__main__':
    main()
