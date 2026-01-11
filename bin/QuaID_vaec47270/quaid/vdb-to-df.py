"""Convert a trimmed nucleotide vdb output to a pandas dataframe and export to CSV-file."""
# Standard library imports
import logging
import os
import argparse as ap
from collections import Counter
# 3rd party imports
import pandas as pd
import dask.dataframe as dd
# QuaID imports
from utils.file_loaders import load_metadata, load_vdb_mutation_data, merge_data

def convert_to_df(_vdb_df, use_dask):
    if use_dask:
        # Beyond MSA sizes of roughly 4-5 million SARS-CoV-2 sequences using Dask is necessary
        vdb_dd = dd.from_pandas(_vdb_df, npartitions=144, sort=True)  # 144 works well on our servers
        vdb_dd = vdb_dd.groupby(['Collection date', 'Pango lineage']).sum()
        vdb_lineage_df = vdb_dd.compute(scheduler='threads')

        vdb_lineage_df = vdb_lineage_df.drop(['Accession ID'], axis=1)
        vdb_lineage_dd = dd.from_pandas(vdb_lineage_df.reset_index(), npartitions=80, sort=True)  # 80 works well on our servers 
        vdb_lineage_dd['SNPs'] = vdb_lineage_dd['SNPs'].map(lambda x: Counter(x.strip().split(' ')))
        vdb_lineage_df = vdb_lineage_dd.compute(scheduler='threads')        
    else:
        vdb_lineage_df = _vdb_df.groupby(['Collection date', 'Pango lineage']).sum()
        
        vdb_lineage_df = vdb_lineage_df.drop(['Accession ID'], axis=1)
        vdb_lineage_df['SNPs'] = vdb_lineage_df['SNPs'].map(lambda x: Counter(x.strip().split(' ')))

    return vdb_lineage_df


def aggregate_counts(_vdb_lineage_df, resample_interval='W-MON'):
    _vdb_lineage_df = _vdb_lineage_df.reset_index()
    _vdb_lineage_df = _vdb_lineage_df.set_index('Collection date')
    _vdb_group_df = _vdb_lineage_df.groupby('Pango lineage')['SNPs'].resample(resample_interval).sum()

    return _vdb_group_df


def main():
    parser = ap.ArgumentParser()
    parser.add_argument('vdb_path', help='path to the aggregated VDB output data')
    parser.add_argument('metadata_path', help='path to the metadata file')
    parser.add_argument('--output', help='output files prefix', 
                        default='vdb_lineage_df')
    parser.add_argument('--output-folder', help='path to the output folder', 
                        default='./vdb_output')
    parser.add_argument('--no-dask', help='do not use Dask for any of the operations [WARNING: only use for small files]',
                        action='store_true')
    parser.add_argument('--occurence-cutoff', help='prevalence cutoff value',
                        type=int, default=0)
    parser.add_argument('--verbose', help='print progress information to stdout',
                        action='store_true')

    args = parser.parse_args()

    metadata_f = args.metadata_path
    vdb_f = args.vdb_path
    use_dask = (not args.no_dask)

    # Configure logging
    logging.basicConfig(level=(logging.DEBUG if args.verbose else logging.ERROR), 
                        format='%(asctime)s [%(name)s %(levelname)s]: %(message)s', 
                        datefmt='%m/%d/%Y %I:%M:%S %p')
    logger = logging.getLogger('VDB2DF')

    # Load mutation information and metadata
    metadata = load_metadata(metadata_f)
    mutations_data = load_vdb_mutation_data(vdb_f)

    # Preprocess data
    logger.info('Merging mutation data with metadata.')
    vdb_df = merge_data(mutations_data, metadata)

    # Convert to DataFrame for downstream analyses
    logger.info('Converting VDB output to DataFrame.')
    vdb_lineage_df = convert_to_df(vdb_df, use_dask)

    vdb_lineage_df.to_csv(f'{args.output_folder}/{args.output}_counts.csv')
    if not os.path.exists(f'{args.output_folder}/{args.output}_counts.csv'):
        logger.error(f"The file {args.output_folder}/{args.output}_counts.csv was not saved.")
        exit(1)

    # Aggregate data to the weekly level
    logger.info('Aggregating mutation data by week.')
    vdb_lineage_df = vdb_lineage_df.set_index('Collection date')
    vdb_week_group_df = aggregate_counts(vdb_lineage_df)

    vdb_week_group_df.to_csv(f'{args.output_folder}/{args.output}_week.csv')
    if not os.path.exists(f'{args.output_folder}/{args.output}_week.csv'):
        logger.error(f"The file {args.output_folder}/{args.output}_week.csv was not saved.")
        exit(1)
    
    logger.info('Output saved succesfully. Exiting...')
    exit(0)


if __name__ == '__main__':
    main()
