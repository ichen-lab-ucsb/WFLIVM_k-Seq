#!/usr/bin/env python3
"""Script for count data preprocessing. Use 'count-data-preprocessing.py -h' for usages and examples"""

from argparse import ArgumentParser, RawTextHelpFormatter
import pandas as pd
from k_seq.data import seq_data, filters, landscape
from k_seq.data.transform import ReactedFractionNormalizer
from pathlib import Path
import sys


def get_seq_in_double_mutants(dataset):
    """Get the sequences only within double mutants to family centers"""
    # filter up to double mutants
    pool_peaks = {
        'pk2.1': 'ATTACCCTGGTCATCGAGTGA',
        'pk2.2': 'ATTCACCTAGGTCATCGGGTG',
        'pk1A.1': 'CTACTTCAAACAATCGGTCTG',
        'pk1B.1': 'CCACACTTCAAGCAATCGGTC',
        'pk3.1': 'AAGTTTGCTAATAGTCGCAAG'
    }

    if not hasattr(dataset.table, 'filtered'):
        dataset.table.filtered = filters.NoAmbiguityFilter.filter(
            target=filters.SeqLengthFilter.filter(min_len=21, max_len=21, target=dataset.table.original, axis=0),
            axis=0
        )

    dataset.pool_peaks = [
        landscape.Peak(seqs=dataset.table.filtered, center_seq=seq,
                       name=name, dist_type='hamming') for name, seq in pool_peaks.items()
    ]
    dataset.pool_peaks_merged = landscape.PeakCollection(peaks=dataset.pool_peaks)

    peak_filter = filters.PeakFilter(max_dist=2,
                                     dist_to_center=dataset.pool_peaks_merged.dist_to_center)
    return peak_filter(dataset.table.filtered).index


def prepare_median_input(input_count, rna_amount, save_to):
    """Preprocess input samples from count data to median input amount (ng)"""
    print('Preparing median input amount from count files...')
    print(f'Loading count files form {input_count}...')
    dataset_input = seq_data.SeqData.from_count_files(
        count_files=input_count,
        pattern_filter='counts',
        name_template='[input_{input_id}]_counts.txt',
        input_sample_name=None,
        sort_by='name',
        note=f"",
        x_values='input_id'
    )

    print('Quantifying amount...')
    input_rna = pd.read_csv(rna_amount, index_col=0).iloc[:, 0][dataset_input.samples]
    input_median = (dataset_input.table.original / dataset_input.table.original.sum(axis=0) * input_rna).median(axis=1).rename('ng')
    input_median.index.name = 'seq'

    print('Filtering sequences up to double mutants...')
    # Only include sequences within double mutants
    seq_in_double_mutants = get_seq_in_double_mutants(dataset_input)
    save_to = Path(save_to)/'input-rna-median-ng.csv'
    input_median[seq_in_double_mutants.values].to_csv(save_to)
    print(f'median input amount saved to {str(save_to)}')


def prepare_seqdata(substrate, reacted_count, rna_amount, input_median, save_to):
    """Preprocess reacted samples from count data to reacted fraction"""

    print(f'Processing dataset for {substrate.upper()}...')

    print(f'Loading count files from {Path(reacted_count).resolve()}...')
    dataset = seq_data.SeqData.from_count_files(
        count_files=reacted_count,
        pattern_filter='counts',
        name_template=substrate[1].upper() + '_[{sub, float}_{input_id}]_counts.txt',
        input_sample_name=None,
        sort_by='name',
        note=f"{substrate.upper()} doped pool k-seq results (up to double mutants)",
        x_values='sub',
        x_unit='M'
    )
    dataset.x_values = dataset.x_values * 1e-6   # uM -> M

    # quantify by total RNA from qPCR
    print('Quantifying amount...')
    total_rna = pd.read_csv(rna_amount, index_col=0).loc[dataset.samples.values, substrate.lower()]
    dataset.add_sample_total(full_table=dataset.table.original,
                             total_amounts=total_rna,
                             unit='ng')

    print('Filtering tables')
    # remove sequences are not 21 length and sequences contains ambiguous nt
    dataset.table.filtered = filters.NoAmbiguityFilter.filter(
        target=filters.SeqLengthFilter.filter(min_len=21, max_len=21, target=dataset.table.original, axis=0),
        axis=0
    )

    seq_in_double = get_seq_in_double_mutants(dataset)

    # import input median over replicates
    all_inputs_ng = pd.read_csv(input_median, index_col=0)['ng']
    amount = dataset.sample_total(dataset.table.filtered).loc[seq_in_double]
    amount['input'] = all_inputs_ng.reindex(amount.index)
    dataset.table.reacted_fraction = ReactedFractionNormalizer(input_samples=['input'])(amount)
    save_to = Path(save_to)/f'{substrate.lower()}-preped.pkl'
    dataset.to_pickle(save_to)
    print(f'Preprocessed {substrate.upper()} data saved to {str(save_to)}')


if __name__ == '__main__':

    parser = ArgumentParser(
        prog="BXO count data preprocessing",
        description="Script to preprocess input (unreacted) or reacted count data for k-Seq\n"
                    "Examples:\n\n"
                    "Preprocess reacted samples from count file to reacted fraction (pickled dataset for k-seq package):\n"
                    "count-data-preprocessing.py \ \n"
                    "    --substrate BFO \ \n"
                    "    --reacted_count /path/to/count/folder \ \n"
                    "    --rna_amount /path/to/rna-ng-merged.csv \ \n"
                    "    --input_median /path/to/input-rna-median-ng.csv \ \n"
                    "    --output /path/to/output/folder\n\n"
                    "Preprocess input samples from count file to input-rna-median-ng.csv (median input amount):\n"
                    "count-data-preprocessing.py \ \n"
                    "    --input_count /path/to/input/count/folder \ \n"
                    "    --rna_amount /path/to/rna-ng-merged.csv \ \n"
                    "    --output /path/to/output/folder",
        formatter_class=RawTextHelpFormatter
    )
    parser.add_argument('--substrate', type=str, default=None,
                        help='BXO substrate. E.g. BFO, BLO, ...')
    parser.add_argument('--reacted_count', type=str, default=None,
                        help='Folder contains count files for reacted samples')
    parser.add_argument('--input_median', type=str, default=None,
                        help='CSV file contains the median input amount of sequences')
    parser.add_argument('--rna_amount', type=str, default=None,
                        help="the csv file contain sample RNA amount measured by qPCR")
    parser.add_argument('--input_count', type=str, default=None,
                        help='Folder contains count files for input samples')
    parser.add_argument('--output', type=str, default='.',
                        help="Directory to save preprocessed results")

    args = parser.parse_args()

    if args.reacted_count is None and args.input_count is None:
        print('Please provide --reacted_count (for preprocessing reacted samples) or '
              '--input_count (for processing input samples)')
        print('-' * 20)
        parser.print_help()
        sys.exit()
    elif args.reacted_count is not None and args.input_count is not None:
        print('Count not specify both --reacted_count and --input_count. Choose one of the tasks')

    if args.input_count is not None:
        prepare_median_input(args.input_count, rna_amount=args.rna_amount, save_to=args.output)
    else:
        prepare_seqdata(args.substrate, args.reacted_count, args.rna_amount, args.input_median, args.output)

