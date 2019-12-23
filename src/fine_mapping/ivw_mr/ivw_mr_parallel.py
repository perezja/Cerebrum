#!/usr/bin/env python3

import argparse
import subprocess
import pandas as pd
import numpy as np
import tempfile
import shutil
import queue
import json
import time
import logging
from common import *
from collections import defaultdict
import os
import re
import sys

def parse_args():

    parser = argparse.ArgumentParser(description='Get pruned list of variants to keep.')
    parser.add_argument('bfile', help='Plink bfile containing variant data for all participants.')
    parser.add_argument('eqtl_sumstats', help='Summary statistics for exposure.') 
    parser.add_argument('outcome_file', help='Outcome file mapping outcome ID to file.')
    parser.add_argument('--prefix', default='out', help='Outfile prefix.')
    parser.add_argument('--chunksize', type=int, default=50, help='Subset of exposures.')
    parser.add_argument('--log-level', default='INFO', choices=['NOTSET','DEBUG','INFO', 'WARNING','CRITICAL',
                                                        'ERROR','CRITICAL'], help='Log level')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory.')

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    args.output_dir = os.path.abspath(args.output_dir)

    log.setLevel(args.log_level)
    log.info(sys.argv)    

    return(args)

def run_ivw_chunk(arg_list, log_fp):

    exc = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ivw_mr_run_chunk.py')
    cmd = ' '.join([exc] + arg_list)

    return(subprocess.Popen(bsub(cmd, 'ivw_mr', log_fp), shell=True))

def main(args):
         
    chunk_file_dir = tempfile.mkdtemp(dir=args.output_dir)
    log_file_dir = os.path.join(args.output_dir, 'logs') 
    mr_outdir = os.path.join(args.output_dir, 'mr')

    if not os.path.exists(log_file_dir):
        os.makedirs(log_file_dir)
    if not os.path.exists(mr_outdir):
        os.makedirs(mr_outdir)

    # get outcomes
    with open(args.outcome_file) as fp:
        outcome_dat = {row.split('\t')[0]: row.split('\t')[1] for row in fp.read().strip().split('\n')}

    print('\rReading eQTL data...', end='', flush=True)
    df = pd.read_csv(args.eqtl_sumstats, sep='\t')

    total_submissions = 0
    chunk_dict = defaultdict(lambda: dict())
    chunk_processes = queue.Queue()

    print('\rPreparing inputs for submission...', end='', flush=True)
    for chunk_name, chunk_df in loop_chunks(df, args.chunksize, chunk_dict) :

        # Write out exposure files constituting a chunk 

        exposure_outfile = os.path.join(args.output_dir, chunk_name +'.tmp')

        chunk_df.to_csv(exposure_outfile, sep='\t', index=False) 

        for outcome_id, outcome_fp in outcome_dat.items():

            chk_id = chunk_name +'_'+ outcome_id
    
            log_outfile = os.path.join(log_file_dir, chk_id + '.log')
    
            arg_list = [args.bfile, exposure_outfile, outcome_fp, outcome_id, '-o', mr_outdir, '--prefix', chk_id] 
            po = run_ivw_chunk(arg_list, log_outfile)

            chunk_processes.put((chk_id, po))
            chunk_dict[chk_id]['log'] = log_outfile

            total_submissions += 1

    print('submitted {} jobs.'.format(total_submissions))
   
if __name__ == '__main__':
    args = parse_args()
    main(args)
