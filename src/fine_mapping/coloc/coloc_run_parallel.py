#!/usr/bin/env python3

import argparse
import subprocess
import pandas as pd
import numpy as np
from collections import defaultdict
import json
import queue
import tempfile
import time
import shutil
import logging
import glob
import re
import os
import sys

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)

log = logging.getLogger(__name__)

def parse_args():

    parser = argparse.ArgumentParser(description='Run TwoSampleMR on a set of exposure files.')
    parser.add_argument('locus_eqtl_dir', help='Directory with locus eQTL data.')
    parser.add_argument('eqtl_freq_file', help='')
    parser.add_argument('outcome_file', help='File mapping outcome ID to outcome file path')
    parser.add_argument('-o', '--output_dir', help='Output directory.')

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    args.output_dir = os.path.abspath(args.output_dir)

    log.setLevel('INFO')
    log.info(sys.argv)    

    return(args)
    
def bsub(cmd, job_name, log_fp, mem=8, gtmp=2, docker_image="apollodorus/eqtl-coloc:v1", queue="research-hpc"):
    """ Creates a bsub command for LSF job submission.
    Args:
        cmd: command to be run.
        queue: queue to submit job
        job_name: name of job
        mem: RAM space to request from scheduler (MB)
        gtmp: temporary space to request from scheduler (GB)
        docker_image: [username]/[repos]:[tag]
    Returns:
        string: bsub command string.
    """

    mem = str(mem)
    gtmp = str(gtmp)
    bsub = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q "+queue \
        + " -J "+job_name \
        + " -M "+mem+"000000" \
        + " -o "+log_fp \
        + " -R 'select[mem>"+mem+"000 && gtmp > "+gtmp+"]" \
        + " rusage[mem="+mem+"000, gtmp="+gtmp+"]'" \
        + " -a 'docker("+docker_image+")'" \
        + " /bin/bash -c '"+cmd+"'" 

    
    return(bsub)

def run_coloc(arg_list, log_fp):

    exc = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_coloc.R')
    cmd = ' '.join([exc] + arg_list)

    with open(log_fp, 'w+') as fp:
        return(subprocess.Popen(bsub(cmd, 'coloc', log_fp), stdout=fp, stderr=subprocess.STDOUT, shell=True))

def merge_chunks(input_dir, outfile):
    """ process to concatenate all finished chunk output files."""

    header = ['locus_name', 'sentinel_snp', 'sentinel_snp_genes', 'outcome', 'nsnps', 'PP.H0.abf', 'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf']
    with open(outfile, mode='w+') as fp:
        fp.write('\t'.join(header)+'\n')

    # find: full path to all files in specified directory
    # tail: cat all but header line
    cmd = "find {} -maxdepth 1 -type f | sort -V | xargs -I{{}} tail -n +2 {{}} >> {}".format(input_dir, outfile) 

    subprocess.check_call(bsub(cmd, 'merge_res'), shell=True)


def main(args):

    log_file_dir = os.path.join(args.output_dir, 'logs') 

    if not os.path.exists(log_file_dir):
        os.makedirs(log_file_dir)

    loci_sumstats = glob.glob(os.path.join(args.locus_eqtl_dir, '*'))

    outcome_data = dict()
    with open(args.outcome_file) as fp:
        outcome_data = {row.split('\t')[0]: row.split('\t')[1] for row in fp.read().strip().split('\n')}
        
    total_jobs = 0
    for k,locus_eqtl_fp in enumerate(loci_sumstats):
        for outcome_id,outcome_fp in outcome_data.items():

            of_prefix = os.path.split(locus_eqtl_fp)[1] + '_' + outcome_id
            outfile = os.path.join(args.output_dir, of_prefix+'.coloc')

            log_outfile = os.path.join(log_file_dir, of_prefix +'.log')

            arg_list = [locus_eqtl_fp, args.eqtl_freq_file, outcome_fp, outcome_id, '-o', outfile] 
            run_coloc(arg_list, log_outfile)
            total_jobs+=1

    print('Submitted {} jobs for bayesian COLOC analysis.'.format(total_jobs))

if __name__ == '__main__':
    args = parse_args()
    main(args)
