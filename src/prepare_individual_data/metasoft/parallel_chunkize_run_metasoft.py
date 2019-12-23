#!/usr/bin/env python3

import argparse
import subprocess
import pandas as pd
import numpy as np
from collections import defaultdict
import logging
import tempfile
import shutil
import queue
import json
import time
import os
import sys

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)

log = logging.getLogger(__name__)

def parse_args():

    parser = argparse.ArgumentParser(description='Get pruned list of variants to keep.')
    parser.add_argument('chunk_fp', help='File listing METASOFT input chunks.')
    parser.add_argument('-f', '--chunk_list', help='File listing input chunk IDs to include.')
    parser.add_argument('-c', '--chunksize', type=int, default=5*(10**6), help='Size of mini-chunks to make out of component chunk files.')
    parser.add_argument('--log-level', default='INFO', choices=['NOTSET','DEBUG','INFO', 'WARNING','CRITICAL',
                                                        'ERROR','CRITICAL'], help='Log level')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory.')

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    args.output_dir = os.path.abspath(args.output_dir)

    log.setLevel('INFO')
    log.info(sys.argv)    

    return(args)

def bsub(cmd, log, job_name, mem=1, gtmp=1, docker_image="apollodorus/brain-eqtl:eqtl-analysis", queue="research-hpc"):
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
    bsub = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false LSB_JOB_REPORT_MAIL=N bsub -q "+queue \
        + " -J "+job_name \
        + " -M "+mem+"000000" \
        + " -o "+log \
        + " -R 'select[mem>"+mem+"000 && gtmp > "+gtmp+"]" \
        + " rusage[mem="+mem+"000, gtmp="+gtmp+"]'" \
        + " -a 'docker("+docker_image+")'" \
        + " /bin/bash -c '"+cmd+"'" 

    
    return(bsub)

def run_metasoft(chunk_fp, chunk_name, outdir, log_fp):

    exc = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_metasoft_chunk.py')
    cmd = " ".join([exc, chunk_fp, '-o', outdir])

    with open(log_fp, 'w+') as fp:
        fp.write('* cmd: '+cmd+'\n')
        po = subprocess.Popen(bsub(cmd, log_fp, 'metasoft'), stdout=fp, stderr=subprocess.STDOUT, shell=True)

def main(args):

    # 1. Get path to all gene files
         
    assert os.path.exists(args.chunk_fp)

    def chunk_id(fp):
       return(os.path.split(fp)[1].split('.')[0].split('_')[-1])

    ms_outdir = os.path.join(args.output_dir, 'metasoft')
    if not os.path.exists(ms_outdir):
        os.makedirs(ms_outdir)

    tmpdir = tempfile.mkdtemp(dir=args.output_dir)
    log_file_dir = os.path.join(args.output_dir, 'logs') 

    if not os.path.exists(log_file_dir):
        os.makedirs(log_file_dir)

    chunk_name = chunk_id(args.chunk_fp) 

    for k,chunk_df in enumerate(pd.read_csv(args.chunk_fp, sep='\t', chunksize=args.chunksize)): 
        print('\rprocessing {}'.format(k), end='', flush=True)

        mini_chunk_name = chunk_name + str(k+1)
        log_outfile = os.path.join(log_file_dir, mini_chunk_name +'.log')
    
        input_chunk_path = os.path.join(tmpdir, 'ms_chunk_'+mini_chunk_name+'.txt') 
        chunk_df.to_csv(input_chunk_path, sep='\t', header=False, index=False)
    
        run_metasoft(input_chunk_path, mini_chunk_name, ms_outdir, log_outfile)
    
    log.info('finished processing chunks.')

#shutil.rmtree(tmpdir)

if __name__ == '__main__':

    args = parse_args()
    main(args)
