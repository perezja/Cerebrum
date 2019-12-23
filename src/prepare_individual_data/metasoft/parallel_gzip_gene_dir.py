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
    parser.add_argument('gene_chunk_dir', help='Directory with study-specific gene chunk directories.')
    parser.add_argument('-c', '--chunksize', type=int, default=75, help='Number of genes to merge per blade on MGI.')
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

def bsub(cmd, job_name, mem=20, gtmp=4, docker_image="apollodorus/bioinf:pr", queue="research-hpc"):
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
    bsub = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -Is -q "+queue \
        + " -J "+job_name \
        + " -M "+mem+"000000" \
        + " -R 'select[mem>"+mem+"000 && gtmp > "+gtmp+"]" \
        + " rusage[mem="+mem+"000, gtmp="+gtmp+"]'" \
        + " -a 'docker("+docker_image+")'" \
        + " /bin/bash -c '"+cmd+"'" 

    
    return(bsub)

def loop_chunks(l, n, d):

    # l: list, n: chunksize, d: defaultdict(lambda: dict())

    for num, i in enumerate(range(0, len(l), n)): 

        chunk_files = l[i:i+n]

        chunk_name = 'chunk_'+str(num+1)

        d[chunk_name] = {'return_code':None, 'genes':chunk_files, 'log':None}

        yield(chunk_name, chunk_files)

def gzip_gene_chunk(gene_list, chunk_name, process_container, process_dict, log_fp):

    exc = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gzip_dir.py')
    cmd = " ".join([exc, gene_list])

    with open(log_fp, 'w+') as fp:
        po = subprocess.Popen(bsub(cmd, 'gzip'), stdout=fp, stderr=subprocess.STDOUT, shell=True)
        process_container.put((chunk_name, po))

def main(args):

    # 1. Get path to all gene files
         
    gene_path_list= [os.path.join(os.path.abspath(args.gene_chunk_dir),i) for i in os.listdir(args.gene_chunk_dir) if i.startswith('ENSG')]

    assert all([os.path.exists(i) for i in gene_path_list])

    # 2. Run parallelization of gene-wise gzipping 

    tmpdir = tempfile.mkdtemp(dir=args.output_dir)
    log_file_dir = os.path.join(args.output_dir, 'logs') 

    chunk_dict = defaultdict(lambda: dict())
    chunk_processes = queue.Queue()

    if not os.path.exists(log_file_dir):
        os.makedirs(log_file_dir)

    for k,(chunk_name, gene_list) in enumerate(loop_chunks(gene_path_list, args.chunksize, chunk_dict)):

        log_outfile = os.path.join(log_file_dir, chunk_name +'.log')
        chunk_dict[chunk_name]['log'] = log_outfile

        gl_path = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)
        with open(gl_path.name, 'w') as fp:
            fp.write('\n'.join(gene_list))

        gzip_gene_chunk(gl_path.name, chunk_name, chunk_processes, chunk_dict, log_outfile)

    n = 0
    while not chunk_processes.empty():
        chunk_name, po = chunk_processes.get()

        while po.poll() is None:
            time.sleep(0.5)

        n+=1
        print('\r{}/{} jobs are returned.'.format(n, k+1, end='', flush=True))
        chunk_dict[chunk_name]['return_code'] = po.returncode

        # job failed
        if po.returncode:
            with open(chunk_dict[chunk_name]['log'], 'r') as fp:
                log.error(fp.read())
                log.error("{} failed with exit code {}".format(chunk_name, po.returncode))

    monitor_json = os.path.join(args.output_dir, 'monitor.json')

    with open(monitor_json, 'w+') as fp:
        json.dump(chunk_dict, fp)

    log.info('wrote log to: {}'.format(monitor_json))
    
    log.info('finished processing chunks.')
    
    shutil.rmtree(tmpdir)

if __name__ == '__main__':

    args = parse_args()
    main(args)
