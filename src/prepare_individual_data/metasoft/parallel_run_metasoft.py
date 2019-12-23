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
    parser.add_argument('chunk_dir', help='Directory METASOFT input chunks.')
    parser.add_argument('-s', '--chunk_blacklist', help='File listing input chunk IDs to NOT include.')
    parser.add_argument('-c', '--chunksize', type=int, default=3, help='Number of chunks to runper blade on MGI.')
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

def bsub(cmd, log, job_name, mem=30, gtmp=4, docker_image="apollodorus/brain-eqtl:eqtl-analysis", queue="research-hpc"):
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

def loop_chunks(l, n, d):

    # l: list, n: chunksize, d: defaultdict(lambda: dict())

    for num, i in enumerate(range(0, len(l), n)): 

        chunk_files = l[i:i+n]

        chunk_name = 'chunk_'+str(num+1)

        d[chunk_name] = {'return_code':None, 'genes':chunk_files, 'log':None}

        yield(chunk_name, chunk_files)

def run_metasoft(chunk_list, chunk_name, outdir, process_container, process_dict, log_fp):

    exc = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_metasoft.py')
    cmd = " ".join([exc, chunk_list, '-o', outdir])

    with open(log_fp, 'w+') as fp:
        fp.write('* cmd: '+cmd+'\n')
        po = subprocess.Popen(bsub(cmd, log_fp, 'metasoft'), stdout=fp, stderr=subprocess.STDOUT, shell=True)
        process_container.put((chunk_name, po))

def main(args):

    # 1. Get path to all gene files
         
    chunk_path_list= [os.path.join(os.path.abspath(args.chunk_dir),i) for i in os.listdir(args.chunk_dir) if i.startswith('ms_chunk')]
    if args.chunk_blacklist:
        with open(args.chunk_blacklist) as fp:
            bl = fp.read().strip().split('\n')

        def chunk_id(fp):
            return(os.path.split(fp)[1].split('.')[0].split('_')[-1])

        pl = len(chunk_path_list)
        chunk_path_list = [i for i in chunk_path_list if chunk_id(i) not in bl]
        print('{}/{} chunks filtered'.format(len(chunk_path_list), pl))

    ms_outdir = os.path.join(args.output_dir, 'metasoft')
    if not os.path.exists(ms_outdir):
        os.makedirs(ms_outdir)

    assert all([os.path.exists(i) for i in chunk_path_list])

    tmpdir = tempfile.mkdtemp(dir=args.output_dir)
    log_file_dir = os.path.join(args.output_dir, 'logs') 

    chunk_dict = defaultdict(lambda: dict())
    chunk_processes = queue.Queue()

    if not os.path.exists(log_file_dir):
        os.makedirs(log_file_dir)

    job_cycle = 0
    for k,(chunk_name, chunk_list) in enumerate(loop_chunks(chunk_path_list, args.chunksize, chunk_dict)):

        log_outfile = os.path.join(log_file_dir, chunk_name +'.log')
        chunk_dict[chunk_name]['log'] = log_outfile

        input_chunk_path = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)
        with open(input_chunk_path.name, 'w') as fp:
            fp.write('\n'.join(chunk_list))

        run_metasoft(input_chunk_path.name, chunk_name, ms_outdir, chunk_processes, chunk_dict, log_outfile)

    print('Submitted {} jobs.'.format(k+1))
    if __name__ == '__main__':

    args = parse_args()
    main(args)
