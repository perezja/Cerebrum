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
    parser.add_argument('chunk_dir', help='Directory with matrix-eqtl assocation data.')
    parser.add_argument('prefix', help='Prefix for group names.')
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

def bsub(cmd, job_name, mem=20, gtmp=4, docker_image="apollodorus/smr:1.02", queue="research-hpc"):
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

def fragment_chunk_group(chunk_fp, chunk_name, gene_outdir, tmpdir, process_container, process_dict, log_fp):

    cmd = "awk -F\"\\t\" 'BEGIN{{OFS=\"\\t\";}} NR>1 {{print > \"{}/\"$2}}' {}".format(gene_outdir, chunk_fp)

    exc_file = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)

    with open(exc_file.name, 'w') as fp:
        fp.write('#!/bin/bash\n' + cmd)
        fp.close()

    os.chmod(exc_file.name, 0o755)

    cmd = exc_file.name

    with open(log_fp, 'w+') as fp:
        fp.write(cmd)
        po = subprocess.Popen(bsub(cmd, 'fragment_chunks'), stdout=fp, stderr=subprocess.STDOUT, shell=True)
        process_container.put((chunk_name, po))

def main(args):
         
    tmpdir = tempfile.mkdtemp(dir=args.output_dir)
    log_file_dir = os.path.join(args.output_dir, 'logs') 

    chunk_dict = defaultdict(lambda: dict())
    chunk_processes = queue.Queue()

    if not os.path.exists(log_file_dir):
        os.makedirs(log_file_dir)

    chunk_files = [os.path.join(os.path.abspath(args.chunk_dir),i) for i in os.listdir(args.chunk_dir) if i.endswith('.assoc.txt')]

    ## DEBUGG
#    chunk_files = [os.path.join(os.path.abspath(args.chunk_dir),i) for i in os.listdir(args.chunk_dir) if i.endswith('_debug.txt')]
    group_dir = os.path.join(args.output_dir, args.prefix)
    if not os.path.exists(group_dir):
        os.makedirs(group_dir)

    for k,chunk_fp in enumerate(chunk_files):

        chunk_name = os.path.split(chunk_fp)[1].split('.')[0]
        log_outfile = os.path.join(log_file_dir, chunk_name +'.log')
        chunk_dict[chunk_name]['log'] = log_outfile

        fragment_chunk_group(chunk_fp, chunk_name, group_dir, tmpdir, chunk_processes, chunk_dict, log_outfile)

    n = 0
    while not chunk_processes.empty():
        chunk_name, po = chunk_processes.get()

        while po.poll() is None:
            time.sleep(0.5)

        n+=1
        print('\r{}/{} jobs are returned.'.format(n, len(chunk_files)), end='', flush=True)
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
