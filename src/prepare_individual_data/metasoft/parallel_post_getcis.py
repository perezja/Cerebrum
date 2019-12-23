#!/usr/bin/env python3
# Author: James A. Perez 

import argparse
import numpy as np
import pandas as pd
from datetime import datetime
import subprocess
import glob
import re
import json
import os

parser = argparse.ArgumentParser(description='Prepare input for Metasoft.')
parser.add_argument('wave_list', help='List of wave directories made during METASOFT processing')
parser.add_argument('bed_file', help='bed file with TSS for genes.')
parser.add_argument('-p', '--pthreshold', default=1, type=float, help='')
parser.add_argument('-o', '--output_dir', default='./', help='Directory for output files.')
args = parser.parse_args() 

def get_chk_files(dn):

    def chk_name(fp):
        return(os.path.split(fp)[1].split('.')[0])
        

    chk_files = [(chk_name(chk_fp),chk_fp.replace('.log', '.txt')) for chk_fp in glob.glob(os.path.join(dn, 'metasoft', '*.log'))]
    assert all([os.path.exists(i[1]) for i in chk_files])

    return(chk_files)

def bsub(cmd, log, job_name, mem=2, gtmp=1, docker_image="apollodorus/bioinf:pr", queue="research-hpc"):
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

def post_process(chunk_fp, bed_file, chunk_name, outdir, log_fp):

    exc = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'post_process_getcis.py')
    cmd = " ".join([exc, chunk_fp, bed_file, chunk_name, '-o', outdir, '-p', str(args.pthreshold)])

    with open(log_fp, 'w+') as fp:
        fp.write('* cmd: '+cmd+'\n')
        po = subprocess.Popen(bsub(cmd, log_fp, 'pandas'), stdout=fp, stderr=subprocess.STDOUT, shell=True)

def main():
    
    with open(args.wave_list) as fp:
        l = fp.read().strip().split('\n')

    chk_list = list()
    for n,dn in enumerate(l):
        chk_list += get_chk_files(dn)

    print('detected {} total chunks.'.format(len(chk_list)))
    
    log_dir = os.path.join(args.output_dir, 'logs')
    outdir = os.path.join(args.output_dir, 'post_ms_chunks')

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    with open(os.path.join(args.output_dir, 'chunk_detected.tmp'), 'w+') as fp:
        fp.write('\n'.join(['\t'.join(i) for i in chk_list]))

    for k,(chk_name,fp) in enumerate(chk_list):

        chk_name=chk_name+'_chunk'+str(k)
        log_fp = os.path.join(log_dir, chk_name+'.log')
        post_process(fp, args.bed_file, chk_name, outdir, log_fp) 

    print('launched all jobs.')

if __name__ == "__main__":
    main()

