#!/usr/bin/env python3
import argparse
import subprocess
import tempfile
import shutil
import json
import os

def parse_arguments():

    parser = argparse.ArgumentParser(description='Run METASOFT.')
    parser.add_argument('chunk_fp', help='Paths to METASOFT input file.')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    return(args)
     
def main(args):

    # chunk files should be gzipped
    prefix = os.path.split(args.chunk_fp)[1].split('.')[0].split('_')[-1]
    print('Running "{}"'.format(prefix))

    cmd = ""

    # jar and pvalue table paths are in the docker image
    cmd += 'java -jar /opt/metasoft/Metasoft.jar'\
        +' -input '+args.chunk_fp\
        +' -pvalue_table /opt/metasoft/HanEskinPvalueTable.txt' \
        +' -output '+os.path.join(args.output_dir, prefix+'.metasoft.txt')\
        +' -mvalue_p_thres 1.0'\
        +' -mvalue_method mcmc'\
        +' -seed 100'\
        +' -log '+os.path.join(args.output_dir, prefix+'.metasoft.log')\
        +' -mvalue'

    subprocess.check_call(cmd, shell=True)

    print('finished chunk processing.') 

if __name__ == '__main__':
    args = parse_arguments()
    main(args)
