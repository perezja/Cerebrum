#!/usr/bin/env python3
import argparse
import subprocess
import tempfile
import shutil
import json
import os

def parse_arguments():

    parser = argparse.ArgumentParser(description='Run METASOFT.')
    parser.add_argument('chunk_file_list', help='File listing paths to METASOFT input chunks')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    return(args)
     
def main(args):

    tmpdir = tempfile.mkdtemp(dir=args.output_dir)
    with open(args.chunk_file_list, 'r') as fp:
        chunk_files = fp.read().strip().split('\n') 

    # chunk files should be gzipped
    for chunk_path in chunk_files: 

        prefix = os.path.split(chunk_path)[1].split('.')[0].split('_')[-1]

        cmd = ""

        input_file = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False) 
        cmd += 'gunzip -c {} | tail -n+2 > {}; '.format(chunk_path, input_file.name)

        # jar and pvalue table paths are in the docker image
        cmd += 'java -jar /opt/metasoft/Metasoft.jar'\
            +' -input '+input_file.name\
            +' -pvalue_table /opt/metasoft/HanEskinPvalueTable.txt' \
            +' -output '+os.path.join(args.output_dir, prefix+'.metasoft.txt')\
            +' -mvalue_p_thres 1.0'\
            +' -mvalue_method mcmc'\
            +' -seed 100'\
            +' -log '+os.path.join(args.output_dir, prefix+'.metasoft.log')\
            +' -mvalue'

        subprocess.check_call(cmd, shell=True)
        os.remove(input_file.name)

    print('finished chunk processing.') 
    shutil.rmtree(tmpdir)

if __name__ == '__main__':
    args = parse_arguments()
    main(args)
