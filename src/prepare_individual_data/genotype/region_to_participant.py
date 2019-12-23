import pandas as pd
import numpy as np
import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Extract participant ids from expression tables.')
parser.add_argument('expression_dir', help='Lookup table linking all participants to a given region.')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

region_files = [i for i in os.listdir(args.expression_dir) if i.endswith('.txt.gz')] 
region_ids = [ i.split('.')[0] for i in region_files]

region_files.sort(key=lambda x: os.path.split(x)[1])
region_ids.sort()

rs_df = list() 
for rid, f in zip(region_ids, region_files):
    print("\r* Extracting pariticipant ids for {}".format(rid), end='', flush=True)
    exp_file = os.path.join(args.expression_dir, f)
    exp_df = pd.read_csv(exp_file, sep='\t', compression='gzip', index_col=0)
    region_samples = exp_df.columns.tolist()
    
    print([x for n,x in enumerate(region_samples) if x in region_samples[:n]]) 
    df = pd.DataFrame({'participant_id':region_samples, 'region_id':np.full(len(region_samples), rid)})
    rs_df.append(df)
rs_df = pd.concat(rs_df)

rs_df.to_csv('region_to_participant_lookup.txt', index=False, sep='\t')
