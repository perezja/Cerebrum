def merge_chunks(chunk_dir_loc, header, outfile):
    """ process to concatenate all finished chunk output files."""

    with open(outfile, mode='wt') as fp:
        fp.write('\t'.join(header)+'\n')

    # find: full path to all files in specified directory
    # tail: cat all but header line
    cmd = "find {} -maxdepth 1 -type f | xargs -I{{}} tail -n+2 {{}} | cat >> {}".format(chunk_dir_loc, outfile) 

    subprocess.check_call(bsub(cmd, 'merging'), shell=True)
