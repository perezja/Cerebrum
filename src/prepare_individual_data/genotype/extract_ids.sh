#!/bin/bash
OUT_DIR=$1
if [ $# -eq 0 ]
then
    echo -e "usage: extract_ids.sh output_dir\n\toutput_dir: output directory to write 'participant_id_list.txt'"; exit;
fi

for i in `find ../data/expression -name '*txt.gz'`;
do id=$(echo $i | cut -d'/' -f 4 | cut -d'.' -f 1); 
  mkdir -p $OUT_DIR;
  zcat $i | awk -F"\t" 'BEGIN {OFS="\n"} NR==1 {$1=""; print substr($0, 2)}' >> "$OUT_DIR/${id}_id_list.txt";
  echo "wrote to: $OUT_DIR/${id}_id_list.txt";
done;   

