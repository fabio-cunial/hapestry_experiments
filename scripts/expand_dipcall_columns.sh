#!/bin/bash
#
# Expands columns $02_dipcall_normed_filtered_{tbi,vcf}$, that contain arrays of
# strings, into separate columns (with one string per column).
#
INPUT_TSV="hprc_y1_sample.tsv"
COLUMN_FROM="12"  # 02_dipcall_normed_filtered_tbi
COLUMN_TO="13"  # 02_dipcall_normed_filtered_vcf

# Expanding the header
head -n 1 ${INPUT_TSV} > header.txt
cut -f 1-$((${COLUMN_FROM} - 1)) header.txt > left.txt
echo -e "02_dipcall_normed_filtered_10_tbi\t02_dipcall_normed_filtered_20_tbi\t02_dipcall_normed_filtered_30_tbi\t02_dipcall_normed_filtered_40_tbi\t02_dipcall_normed_filtered_50_tbi\t02_dipcall_normed_filtered_10_vcf\t02_dipcall_normed_filtered_20_vcf\t02_dipcall_normed_filtered_30_vcf\t02_dipcall_normed_filtered_40_vcf\t02_dipcall_normed_filtered_50_vcf" > center.txt
cut -f $((${COLUMN_TO} + 1))- header.txt > right.txt
paste left.txt center.txt right.txt > header_new.txt
rm -f left.txt center.txt right.txt

# Expanding the body
tail -n +2 ${INPUT_TSV} > body.txt
cut -f 1-$((${COLUMN_FROM} - 1)) body.txt > left.txt
cut -f ${COLUMN_FROM}-${COLUMN_TO} body.txt | tr -d '[' | tr -d ']' | tr ',' '\t' > center.txt
cut -f $((${COLUMN_TO} + 1))- body.txt > right.txt
paste left.txt center.txt right.txt > body_new.txt
rm -f left.txt center.txt right.txt

# Concatenating
cat header_new.txt body_new.txt > out.tsv
rm -f header_new.txt body_new.txt
