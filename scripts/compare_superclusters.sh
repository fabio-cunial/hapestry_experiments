#!/bin/bash
#
HAPESTRY_SUPERCLUSTERS_TSV="50bp_32x_hapestry_HG002_superclusters.tsv"
KANPIG_SUPERCLUSTERS_TSV="50bp_32x_kanpig_HG002_superclusters.tsv"

set -euxo pipefail


# Computing all overlaps between superclusters
cat ${HAPESTRY_SUPERCLUSTERS_TSV} | tail -n +2 | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
    if ($10<$11) min=$10; \
    else min=$11; \
    printf("%s\t%d\t%d\t%d\n",$1,$3,$4,min); \
}' > hapestry.bed
cat ${KANPIG_SUPERCLUSTERS_TSV} | tail -n +2 | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
    if ($10<$11) min=$10; \
    else min=$11; \
    printf("%s\t%d\t%d\t%d\n",$1,$3,$4,min); \
}' > kanpig.bed
N_RECORDS_HAPESTRY=$(wc -l < hapestry.bed)
N_RECORDS_KANPIG=$(wc -l < kanpig.bed)

# Building clusters of overlapping or adjacent records
java CompareSuperclusters hapestry.bed ${N_RECORDS_HAPESTRY} kanpig.bed ${N_RECORDS_KANPIG} > hapestry_kanpig.tsv
java CompareSuperclusters kanpig.bed ${N_RECORDS_KANPIG} hapestry.bed ${N_RECORDS_HAPESTRY} > tmp.tsv
cut -f 1-3 hapestry_kanpig.tsv > tmp1.txt
cut -f 1-3 tmp.tsv > tmp2.txt
diff --brief tmp1.txt tmp2.txt
rm -f tmp.tsv tmp1.txt tmp2.txt
cut -f 4,5 hapestry_kanpig.tsv > hapestry_kanpig_distances.tsv
rm -f hapestry_kanpig.tsv







#bedtools intersect -wo -a hapestry.bed -b kanpig.bed > hapestry_kanpig_overlaps.bed
#bedtools intersect -wo -a kanpig.bed -b hapestry.bed > kanpig_hapestry_overlaps.bed
#rm -f hapestry.bed kanpig.bed
# # Computing best matches between superclusters
# java GetBestMatches hapestry_kanpig_overlaps.bed > hapestry_kanpig_best.bed
# java GetBestMatches kanpig_hapestry_overlaps.bed > kanpig_hapestry_best.bed
# cut -f 1,2,3,5,6,7 hapestry_kanpig_best.bed > hk.bed
# awk 'BEGIN { IFS="\t"; OFS="\t"; } {print $5, $6, $7, $1, $2, $3}' kanpig_hapestry_best.bed > kh.bed
# diff --brief hk.bed kh.bed
# rm -f hk.bed kh.bed
# N_MATCHING_SUPERCLUSTERS_HAPESTRY=$(wc -l < hapestry_kanpig_best.bed)
# N_TOTAL_SUPERCLUSTERS_HAPESTRY=$(wc -l < ${HAPESTRY_SUPERCLUSTERS_TSV})
# echo "${N_MATCHING_SUPERCLUSTERS_HAPESTRY}/${N_TOTAL_SUPERCLUSTERS_HAPESTRY} superclusters in hapestry match those in kanpig"
# N_MATCHING_SUPERCLUSTERS_KANPIG=$(wc -l < kanpig_hapestry_best.bed)
# N_TOTAL_SUPERCLUSTERS_KANPIG=$(wc -l < ${KANPIG_SUPERCLUSTERS_TSV})
# echo "${N_MATCHING_SUPERCLUSTERS_KANPIG}/${N_TOTAL_SUPERCLUSTERS_KANPIG} superclusters in kanpig match those in hapestry"
#
# # Printing the edit distances of each pair of superclusters
# cut -f 4,8 hapestry_kanpig_best.bed > hapestry_kanpig_distances.tsv
