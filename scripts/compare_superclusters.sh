#!/bin/bash
#
HAPESTRY_SUPERCLUSTERS_TSV="/Users/fcunial/Downloads/50bp_32x_hapestry_HG002_superclusters.tsv"
KANPIG_SUPERCLUSTERS_TSV="/Users/fcunial/Downloads/50bp_32x_kanpig_HG002_superclusters.tsv"


set -euxo pipefail


# Extracting supercluster intervals and distances
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
N_SUPERCLUSTERS_HAPESTRY=$(wc -l < hapestry.bed)
N_SUPERCLUSTERS_KANPIG=$(wc -l < kanpig.bed)

# Building clusters of overlapping or adjacent records
java CompareSuperclusters hapestry.bed ${N_SUPERCLUSTERS_HAPESTRY} kanpig.bed ${N_SUPERCLUSTERS_KANPIG} > hapestry_kanpig.tsv

# Basic consistency check
java CompareSuperclusters kanpig.bed ${N_SUPERCLUSTERS_KANPIG} hapestry.bed ${N_SUPERCLUSTERS_HAPESTRY} > tmp.tsv
cut -f 1-3 hapestry_kanpig.tsv > tmp1.txt
cut -f 1-3 tmp.tsv > tmp2.txt
diff --brief tmp1.txt tmp2.txt
rm -f tmp.tsv tmp1.txt tmp2.txt

# Extracting pairs of distances for plotting
cut -f 4,5 hapestry_kanpig.tsv > hapestry_kanpig_distances.tsv
