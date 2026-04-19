#!/bin/bash
#
ROOT_DIR="/Users/fcunial/Downloads/HAPESTRY/experimental_section_1/experimental_section_1_1/grapheval_per_sample/50bp_32x_HG002/chr1_evaluation"

set -euo pipefail

#find "${ROOT_DIR}" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; > dirs.txt
while read DIRECTORY; do
    AVG_COVERAGE_HAPESTRY=$(cat ${ROOT_DIR}/${DIRECTORY}/hapestry/nodes.csv | awk ' \
        BEGIN {  FS=","; OFS=","; numerator=0; denominator=0; } { \
            if ($3==0) { numerator+=$5; denominator+=1; } \
        } \
        END { \
            if (denominator>0) printf("%f\n",numerator/denominator);
            else printf("0\n"); \
        }')
    AVG_IDENTITY_HAPESTRY=$(cat ${ROOT_DIR}/${DIRECTORY}/hapestry/nodes.csv | awk ' \
        BEGIN {  FS=","; OFS=","; numerator=0; denominator=0; } { \
            if ($3==0) { numerator+=$6; denominator+=1; } \
        } \
        END { \
            if (denominator>0) printf("%f\n",numerator/denominator);
            else printf("0\n"); \
        }')

    AVG_COVERAGE_KANPIG=$(cat ${ROOT_DIR}/${DIRECTORY}/truvari_1_kanpig/nodes.csv | awk ' \
        BEGIN {  FS=","; OFS=","; numerator=0; denominator=0; } { \
            if ($3==0) { numerator+=$5; denominator+=1; } \
        } \
        END { \
            if (denominator>0) printf("%f\n",numerator/denominator);
            else printf("0\n"); \
        }')
    AVG_IDENTITY_KANPIG=$(cat ${ROOT_DIR}/${DIRECTORY}/truvari_1_kanpig/nodes.csv | awk ' \
        BEGIN {  FS=","; OFS=","; numerator=0; denominator=0; } { \
            if ($3==0) { numerator+=$6; denominator+=1; } \
        } \
        END { \
            if (denominator>0) printf("%f\n",numerator/denominator);
            else printf("0\n"); \
        }')

    if (( $(bc -l <<< "${AVG_COVERAGE_HAPESTRY} < ${AVG_COVERAGE_KANPIG}") )); then
        echo "${DIRECTORY} COVERAGE HAPESTRY: ${AVG_COVERAGE_HAPESTRY} KANPIG: ${AVG_COVERAGE_KANPIG}"
    fi
    if (( $(bc -l <<< "${AVG_IDENTITY_HAPESTRY} < ${AVG_IDENTITY_KANPIG}") )); then
        echo "${DIRECTORY} IDENTITY HAPESTRY: ${AVG_IDENTITY_HAPESTRY} KANPIG: ${AVG_IDENTITY_KANPIG}"
    fi
done < dirs.txt
#rm -f dirs.txt
