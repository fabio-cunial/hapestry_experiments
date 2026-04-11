#!/bin/bash
#
SAMPLES="HG002 HG00438 HG00733 HG01123 HG01891 HG02080 HG02572 HG03098 HG03492 HG03516"
REMOTE_ROOT="gs://fc-a90ab401-9c4b-43d1-b891-f0410c667ff2/submissions/"
LOCAL_ROOT=""  # Absolute path

set -euxo pipefail


# -------------------------- Steps of the pipeline -----------------------------

function DownloadSuperclusters() {
    CURRENT_DIR="${LOCAL_ROOT}/10bp/8x/kanpig/"
    rm -rf ${CURRENT_DIR} ; mkdir -p ${CURRENT_DIR} ; cd ${CURRENT_DIR}
    for SUBMISSION_ID in 1b4c026f-b753-40da-9d17-401b342008f6  60fcf740-4ed8-4fd6-aebc-b5288fd77001  ea3d9c0c-d121-4ce2-941d-982d47d57375 ; do
        gcloud storage cp ${REMOTE_ROOT}/${SUBMISSION_ID}/VcfdistEvaluationPrime/'*/call-Impl/*_superclusters.tsv' . --gzip-encoded
    done

    CURRENT_DIR="${LOCAL_ROOT}/10bp/32x/kanpig/"
    rm -rf ${CURRENT_DIR} ; mkdir -p ${CURRENT_DIR} ; cd ${CURRENT_DIR}
    for SUBMISSION_ID in c3bf0950-d2af-4c2b-8101-8acedd8863f3  9508d7a6-38b2-4e95-9906-cc263e853e0a  88fc89ae-2b58-4981-a6d5-ed1ac590f997 ; do
        gcloud storage cp ${REMOTE_ROOT}/${SUBMISSION_ID}/VcfdistEvaluationPrime/'*/call-Impl/*_superclusters.tsv' .
    done

    CURRENT_DIR="${LOCAL_ROOT}/50bp/8x/kanpig/"
    rm -rf ${CURRENT_DIR} ; mkdir -p ${CURRENT_DIR} ; cd ${CURRENT_DIR}
    for SUBMISSION_ID in 324acd46-5790-428b-81b6-af2a755817eb  4cb13e3b-76dd-4a89-aa95-9f3ad68b95c8  fad42427-1eb3-4967-aea9-4d2290abbf3f ; do
        gcloud storage cp ${REMOTE_ROOT}/${SUBMISSION_ID}/VcfdistEvaluationPrime/'*/call-Impl/*_superclusters.tsv' .
    done

    CURRENT_DIR="${LOCAL_ROOT}/50bp/32x/kanpig/"
    rm -rf ${CURRENT_DIR} ; mkdir -p ${CURRENT_DIR} ; cd ${CURRENT_DIR}
    for SUBMISSION_ID in 2ebc8a6b-931a-4c83-be6d-caf84b0afac9  06f47761-9a5b-4326-bf8e-83da8182f770  99dd40fa-e637-4f2c-9ef7-d39dd1df485c ; do
        gcloud storage cp ${REMOTE_ROOT}/${SUBMISSION_ID}/VcfdistEvaluationPrime/'*/call-Impl/*_superclusters.tsv' .
    done

    CURRENT_DIR="${LOCAL_ROOT}/10bp/8x/hapestry/"
    rm -rf ${CURRENT_DIR} ; mkdir -p ${CURRENT_DIR} ; cd ${CURRENT_DIR}
    for SUBMISSION_ID in 2ef5c92e-8000-4cd8-ae12-8f956149d451  7188be2b-a9ab-49ee-a4cf-d4e24803087d  519269d2-4e06-4cdf-a09f-07544c9b47f4 ; do
        gcloud storage cp ${REMOTE_ROOT}/${SUBMISSION_ID}/VcfdistEvaluationPrime/'*/call-Impl/*_superclusters.tsv' .
    done

    CURRENT_DIR="${LOCAL_ROOT}/10bp/32x/hapestry/"
    rm -rf ${CURRENT_DIR} ; mkdir -p ${CURRENT_DIR} ; cd ${CURRENT_DIR}
    for SUBMISSION_ID in 1602f71d-51d4-4ffb-b947-043bb207f364  22df22ef-c20b-4a74-8e7f-6439e9c1f0e8  5b2bbff3-4f03-4371-aed5-680935c12355 ; do
        gcloud storage cp ${REMOTE_ROOT}/${SUBMISSION_ID}/VcfdistEvaluationPrime/'*/call-Impl/*_superclusters.tsv' .
    done

    CURRENT_DIR="${LOCAL_ROOT}/50bp/8x/hapestry/"
    rm -rf ${CURRENT_DIR} ; mkdir -p ${CURRENT_DIR} ; cd ${CURRENT_DIR}
    for SUBMISSION_ID in 5ec3fe6b-8afc-4572-aa88-2e525481d44a  cc3b1878-54d7-4e47-88e6-c78a4c15b670 ; do
        gcloud storage cp ${REMOTE_ROOT}/${SUBMISSION_ID}/VcfdistEvaluationPrime/'*/call-Impl/*_superclusters.tsv' .
    done

    CURRENT_DIR="${LOCAL_ROOT}/50bp/32x/hapestry/"
    rm -rf ${CURRENT_DIR} ; mkdir -p ${CURRENT_DIR} ; cd ${CURRENT_DIR}
    for SUBMISSION_ID in b83ef593-f15e-4ecb-a22f-0abb002f3a9e  915ebb16-dadd-4790-b530-6b641f128fa5 ; do
        gcloud storage cp ${REMOTE_ROOT}/${SUBMISSION_ID}/VcfdistEvaluationPrime/'*/call-Impl/*_superclusters.tsv' .
    done
    
    cd ${LOCAL_ROOT}
}


function CompareSuperclusters() {
    local SAMPLE_ID=$1
    local MIN_SV_LENGTH=$2
    local COVERAGE=$3
    local FRACTIONS_FILE=$4  # The procedure appends to it
    
    HAPESTRY_SUPERCLUSTERS_TSV="${LOCAL_ROOT}/${MIN_SV_LENGTH}/${COVERAGE}/hapestry/${SAMPLE_ID}_superclusters.tsv"
    KANPIG_SUPERCLUSTERS_TSV="${LOCAL_ROOT}/${MIN_SV_LENGTH}/${COVERAGE}/kanpig/${SAMPLE_ID}_superclusters.tsv"
    
    # Extracting supercluster intervals and distances
    cat ${HAPESTRY_SUPERCLUSTERS_TSV} | tail -n +2 | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
        if ($10<$11) min=$10; \
        else min=$11; \
        printf("%s\t%d\t%d\t%d\n",$1,$3,$4,min); \
    }' > ${SAMPLE_ID}_hapestry.bed
    cat ${KANPIG_SUPERCLUSTERS_TSV} | tail -n +2 | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
        if ($10<$11) min=$10; \
        else min=$11; \
        printf("%s\t%d\t%d\t%d\n",$1,$3,$4,min); \
    }' > ${SAMPLE_ID}_kanpig.bed
    N_SUPERCLUSTERS_HAPESTRY=$(wc -l < ${SAMPLE_ID}_hapestry.bed)
    N_SUPERCLUSTERS_KANPIG=$(wc -l < ${SAMPLE_ID}_kanpig.bed)

    # Building clusters of overlapping or adjacent records
    java CompareSuperclusters ${SAMPLE_ID}_hapestry.bed ${N_SUPERCLUSTERS_HAPESTRY} ${SAMPLE_ID}_kanpig.bed ${N_SUPERCLUSTERS_KANPIG} 1> ${SAMPLE_ID}_hapestry_kanpig.tsv 2>> ${FRACTIONS_FILE}

    # Basic consistency check
    java CompareSuperclusters ${SAMPLE_ID}_kanpig.bed ${N_SUPERCLUSTERS_KANPIG} ${SAMPLE_ID}_hapestry.bed ${N_SUPERCLUSTERS_HAPESTRY} > tmp.tsv
    cut -f 1-3 hapestry_kanpig.tsv > tmp1.txt
    cut -f 1-3 tmp.tsv > tmp2.txt
    diff --brief tmp1.txt tmp2.txt
    rm -f tmp.tsv tmp1.txt tmp2.txt

    # Extracting pairs of distances for plotting
    cut -f 4,5 ${SAMPLE_ID}_hapestry_kanpig.tsv > ${SAMPLE_ID}_hapestry_kanpig_distances.tsv
}



# ------------------------------- Main program ---------------------------------

DownloadSuperclusters 

for SAMPLE_ID in ${SAMPLES}; do
    FRACTIONS_FILE="${LOCAL_ROOT}/10bp/8x/fractions.txt"
    rm -f ${FRACTIONS_FILE}
    CompareSuperclusters ${SAMPLE_ID} 10bp 8x ${FRACTIONS_FILE}
done

for SAMPLE_ID in ${SAMPLES}; do
    FRACTIONS_FILE="${LOCAL_ROOT}/10bp/32x/fractions.txt"
    rm -f ${FRACTIONS_FILE}
    CompareSuperclusters ${SAMPLE_ID} 10bp 32x ${FRACTIONS_FILE}
done

for SAMPLE_ID in ${SAMPLES}; do
    FRACTIONS_FILE="${LOCAL_ROOT}/50bp/8x/fractions.txt"
    rm -f ${FRACTIONS_FILE}
    CompareSuperclusters ${SAMPLE_ID} 50bp 8x ${FRACTIONS_FILE}
done

for SAMPLE_ID in ${SAMPLES}; do
    FRACTIONS_FILE="${LOCAL_ROOT}/50bp/32x/fractions.txt"
    rm -f ${FRACTIONS_FILE}
    CompareSuperclusters ${SAMPLE_ID} 50bp 32x ${FRACTIONS_FILE}
done
