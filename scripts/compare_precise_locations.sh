#!/bin/bash
#
# A version of `get_hap_consensus_alignment.sh` (from the hapestry repo) that
# picks the lowest scoring haplotypes of an SV merging tool from a graph
# evaluation run, reconstructs the corresponding haplotypes for several SV
# merging tools using their output VCFs, aligns them to the reference, and
# creates an IGV script to visualize the region. 
#
# This is useful e.g. for comparing kanpig to hapestry in terms of the precise
# location of the calls.
#
INPUT_DIR="/Users/fcunial/Downloads/trgt_track/32x"
OUTPUT_DIR="/Users/fcunial/Downloads/IGV_WINDOW"
REFERENCE_FA="/Users/fcunial/Downloads/trgt_track/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
TANDEMS_BED="/Users/fcunial/Downloads/trgt_track/human_GRCh38_no_alt_analysis_set.trf.bed"
GENES_BED="/Users/fcunial/Downloads/trgt_track/refseq_all_genes.bed.gz"

SLACK_BP="100000"
N_MINIMAP_THREADS="8"

BCFTOOLS_COMMAND="/Users/fcunial/git/bcftools-1.22/bcftools"
SAMTOOLS_COMMAND="/Users/fcunial/git/samtools-1.22/samtools"
export BCFTOOLS_PLUGINS="/Users/fcunial/git/bcftools-1.22/plugins"

SHOW_BCFTOOLS="1"
SHOW_TRUVARI="1"

set -euxo pipefail




# ------------------------------- Functions ------------------------------------

function getKanpigVcf() {
    local WINDOW_DOTS=$1
    local WINDOW_UNDERSCORE=$2
    local SAMPLE=$3
    
    ${BCFTOOLS_COMMAND} view -s ${SAMPLE} --regions-file regions.txt --regions-overlap variant ${INPUT_DIR}/03_truvari_mode1_kanpig.vcf.gz ${WINDOW_DOTS} | ${BCFTOOLS_COMMAND} view --include 'GT="alt" && ALT!="*"' > ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/kanpig_tmp.vcf
    ${BCFTOOLS_COMMAND} view --header-only ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/kanpig_tmp.vcf > ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_kanpig.vcf
    ${BCFTOOLS_COMMAND} view --no-header ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/kanpig_tmp.vcf | awk '{ p=index($10,":"); printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,"GT",substr($10,1,p-1)); }' >> ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_kanpig.vcf
    bgzip ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_kanpig.vcf
    rm -f ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/kanpig_tmp.vcf
}


function getConsensus() {
    local CHR=$1
    local WINDOW_UNDERSCORE=$2
    local WINDOW_PADDED=$3
    local SAMPLE=$4
    local TOOL_ID=$5
    
    # The following do not seem to be needed
    #python3 ./resolve_light.py ${OUTPUT_DIR}/${SAMPLE}/${TOOL_ID}.vcf.gz ${REFERENCE_FA} | bgzip > ${OUTPUT_DIR}/${SAMPLE}/${TOOL_ID}_resolved.vcf.gz
    #tabix -f ${OUTPUT_DIR}/${SAMPLE}/${TOOL_ID}_resolved.vcf.gz
    #${BCFTOOLS_COMMAND} +fill-from-fasta --output-type z ${OUTPUT_DIR}/${SAMPLE}/${TOOL_ID}_resolved.vcf.gz -- -c REF -f ${REFERENCE_FA} > ${OUTPUT_DIR}/${SAMPLE}/${TOOL_ID}_filled.vcf.gz
    #tabix -f ${OUTPUT_DIR}/${SAMPLE}/${TOOL_ID}_filled.vcf.gz
    
    # There is a bug in v1.22. Using an older version.
    bcftools consensus --haplotype 1 --fasta-ref ${CHR}.fa ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_${TOOL_ID}.vcf.gz > ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${TOOL_ID}_hap1.consensus &
    bcftools consensus --haplotype 2 --fasta-ref ${CHR}.fa ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_${TOOL_ID}.vcf.gz > ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${TOOL_ID}_hap2.consensus &
    wait
    ${SAMTOOLS_COMMAND} faidx ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${TOOL_ID}_hap1.consensus ${WINDOW_PADDED} > ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${TOOL_ID}_hap1.fa &
    ${SAMTOOLS_COMMAND} faidx ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${TOOL_ID}_hap2.consensus ${WINDOW_PADDED} > ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${TOOL_ID}_hap2.fa &
    wait
    rm -f ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${TOOL_ID}_hap*.consensus*
}


function printIgvScript() {
    local INPUT_DIR=$1
    local SAMPLE=$2
    local WINDOW=$3
    
    rm -rf ${INPUT_DIR}/igv.batch
    echo "new" >> ${INPUT_DIR}/igv.batch
    echo "genome hg38" >> ${INPUT_DIR}/igv.batch
    echo "goto ${WINDOW}" >> ${INPUT_DIR}/igv.batch
    echo "snapshotDirectory ${INPUT_DIR}" >> ${INPUT_DIR}/igv.batch
    
    # Top part
    echo "load ${GENES_BED}" >> ${INPUT_DIR}/igv.batch
    echo "load ${TANDEMS_BED}" >> ${INPUT_DIR}/igv.batch
    echo "load ${INPUT_DIR}/${SAMPLE}_dipcall.bam" >> ${INPUT_DIR}/igv.batch
    echo "load ${INPUT_DIR}/${SAMPLE}_hapestry.bam" >> ${INPUT_DIR}/igv.batch
    echo "load ${INPUT_DIR}/${SAMPLE}_kanpig.bam" >> ${INPUT_DIR}/igv.batch
    
    # Bottom part
    echo "load ${INPUT_DIR}/${SAMPLE}_dipcall.vcf.gz" >> ${INPUT_DIR}/igv.batch
    echo "load ${INPUT_DIR}/${SAMPLE}_hapestry.vcf.gz" >> ${INPUT_DIR}/igv.batch
    echo "load ${INPUT_DIR}/${SAMPLE}_kanpig.vcf.gz" >> ${INPUT_DIR}/igv.batch
    if [ ${SHOW_TRUVARI} -eq 1 ]; then
        echo "load ${INPUT_DIR}/${SAMPLE}_truvari_mode1.vcf.gz" >> ${INPUT_DIR}/igv.batch
    fi
    if [ ${SHOW_BCFTOOLS} -eq 1 ]; then
        echo "load ${INPUT_DIR}/${SAMPLE}_bcftools.vcf.gz" >> ${INPUT_DIR}/igv.batch
    fi
    
    echo "snapshot" >> ${INPUT_DIR}/igv.batch
}




# ------------------------------ Main program ----------------------------------

# Splitting and indexing the reference, if needed.
if [ -e chr1.fa ]; then 
    :
else
    for CHR in $(seq 1 22) X Y; do
        (${SAMTOOLS_COMMAND} faidx ${REFERENCE_FA} chr${CHR} > chr${CHR}.fa; minimap2 -t 1 -x asm20 -d chr${CHR}.mmi chr${CHR}.fa) &
    done
    wait
fi

# For every window, processing the top-3 samples with lowest kanpig identity.
rm -rf ${OUTPUT_DIR}/*
for CHROMOSOME_DIR in $(find ${INPUT_DIR} -type d -name "*_evaluation" -mindepth 1 -maxdepth 1 | sort --version-sort); do
    for WINDOW_DIR in $(find ${CHROMOSOME_DIR} -type d -mindepth 1 -maxdepth 1 | sort --version-sort); do
        WINDOW_UNDERSCORE=$(basename ${WINDOW_DIR})
        rm -rf ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}
        mkdir -p ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}
        WINDOW_DOTS=$(basename ${WINDOW_DIR} | tr '_' ':')
        echo ${WINDOW_DOTS} | tr ':' '\t' | tr '-' '\t' > regions.txt
        CHR=$(cut -f 1 regions.txt)
        FROM=$(cut -f 2 regions.txt)
        TO=$(cut -f 3 regions.txt)
        WINDOW_PADDED="${CHR}:$(( ${FROM} - ${SLACK_BP} ))-$(( ${TO} + ${SLACK_BP} ))"
        LOW_IDENTITY_SAMPLES=$(tail -n +2 ${WINDOW_DIR}/truvari_kanpig_50bp_32x/haps.csv | sort --numeric-sort -t , -k 6 | head -n 3 | cut -d , -f 1 | awk '{ p=index($0,"#"); printf("%s\n",substr($0,1,p-1)); }' | sort | uniq)
        for SAMPLE in ${LOW_IDENTITY_SAMPLES}; do
            rm -rf ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}
            mkdir -p ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}
            
            # Creating the VCFs
            ${BCFTOOLS_COMMAND} view -s ${SAMPLE} --regions-file regions.txt --regions-overlap variant ${INPUT_DIR}/../dipcall_50.vcf.gz ${WINDOW_DOTS} | ${BCFTOOLS_COMMAND} view --include 'GT="alt" && ALT!="*"' --output-type z > ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_dipcall.vcf.gz &
            if [ ${SHOW_BCFTOOLS} -eq 1 ]; then
                # bcftools of all samples
                ${BCFTOOLS_COMMAND} view --regions-file regions.txt --regions-overlap variant --include 'SVLEN<=10000 && SVLEN>=-10000' ${INPUT_DIR}/03_bcftools_merge.vcf.gz ${WINDOW_DOTS} | ${BCFTOOLS_COMMAND} view --include 'GT="alt" && ALT!="*"' --output-type z > ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_bcftools.vcf.gz &
            fi
            if [ ${SHOW_TRUVARI} -eq 1 ]; then
                ${BCFTOOLS_COMMAND} view -s ${SAMPLE} --regions-file regions.txt --regions-overlap variant --include 'SVLEN<=10000 && SVLEN>=-10000' ${INPUT_DIR}/03_truvari_mode1.vcf.gz ${WINDOW_DOTS} | ${BCFTOOLS_COMMAND} view --include 'GT="alt" && ALT!="*"' --output-type z > ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_truvari_mode1.vcf.gz &
            fi
            ${BCFTOOLS_COMMAND} view -s ${SAMPLE} --regions-file regions.txt --regions-overlap variant ${INPUT_DIR}/hapestry_on_windows.vcf.gz ${WINDOW_DOTS} | ${BCFTOOLS_COMMAND} view --include 'GT="alt" && ALT!="*"' --output-type z > ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_hapestry.vcf.gz &
            getKanpigVcf ${WINDOW_DOTS} ${WINDOW_UNDERSCORE} ${SAMPLE} &            
            wait
            tabix -f ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_dipcall.vcf.gz &
            if [ ${SHOW_BCFTOOLS} -eq 1 ]; then
                tabix -f ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_bcftools.vcf.gz &
            fi
            if [ ${SHOW_TRUVARI} -eq 1 ]; then
                tabix -f ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_truvari_mode1.vcf.gz &
            fi
            tabix -f ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_kanpig.vcf.gz &
            tabix -f ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_hapestry.vcf.gz &
            wait

            # Using the phased VCFs to create aligned haplotypes
            getConsensus ${CHR} ${WINDOW_UNDERSCORE} ${WINDOW_PADDED} ${SAMPLE} dipcall &
            getConsensus ${CHR} ${WINDOW_UNDERSCORE} ${WINDOW_PADDED} ${SAMPLE} kanpig &
            getConsensus ${CHR} ${WINDOW_UNDERSCORE} ${WINDOW_PADDED} ${SAMPLE} hapestry &
            wait
            minimap2 -t ${N_MINIMAP_THREADS} -a -L ${CHR}.fa ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/dipcall_hap1.fa ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/dipcall_hap2.fa | ${SAMTOOLS_COMMAND} sort --threads ${N_MINIMAP_THREADS} -O bam > ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_dipcall.bam &
            minimap2 -t ${N_MINIMAP_THREADS} -a -L ${CHR}.fa ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/kanpig_hap1.fa ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/kanpig_hap2.fa | ${SAMTOOLS_COMMAND} sort --threads ${N_MINIMAP_THREADS} -O bam > ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_kanpig.bam &
            minimap2 -t ${N_MINIMAP_THREADS} -a -L ${CHR}.fa ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/hapestry_hap1.fa ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/hapestry_hap2.fa | ${SAMTOOLS_COMMAND} sort --threads ${N_MINIMAP_THREADS} -O bam > ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_hapestry.bam &
            wait
            rm -f ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/*_hap*.fa
            ${SAMTOOLS_COMMAND} index -@ ${N_MINIMAP_THREADS} ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_dipcall.bam &
            ${SAMTOOLS_COMMAND} index -@ ${N_MINIMAP_THREADS} ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_kanpig.bam &
            ${SAMTOOLS_COMMAND} index -@ ${N_MINIMAP_THREADS} ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE}/${SAMPLE}_hapestry.bam &
            wait
            
            printIgvScript ${OUTPUT_DIR}/${WINDOW_UNDERSCORE}/${SAMPLE} ${SAMPLE} ${WINDOW_DOTS}
        done
    done
done
