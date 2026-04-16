version 1.0


#
workflow VcfdistEvaluationCompareHapestryKanpig {
    input {
        String sample_id
        String remote_input_dir
        String coverage_id
        File sample_dipcall_vcf_gz
        File sample_dipcall_tbi
        
        Int min_sv_length
        String? region = "chr1"
        File supercluster_components_bed
        String vcfdist_extra_args = ""
        
        File reference_fa
        File reference_fai

        String docker_image = "fcunial/hapestry_experiments"
    }
    parameter_meta {
        sample_dipcall_vcf_gz: "Assumed to contain SNPs and short INDELs, but no multiallelics."
    }
    
    call Impl { 
        input:
            sample_id = sample_id,
            remote_input_dir = remote_input_dir,
            coverage_id = coverage_id,
            sample_dipcall_vcf_gz = sample_dipcall_vcf_gz,
            sample_dipcall_tbi = sample_dipcall_tbi,
            
            min_sv_length = min_sv_length,
            region = region,
            supercluster_components_bed = supercluster_components_bed,
            vcfdist_extra_args = vcfdist_extra_args,
            
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            
            docker_image = docker_image
    }
}


# Performance on an 8-core, 32GB VM, chr1 only, 50bp 32x:
#
# COMMAND           CPU         RAM         TIME        COST
# vcfdist           
#
task Impl {
    input {
        String sample_id
        String remote_input_dir
        String coverage_id
        File sample_dipcall_vcf_gz
        File sample_dipcall_tbi
        
        Int min_sv_length
        String? region
        File supercluster_components_bed
        String vcfdist_extra_args

        File reference_fa
        File reference_fai

        String docker_image
        Int n_cpu = 8
        Int ram_size_gb = 32
        Int disk_size_gb = 50
        Int preemptible_number = 0
    }

    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        cat << 'END' > vcfdist_on_chunk.sh
#!/bin/bash
set -euxo pipefail

HAPESTRY_VCF_GZ=$1
KANPIG_VCF_GZ=$2
TRUTH_VCF_GZ=$3
EFFECTIVE_RAM_GB=$4
ID=$5

rm -f ${ID}.out
while read -u 3 ROW; do
    echo ${ROW} | tr ',' '\t' > ${ID}.bed
    CHROM=$(echo ${ROW} | cut -f 1)
    START=$(echo ${ROW} | cut -f 2)
    END=$(echo ${ROW} | cut -f 3)
    
    # Truth
    bcftools view --output-type z ${TRUTH_VCF_GZ} ${CHROM}:${START}-${END} --output ${ID}_truth.vcf.gz
    bcftools index -f -t ${ID}_truth.vcf.gz
    
    # Hapestry
    bcftools view --output-type z ${HAPESTRY_VCF_GZ} ${CHROM}:${START}-${END} --output ${ID}_query.vcf.gz
    bcftools index -f -t ${ID}_query.vcf.gz
    vcfdist ${ID}_query.vcf.gz ${ID}_truth.vcf.gz reference.fa \
        --sv-threshold ~{min_sv_length} \
        \
        --largest-variant 10000 \
        --max-supercluster-size 10002 \
        \
        --max-threads 1 --max-ram ${EFFECTIVE_RAM_GB} --verbosity 1 \
        --realign-query --realign-truth \
        --bed ${ID}.bed \
        --distance \
        ~{vcfdist_extra_args} \
        --prefix ${ID}_
    HAPESTRY_STRING=$(grep 'ALL' ${ID}_distance-summary.tsv | grep 'BEST')
    rm -f ${ID}_* 
    
    # Kanpig
    bcftools view --output-type z ${KANPIG_VCF_GZ} ${CHROM}:${START}-${END} --output ${ID}_query.vcf.gz
    bcftools index -f -t ${ID}_query.vcf.gz
    vcfdist ${ID}_query.vcf.gz ${ID}_truth.vcf.gz reference.fa \
        --sv-threshold ~{min_sv_length} \
        \
        --largest-variant 10000 \
        --max-supercluster-size 10002 \
        \
        --max-threads 1 --max-ram ${EFFECTIVE_RAM_GB} --verbosity 1 \
        --realign-query --realign-truth \
        --bed ${ID}.bed \
        --distance \
        ~{vcfdist_extra_args} \
        --prefix ${ID}_
    KANPIG_STRING=$(grep 'ALL' ${ID}_distance-summary.tsv | grep 'BEST')
    rm -f ${ID}_* 
    
    echo -e "${HAPESTRY_STRING}\t${KANPIG_STRING}" >> ${ID}.out
    rm -f ${ID}.bed
done 3< chunk_${ID}
END
        chmod +x vcfdist_on_chunk.sh
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        # Preparing the input VCFs
        ${TIME_COMMAND} gcloud storage cp ~{remote_input_dir}/~{min_sv_length}bp/~{coverage_id}/hapestry/~{sample_id}_extracted.vcf.gz ./hapestry.vcf.gz
        ${TIME_COMMAND} gcloud storage cp ~{remote_input_dir}/~{min_sv_length}bp/~{coverage_id}/hapestry/~{sample_id}_extracted.vcf.gz.tbi ./hapestry.vcf.gz.tbi
        ${TIME_COMMAND} gcloud storage cp ~{remote_input_dir}/~{min_sv_length}bp/~{coverage_id}/kanpig/~{sample_id}_extracted.vcf.gz ./kanpig.vcf.gz
        ${TIME_COMMAND} gcloud storage cp ~{remote_input_dir}/~{min_sv_length}bp/~{coverage_id}/kanpig/~{sample_id}_extracted.vcf.gz.tbi ./kanpig.vcf.gz.tbi
        mv ~{sample_dipcall_vcf_gz} dipcall.vcf.gz
        mv ~{sample_dipcall_tbi} dipcall.vcf.gz.tbi
        mv ~{reference_fa} reference.fa
        mv ~{reference_fai} reference.fa.fai
        if ~{defined(region)}
        then
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --output-type z hapestry.vcf.gz ~{region} --output query.vcf.gz
            rm -f hapestry.vcf.gz* ; mv query.vcf.gz hapestry.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t hapestry.vcf.gz
            
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --output-type z kanpig.vcf.gz ~{region} --output query.vcf.gz
            rm -f kanpig.vcf.gz* ; mv query.vcf.gz kanpig.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t kanpig.vcf.gz
            
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --output-type z dipcall.vcf.gz ~{region} --output truth.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t truth.vcf.gz
        else
            mv dipcall.vcf.gz truth.vcf.gz
            mv dipcall.vcf.gz.tbi truth.vcf.gz.tbi
        fi
        
        # Evaluating every window in parallel
        cat ~{supercluster_components_bed} | tr '\t' ',' > superclusters.csv
        rm -f ~{supercluster_components_bed}
        N_COMPONENTS=$(wc -l < superclusters.csv)
        N_COMPONENTS_PER_THREAD=$(( ${N_COMPONENTS} / ${N_THREADS} ))
        split -d -a 2 -l ${N_COMPONENTS_PER_THREAD} superclusters.csv chunk_
        N_FILES=$(ls chunk_* | wc -l)
        ls chunk_* | sort -V | cut -c 7- > list.txt
        ${TIME_COMMAND} xargs --arg-file=list.txt --max-lines=1 --max-procs=${N_THREADS} ./vcfdist_on_chunk.sh hapestry.vcf.gz kanpig.vcf.gz truth.vcf.gz ${EFFECTIVE_RAM_GB}
        ls -laht 1>&2
        df -h 1>&2
        
        # Concatenating results
        rm -f out.tsv
        while read -u 3 ID; do
            cat ${ID}.out >> out.tsv
        done 3< list.txt
    >>>

    output {
        File out_tsv = "out.tsv"
    }

    runtime {
        docker: docker_image
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: ram_size_gb + " GiB"
        cpu: n_cpu
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
