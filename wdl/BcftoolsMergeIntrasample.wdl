version 1.0


#
workflow BcftoolsMergeIntrasample {
    input {
        String sample_id
        Array[File] sample_vcf_gz
        Array[File] sample_tbi
        Int ram_gb = 8
        Int compression_level = 1
    }
    parameter_meta {
    }
    
    call IntrasampleMerge {
        input:
            sample_id = sample_id,
            sample_vcf_gz = sample_vcf_gz,
            sample_tbi = sample_tbi,
            ram_gb = ram_gb,
            compression_level = compression_level
    }
    
    output {
        File output_vcf_gz = IntrasampleMerge.output_vcf_gz
        File output_tbi = IntrasampleMerge.output_tbi
    }
}


# Remark: we use $bcftools merge$ instead of $bcftools concat$, since we must
# collapse identical calls made by different callers (otherwise they would
# remain in the inter-sample VCF, since the inter-sample $bcftools merge$ does
# not collapse records from the same sample).
#
task IntrasampleMerge {
    input {
        String sample_id
        Array[File] sample_vcf_gz
        Array[File] sample_tbi
        Int ram_gb
        Int compression_level
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(sample_vcf_gz, "GB")) + 10
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    Int n_vcfs = length(sample_vcf_gz)
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_gb} - 2 ))
        
        # Merging all single-caller VCFs
        INPUT_FILES=~{sep=',' sample_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            echo ${INPUT_FILE} >> list.txt
        done
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list list.txt --output-type z > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        
        # Removing multiallelic records, if any are generated during the merge.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp1.vcf.gz > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        rm -f tmp1.vcf.gz*
        
        # Picking the first nonzero GT of every call
        bcftools view --header-only tmp2.vcf.gz > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > ~{sample_id}_bcftools_merge.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{sample_id}" >> ~{sample_id}_bcftools_merge.vcf
        ${TIME_COMMAND} bcftools view --no-header tmp2.vcf.gz | awk '{ \
            printf("%s",$1); \
            for (i=2; i<=8; i++) printf("\t%s",$i); \
            gt="0/1"; \
            for (i=9+1; i<=9+~{n_vcfs}; i++) { \
                value=substr($i,1,3); \
                if (value=="0/1" || value=="1/0" || value=="1/1" || value=="0|1" || value=="1|0" || value=="1|1" || value=="1") { \
                    gt=value; \
                    break; \
                } \
            } \
            printf("\tGT\t%s\n",gt); \
        }' >> ~{sample_id}_bcftools_merge.vcf
        
        # Outputting
        ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level ~{compression_level} ~{sample_id}_bcftools_merge.vcf
        tabix -f ~{sample_id}_bcftools_merge.vcf.gz
        ls -laht; tree
    >>>

    output {
        File output_vcf_gz = work_dir + "/" + sample_id + "_bcftools_merge.vcf.gz"
        File output_tbi = work_dir + "/" + sample_id + "_bcftools_merge.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 2
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
