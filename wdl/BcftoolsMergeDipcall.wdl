version 1.0


#
workflow BcftoolsMergeDipcall {
    input {
        Array[String] sample_id
        Array[File] sample_vcf_gz
        Array[File] sample_tbi
        Int ram_gb = 32
    }
    parameter_meta {
    }
    
    call BcftoolsMergeDipcallImpl {
        input:
            sample_id = sample_id,
            sample_vcf_gz = sample_vcf_gz,
            sample_tbi = sample_tbi,
            ram_gb = ram_gb
    }
    
    output {
        File output_vcf_gz = BcftoolsMergeDipcallImpl.output_vcf_gz
        File output_tbi = BcftoolsMergeDipcallImpl.output_tbi
    }
}


task BcftoolsMergeDipcallImpl {
    input {
        Array[String] sample_id
        Array[File] sample_vcf_gz
        Array[File] sample_tbi
        Int ram_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    Int disk_size_gb = 10*ceil(size(sample_vcf_gz, "GB")) + 100
    Int n_files = length(sample_vcf_gz)
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        INPUT_FILES=~{sep=',' sample_vcf_gz}
        SAMPLE_IDS=~{sep=',' sample_id}
        rm -f list.txt
        for i in $(seq 1 ~{n_files}); do
            # - Enforcing the right sample name
            INPUT_FILE=$(echo ${INPUT_FILES} | cut -d , -f ${i})
            SAMPLE_ID=$(echo ${SAMPLE_IDS} | cut -d , -f ${i})
            echo ${SAMPLE_ID} > samples.txt
            bcftools reheader --samples samples.txt ${INPUT_FILE} > tmp1.vcf.gz
            rm -f samples.txt
            tabix -f tmp1.vcf.gz
            # - Removing multiallelic records
            bcftools norm --multiallelics - --output-type z tmp1.vcf.gz > ${SAMPLE_ID}.vcf.gz
            rm -f tmp1.vcf.gz*
            tabix -f ${SAMPLE_ID}.vcf.gz
            echo ${SAMPLE_ID}.vcf.gz >> list.txt
        done
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list list.txt --output-type z > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        
        # - Removing multiallelic records one last time
        bcftools norm --multiallelics - --output-type z tmp1.vcf.gz > merged.vcf.gz
        rm -f tmp1.vcf.gz*
        tabix -f merged.vcf.gz
        ls -laht; tree
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/merged.vcf.gz"
        File output_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 4
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
