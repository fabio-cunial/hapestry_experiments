version 1.0


#
workflow BcftoolsMergeIntersample {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
        Int ram_gb = 200
    }
    parameter_meta {
    }
    
    call InterSampleMerge {
        input:
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi,
            ram_gb = ram_gb
    }
    
    output {
        File output_vcf_gz = InterSampleMerge.output_vcf_gz
        File output_tbi = InterSampleMerge.output_tbi
    }
}



# Performance on 1074 AoU 8x samples:
# COMMAND           RUNTIME     N_CPUS      MAX_RSS
# bcftools merge    2.5h        2           140G
# bcftools norm     40m         3           80G
#
task InterSampleMerge {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
        Int ram_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(input_vcf_gz, "GB")) + 100
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_gb} - 2 ))
        
        INPUT_FILES=~{sep=',' input_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            echo ${INPUT_FILE} >> list.txt
        done
        ls -laht; tree; df -h
        # $--info-rules -$ disables default rules, and it is used just to avoid
        # the following error:
        # Only fixed-length vectors are supported with -i sum:AC
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --info-rules - --file-list list.txt --output-type z > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        ls -laht; tree; df -h
        
        # Removing multiallelic records, if any are generated during the merge.
        # This is just an extra safeguard and might be dropped if $bcftools
        # merge --merge none$ always behaves correctly.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp1.vcf.gz > merged.vcf.gz
        tabix -f merged.vcf.gz
        rm -f tmp1.vcf.gz*
        ls -laht; tree; df -h
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
