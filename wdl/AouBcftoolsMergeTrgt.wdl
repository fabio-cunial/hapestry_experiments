version 1.0


#
workflow AouBcftoolsMergeTrgt {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
    }
    parameter_meta {
    }
    
    scatter (i in range(length(input_vcf_gz))) {
        call Resolve {
            input:
                input_vcf_gz = input_vcf_gz[i],
                input_tbi = input_tbi[i]
        }
    }
    call Merge {
        input:
            input_vcf_gz = Resolve.output_vcf_gz,
            input_tbi = Resolve.output_tbi
    }
    
    output {
        File output_vcf_gz = Merge.output_vcf_gz
        File output_tbi = Merge.output_tbi
    }
}


task Resolve {
    input {
        File input_vcf_gz
        File input_tbi
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
        
        # - Sorting
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z ~{input_vcf_gz} > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        
        # - Removing multiallelic records, which might be created by TRGT.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp1.vcf.gz > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        rm -f tmp1.vcf.gz*
        
        # - Discarding records with missing or ref GT
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --exclude 'GT="mis" || GT="0/0"' --output-type z tmp2.vcf.gz > resolved.vcf.gz
        tabix -f resolved.vcf.gz
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/resolved.vcf.gz"
        File output_tbi = work_dir + "/resolved.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 2
        memory: "8GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


# Performance on 1074 AoU 8x samples:
#
# COMMAND           RUNTIME     N_CPUS      MAX_RSS
# bcftools merge    
#
task Merge {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
        
        Int ram_gb = 256
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil(size(input_vcf_gz, "GB")) + 100
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
            echo ${INPUT_FILE}.vcf.gz >> list.txt
        done
        
        # $--info-rules -$ disables default rules, and it is used just to avoid
        # the following error:
        # Only fixed-length vectors are supported with -i sum:AC
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --info-rules - --file-list list.txt --output-type z > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        ls -laht
        
        # Removing multiallelic records, if any are generated during the merge.
        # This is just an extra safeguard and might be dropped if $bcftools
        # merge --merge none$ always behaves correctly.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp1.vcf.gz > merged.vcf.gz
        tabix -f merged.vcf.gz
        ls -laht
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
