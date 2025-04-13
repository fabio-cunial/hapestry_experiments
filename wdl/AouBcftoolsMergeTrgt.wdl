version 1.0


# Resolves and merges the raw, single-sample TRGT calls in AoU. Merges the
# resulting cohort VCF with an SV-only cohort VCF, keeping only SV calls that
# overlap with the complement of the TRGT intervals.
#
workflow AouBcftoolsMergeTrgt {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
        
        File sv_merge_vcf_gz
        File sv_merge_tbi
        
        File samples_file
        File reference_fai
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
    call Concat {
        input:
            sv_merge_vcf_gz = sv_merge_vcf_gz,
            sv_merge_tbi = sv_merge_tbi,
            trgt_merge_vcf_gz = Merge.output_vcf_gz,
            trgt_merge_tbi = Merge.output_tbi,
            samples_file = samples_file,
            reference_fai = reference_fai
    }
    
    output {
        File output_vcf_gz = Concat.output_vcf_gz
        File output_tbi = Concat.output_tbi
    }
}


task Resolve {
    input {
        File input_vcf_gz
        File input_tbi
        
        Int ram_gb = 4
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
        
        # Sorting
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ~{input_vcf_gz} > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        
        # Removing multiallelic records, which might be created by TRGT.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp1.vcf.gz > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        rm -f tmp1.vcf.gz*
        
        # Discarding records with missing or ref GT
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
        memory: ram_gb + "GB"
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
        
        INPUT_FILES=~{sep=',' input_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            echo ${INPUT_FILE} >> list.txt
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


#
task Concat {
    input {
        File sv_merge_vcf_gz
        File sv_merge_tbi
        File trgt_merge_vcf_gz
        File trgt_merge_tbi
        
        File samples_file
        File reference_fai
        
        Int n_cpu = 4
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil( size(sv_merge_vcf_gz, "GB") + size(trgt_merge_vcf_gz, "GB") ) + 100
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
        
        # Removing SVs inside TRGT intervals
        ${TIME_COMMAND} bedtools complement -i ~{trgt_merge_vcf_gz} -g ~{reference_fai} > not_trgt.bed
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --regions-file not_trgt.bed --regions-overlap variant --output-type z ~{sv_merge_vcf_gz} > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        rm -f ~{sv_merge_vcf_gz}
        
        # Ensuring that samples have the same order in both files
        ${TIME_COMMAND} bcftools view --samples-file ~{samples_file} --output-type z tmp1.vcf.gz > tmp2.vcf.gz
        tabix -f tmp1.vcf.gz
        ${TIME_COMMAND} bcftools view --samples-file ~{samples_file} --output-type z ~{trgt_merge_vcf_gz} > tmp3.vcf.gz
        tabix -f tmp3.vcf.gz
        rm -f ~{trgt_merge_vcf_gz}
        
        # Combining SV and TRGT calls
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --output-type z tmp2.vcf.gz tmp3.vcf.gz > out.vcf.gz
        tabix -f out.vcf.gz
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/out.vcf.gz"
        File output_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}