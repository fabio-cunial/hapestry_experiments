version 1.0


#
workflow TrgtMerge {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
        File reference_fa
        File reference_fai
        Int n_cpu
        Int ram_gb
    }
    
    call TrgtMergeImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            n_cpu = n_cpu,
            ram_gb = ram_gb
    }
    output {
        File output_vcf_gz = TrgtMergeImpl.output_vcf_gz
        File output_tbi = TrgtMergeImpl.output_tbi
    }
}


#
task TrgtMergeImpl {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
        File reference_fa
        File reference_fai
        Int n_cpu
        Int ram_gb
    }

    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    Int disk_size_gb = 50 + 4*ceil(size(input_vcf_gz,"GB")) + ceil(size(reference_fa,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        
        INPUT_FILES=~{sep=',' input_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        ${TIME_COMMAND} ~{docker_dir}/trgt merge --vcf ${INPUT_FILES} --genome ~{reference_fa} --output-type z --output merged.trgt.vcf.gz
        ls -laht
        tree
    >>>
    output {
        File output_vcf_gz = work_dir + "/merged.trgt.vcf.gz"
        File output_tbi = work_dir + "/merged.trgt.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
