version 1.0


#
workflow TrgtMerge {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
        File reference_fa
        File reference_fai
        Int n_cpu = 1
        Int ram_gb = 8
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


# Performance on a VM with 4 cores and 128GB of RAM:
#
# COVERAGE              CPU     RAM     TIME
# 32x, 3 samples        100%    60M     3m
# 32x, 12 samples       100%    200M    10m
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
        
        # - Running TRGT
        INPUT_FILES=~{sep=',' input_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        ${TIME_COMMAND} ~{docker_dir}/trgt merge --vcf ${INPUT_FILES} --genome ~{reference_fa} --output-type z --output tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        
        # - Discarding records with AF=0
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --exclude 'COUNT(GT!="mis" && GT!="0/0")=0' --output-type z tmp1.vcf.gz > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        rm -f tmp1.vcf*
        
        # - Removing multiallelic records, which are both in the input files and
        #   are created by trgt merge.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp2.vcf.gz > tmp3.vcf.gz
        tabix -f tmp3.vcf.gz
        rm -f tmp2.vcf*
        
        # Outputting
        mv tmp3.vcf.gz merged.trgt.vcf.gz
        mv tmp3.vcf.gz.tbi merged.trgt.vcf.gz.tbi
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
