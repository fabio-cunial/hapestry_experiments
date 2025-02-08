version 1.0


#
workflow Trgt {
    input {
        String sample_id
        String sex
        File input_bam
        File input_bai
        File ref_fa
        File ref_fai
        File repeat_catalog
        String repeat_catalog_name
        Int n_cpu = 64
        Int ram_gb = 64
    }
    parameter_meta {
        sex: "M or F"
    }
    
    call TrgtImpl {
        input:
            sample_id = sample_id,
            sex = sex,
            input_bam = input_bam,
            input_bai = input_bai,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            repeat_catalog = repeat_catalog,
            repeat_catalog_name = repeat_catalog_name,
            n_cpu = n_cpu,
            ram_gb = ram_gb
    }
    output {
        File output_vcf_gz = TrgtImpl.output_vcf_gz
        File output_tbi = TrgtImpl.output_tbi
    }
}


# Performance on a VM with 64 cores and 64GB of RAM:
#
# COVERAGE  CPU     RAM     TIME
# 32x       ???%    ??G    ??
#
task TrgtImpl {
    input {
        String sample_id
        String sex
        File input_bam
        File input_bai
        File ref_fa
        File ref_fai
        File repeat_catalog
        String repeat_catalog_name
        Int n_cpu
        Int ram_gb
    }
    parameter_meta {
        sex: "M or F"
    }

    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    Int disk_size_gb = 100 + 5*ceil(size(input_bam,"GB"))
    String output_name = sample_id + ".trgt." + repeat_catalog_name

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        EFFECTIVE_MEM_GB=$(( ~{ram_gb} - 2 ))
        
        if [ ~{sex} = "F" ]; then
            KARYOTYPE="XX"
        else
            KARYOTYPE="XY"
        fi
        ${TIME_COMMAND} trgt genotype \
            --sample-name ~{sample_id} \
            --genome ~{ref_fa} \
            --repeats ~{repeat_catalog} \
            --reads ~{input_bam} \
            --threads ${N_THREADS} \
            --output-prefix tmp1 \
            --karyotype ${KARYOTYPE}

        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp1.vcf.gz > ~{output_name}.vcf.gz
        tabix -f ~{output_name}.vcf.gz
    >>>
    output {
        File output_vcf_gz = work_dir + "/" + output_name + ".vcf.gz"
        File output_tbi = work_dir + "/" + output_name + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
