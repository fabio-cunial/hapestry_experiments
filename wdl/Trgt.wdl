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
        File output_resolved_vcf_gz = TrgtImpl.output_resolved_vcf_gz
        File output_resolved_tbi = TrgtImpl.output_resolved_tbi
    }
}


# Remark: in addition to the raw output VCF from TRGT, the task outputs also a
# resolved version where ALT=. records are discarded and multiallelic records
# are split into multiple records. 
#
# Remark: both the resolved and the raw VCFs are sorted.
#
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
        
        # - Running TRGT
        if [ ~{sex} = "F" ]; then
            KARYOTYPE="XX"
        else
            KARYOTYPE="XY"
        fi
        ${TIME_COMMAND} ~{docker_dir}/trgt genotype \
            --threads ${N_THREADS} \
            --disable-bam-output \
            --sample-name ~{sample_id} \
            --karyotype ${KARYOTYPE} \
            --genome ~{ref_fa} \
            --repeats ~{repeat_catalog} \
            --reads ~{input_bam} \
            --output-prefix tmp1
        
        # - Sorting
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp1.vcf.gz > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        rm -f tmp1.vcf*
        
        # - Discarding ALT=. records, which mean GT=0/0.
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --exclude 'ALT="."' --output-type z tmp2.vcf.gz > tmp3.vcf.gz
        tabix -f tmp3.vcf.gz
        
        # - Removing multiallelic records, which are created by TRGT.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp3.vcf.gz > tmp4.vcf.gz
        tabix -f tmp4.vcf.gz
        rm -f tmp3.vcf*
        
        # Outputting
        mv tmp2.vcf.gz ~{output_name}.vcf.gz
        mv tmp2.vcf.gz.tbi ~{output_name}.vcf.gz.tbi
        mv tmp4.vcf.gz ~{output_name}.resolved.vcf.gz
        mv tmp4.vcf.gz.tbi ~{output_name}.resolved.vcf.gz.tbi
    >>>
    output {
        File output_vcf_gz = work_dir + "/" + output_name + ".vcf.gz"
        File output_tbi = work_dir + "/" + output_name + ".vcf.gz.tbi"
        File output_resolved_vcf_gz = work_dir + "/" + output_name + ".resolved.vcf.gz"
        File output_resolved_tbi = work_dir + "/" + output_name + ".resolved.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
