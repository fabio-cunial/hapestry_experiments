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
        Int n_cpu = 32
        Int ram_gb = 16
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
# resolved version where records with missing or ref GT are discarded, and
# multiallelic records are split into multiple records. 
#
# Remark: both the resolved and the raw VCFs are sorted.
#
# Performance on a VM with 64 cores and 64GB of RAM:
#
# COVERAGE  CPU     RAM     TIME
# 32x       3100%   11G     1h46m
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
        Int n_cpu
        Int ram_gb
    }
    parameter_meta {
        sex: "M or F"
    }

    String docker_dir = "/hapestry"
    Int disk_size_gb = 50 + 2*ceil(size(input_bam,"GB")) + 2*ceil(size(repeat_catalog,"GB"))

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
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
        
        # - Discarding records with missing or ref GT
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --exclude 'GT="mis" || GT="0/0"' --output-type z tmp2.vcf.gz > tmp3.vcf.gz
        tabix -f tmp3.vcf.gz
        
        # - Removing multiallelic records, which are created by TRGT.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp3.vcf.gz > tmp4.vcf.gz
        tabix -f tmp4.vcf.gz
        rm -f tmp3.vcf*
        
        # Outputting
        mv tmp2.vcf.gz ~{sample_id}.trgt.vcf.gz
        mv tmp2.vcf.gz.tbi ~{sample_id}.trgt.vcf.gz.tbi
        mv tmp4.vcf.gz ~{sample_id}.trgt.resolved.vcf.gz
        mv tmp4.vcf.gz.tbi ~{sample_id}.trgt.resolved.vcf.gz.tbi
    >>>
    output {
        File output_vcf_gz = work_dir + "/" + sample_id + ".trgt.vcf.gz"
        File output_tbi = work_dir + "/" + sample_id + ".trgt.vcf.gz.tbi"
        File output_resolved_vcf_gz = work_dir + "/" + sample_id + ".trgt.resolved.vcf.gz"
        File output_resolved_tbi = work_dir + "/" + sample_id + ".trgt.resolved.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
