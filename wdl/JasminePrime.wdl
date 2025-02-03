version 1.0

# In intra-sample jasmine we use:
#
# --allow_intrasample min_seq_id=0.9
#
# In particular, `--allow_intrasample` is necessary, since the bcftools merge
# VCF in input contains just one sample column and variants from different
# callers. `min_seq_id=0.9` mimics what is done in a production pipeline like
# AoU, but in jasmine it only works with INS.
#
# When we do inter-sample jasmine on the VCFs emitted by the intra-sample
# jasmine above, we just use:
#
# min_seq_id=0.7
#
# In particular, `--allow_intrasample` is not enabled, since we assume that
# intra-sample merging has already been performed. `min_seq_id=0.7` mimics what
# is done in a production pipeline like AoU. Note that `min_seq_id` defaults to
# zero, which means that the sequence similarity of INS is not taken into
# account.
#
# When we do inter-sample jasmine on the bcftools merge of all raw calls from
# all callers, we use:
#
# --allow_intrasample min_seq_id=0.7
#
# Once again, `--allow_intrasample` is needed, since variants from the same
# sample might come from different callers, and `min_seq_id=0.7` mimics the AoU
# setting above, even though this is a different input VCF.
#
workflow JasminePrime {
    input {
        String sample_id
        Array[File] input_vcf_gz
        Array[File] input_tbi
        String jasmine_params = " "
        Int n_cpu
        Int ram_gb
    }
    parameter_meta {
    }
    
    call JasminePrimeImpl {
        input:
            sample_id = sample_id,
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi,
            jasmine_params = jasmine_params,
            n_cpu = n_cpu,
            ram_gb = ram_gb
    }
    output {
        File output_vcf_gz = JasminePrimeImpl.output_vcf_gz
        File output_tbi = JasminePrimeImpl.output_tbi
    }
}


#
task JasminePrimeImpl {
    input {
        String sample_id
        Array[File] input_vcf_gz
        Array[File] input_tbi
        String jasmine_params
        Int n_cpu
        Int ram_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(input_vcf_gz,"GB")) + 100
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        EFFECTIVE_MEM_GB=$(( ~{ram_gb} - 2 ))
        
        INPUT_FILES=~{sep=',' input_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        i="0"
        for INPUT_FILE in ${INPUT_FILES}; do
            i=$(( ${i} + 1 ))
            gunzip -c ${INPUT_FILE} > ${i}.vcf
            echo ${i}.vcf >> list.txt
        done
        ${TIME_COMMAND} java -cp /opt/conda/bin/jasmine.jar -Xms${EFFECTIVE_MEM_GB}G -Xmx${EFFECTIVE_MEM_GB}G Main threads=${N_THREADS} --output_genotypes ~{jasmine_params} file_list=list.txt out_file=~{sample_id}.jasmine.vcf
        bcftools view --no-header ~{sample_id}.jasmine.vcf | head -n 10
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z ~{sample_id}.jasmine.vcf > ~{sample_id}.jasmine.vcf.gz
        tabix -f ~{sample_id}.jasmine.vcf.gz
    >>>
    output {
        File output_vcf_gz = work_dir + "/" + sample_id + ".jasmine.vcf.gz"
        File output_tbi = work_dir + "/" + sample_id + ".jasmine.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
