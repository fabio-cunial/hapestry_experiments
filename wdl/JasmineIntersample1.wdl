version 1.0

# When we do inter-sample jasmine on the VCFs emitted by intra-sample jasmine,
# we just use:
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
workflow JasmineIntersample1 {
    input {
        String sample_id
        Array[File] input_vcf_gz
        Array[File] input_tbi
        String jasmine_params = " "
        Int n_cpu
        Int ram_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    call JasmineImpl {
        input:
            sample_id = sample_id,
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi,
            jasmine_params = jasmine_params,
            n_cpu = n_cpu,
            ram_gb = ram_gb,
            disk_size_gb = disk_size_gb
    }
    output {
        File output_vcf_gz = JasmineImpl.output_vcf_gz
        File output_tbi = JasmineImpl.output_tbi
    }
}


#
task JasmineImpl {
    input {
        String sample_id
        Array[File] input_vcf_gz
        Array[File] input_tbi
        String jasmine_params
        Int n_cpu
        Int ram_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
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
        JAVA_PATH="/usr/bin/java"  # Using the latest JRE
        


        INPUT_FILES=~{sep=',' input_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        i="0"
        for INPUT_FILE in ${INPUT_FILES}; do
            i=$(( ${i} + 1 ))
            gunzip -c ${INPUT_FILE} > ${i}.vcf
            echo ${i}.vcf >> list.txt
        done
        ${TIME_COMMAND} ${JAVA_PATH} -jar /opt/conda/bin/jasmine.jar -Xms${EFFECTIVE_MEM_GB}G -Xmx${EFFECTIVE_MEM_GB}G threads=${N_THREADS} --output_genotypes ~{jasmine_params} file_list=list.txt out_file=tmp2.vcf
        
        # - Removing a suffix of the INFO field added by Jasmine, since it makes
        #   bcftools sort crash.
        bcftools view --header-only --no-version tmp2.vcf > tmp3.vcf
        N_ROWS=$(wc -l < tmp3.vcf)
        tail -n +$(( ${N_ROWS} + 1 )) tmp2.vcf > body.txt
        rm -f tmp2.vcf
        ${TIME_COMMAND} cat body.txt | awk '{ \
            pattern="ALLVARS_EXT"; \
            \
            printf("%s",$1); \
            for (i=2; i<8; i++) printf("\t%s",$i); \
            i=match($8,pattern); \
            if (i!=0) printf("\t%s",substr($8,1,i-1)); \
            else printf("\t%s",$8); \
            for (i=9; i<=NF; i++) printf("\t%s",$i); \
            printf("\n"); \
        }' >> tmp3.vcf
        rm -f body.txt
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp3.vcf > ~{sample_id}.jasmine.vcf.gz
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
