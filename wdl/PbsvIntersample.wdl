version 1.0


#
workflow PbsvIntersample {
    input {
        Array[File] svsig
        Int min_sv_length
        File reference_fa
        Int n_cores = 32
        Int mem_gb = 128
    }

    call PbsvIntersampleImpl {
        input:
            svsig = svsig,
            min_sv_length = min_sv_length,
            reference_fa = reference_fa,
            n_cores = n_cores,
            mem_gb = mem_gb
    }

    output {
         File output_vcf_gz = PbsvIntersampleImpl.vcf_gz
         File output_tbi = PbsvIntersampleImpl.tbi
    }
}


task PbsvIntersampleImpl {
    input {
        Array[File] svsig
        Int min_sv_length
        File reference_fa
        Int n_cores
        Int mem_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 2*ceil(size(svsig, "GB")) + ceil(size(reference_fa, "GB")) + 100
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        REGIONS=""
        for i in $(seq 1 22) X Y M; do
            REGIONS="${REGIONS} chr${i}"
        done
        
        INPUT_FILES=~{sep=',' svsig}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        for INPUT_FILE in ${INPUT_FILES}; do
            tabix -f ${INPUT_FILE}
            mv ${INPUT_FILE}* .
        done
        for REGION in ${REGIONS}; do
            ${TIME_COMMAND} pbsv call \
                --num-threads 1 \
                --ccs \
                --min-sv-length ~{min_sv_length} \
                --region ${REGION} \
                ~{reference_fa} *.${REGION}.svsig.gz ${REGION}.pbsv_joint.vcf &
        done
        wait
        ls *.pbsv_joint.vcf > list.txt
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --file-list list.txt --output-type z > pbsv_joint.vcf.gz
        tabix -f pbsv_joint.vcf.gz
    >>>

    output {
        File vcf_gz = work_dir + "/pbsv_joint.vcf.gz"
        File tbi = work_dir + "/pbsv_joint.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
