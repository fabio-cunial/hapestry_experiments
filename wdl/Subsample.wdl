version 1.0


# Downloads BAM files up to a given coverage, randomizing both the order in
# which the BAMs are downloaded, and the order of the reads in the
# concatenation of all downloaded BAMs.
#
# Remark: the output FASTQ file for coverage $i+1$ is not necessarily a
# superset of the output FASTQ file for coverage $i$.
#
# Performance on HPRC's HG00621. Requested coverages: 4,8,16,32.
# 16 physical cores, 128 GB RAM, 500 GB HDD.
#
# STEP                  TIME            CPU         RAM
# samtools fastq        10m             80%         13M
# seqkit stats          1m              100%        30M
# seqkit scat           2h              130%        4G
# total                 6h30m
#
workflow Subsample {
    input {
        String sample_id
        Array[String] bam_addresses
        String coverages = "4,8,16,32"
        String remote_dir
        String billing_project = "broad-firecloud-dsde-methods"
        Int haploid_genome_length_gb = 3
        Int n_cores = 8
        Int mem_gb = 32
        Int disk_size_gb = 500
    }
    parameter_meta {
        bam_addresses: "Can be .bam, .fastq, .fastq.gz"
        coverages: "Comma-separated"
        remote_dir: "Output directory in a remote bucket"
        n_cores: "At least one per coverage"
    }
    
    call SubsampleImpl {
        input:
            sample_id = sample_id,
            bam_addresses = bam_addresses,
            coverages = coverages,
            remote_dir = remote_dir,
            billing_project = billing_project,
            haploid_genome_length_gb = haploid_genome_length_gb,
            n_cores = n_cores,
            mem_gb = mem_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
    }
}


task SubsampleImpl {
    input {
        String sample_id
        Array[String] bam_addresses
        String coverages
        String remote_dir
        String billing_project
        Int haploid_genome_length_gb
        Int n_cores
        Int mem_gb
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
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TOTAL_RAM_KB=$(cat /proc/meminfo | grep MemTotal | awk '{print $2}')
        TOTAL_RAM_GB=$(( ${TOTAL_RAM_KB} / 1000000 ))
        RAM_PER_THREAD_BYTES=$(( (1000000000*( ${TOTAL_RAM_GB} - 10)) / ${N_THREADS} ))
        HAPLOID_GENOME_LENGTH_GB=$(( ~{haploid_genome_length_gb} * 1000000000 ))
        df -h
        
        MAX_COVERAGE=~{coverages}
        MAX_COVERAGE=${MAX_COVERAGE##*,}
        
        # 1. Randomizing the order of the BAMs and downloading enough of them
        LIST_FILE=~{write_lines(bam_addresses)}
        shuf ${LIST_FILE} > randomized.txt
        cat randomized.txt
        TARGET_N_CHARS=$(( ${MAX_COVERAGE} * ${HAPLOID_GENOME_LENGTH_GB} ))
        TOTAL_N_CHARS="0"; TOTAL_N_READS="0"; FILE_ID="0";
        while read ADDRESS; do
            FILE_ID=$(( ${FILE_ID} + 1 ))
            SUCCESS="0"
            if [[ ${ADDRESS} == gs://* ]]; then
                SUCCESS=$(gsutil -u ~{billing_project} cp ${ADDRESS} . && echo 1 || echo 0)
            elif [[ ${ADDRESS} == s3://* ]]; then
                aws s3 --no-sign-request cp ${ADDRESS} .
                SUCCESS="1"
            else
                SUCCESS=$(wget ${ADDRESS} && echo 1 || echo 0)
            fi
            if [[ ${SUCCESS} -eq 0 ]]; then
                continue
            fi
            FILE_NAME=$(basename ${ADDRESS})
            if [[ ${FILE_NAME} == *.bam ]]; then
                ${TIME_COMMAND} samtools fastq -@ ${N_THREADS} -n ${FILE_NAME} | pigz --processes ${N_THREADS} --fast --to-stdout > ${FILE_ID}.fastq.gz
            elif [[ ${FILE_NAME} == *.fastq.gz ]]; then
                mv ${FILE_NAME} ${FILE_ID}.fastq.gz
            elif [[ ${FILE_NAME} == *.fastq ]]; then
                ${TIME_COMMAND} pigz --processes ${N_THREADS} --fast --to-stdout ${FILE_NAME} > ${FILE_ID}.fastq.gz
            fi
            rm -f ${FILE_NAME}
            ${TIME_COMMAND} ~{docker_dir}/seqkit stats --threads ${N_THREADS} --tabular ${FILE_ID}.fastq.gz > stats.txt
            N_CHARS=$(cut -f 5 stats.txt | tail -n 1)
            TOTAL_N_CHARS=$(( ${TOTAL_N_CHARS} + ${N_CHARS} ))
            N_READS=$(cut -f 4 stats.txt | tail -n 1)
            TOTAL_N_READS=$(( ${TOTAL_N_READS} + ${N_READS} ))
            if [[ ${TOTAL_N_CHARS} -gt ${TARGET_N_CHARS} ]]; then
                break
            fi
            df -h
        done < randomized.txt
        mkdir ./fastqs
        mv *.fastq.gz ./fastqs
        if [ ${FILE_ID} -eq 1 ]; then
            mv ./fastqs/1.fastq.gz tmp1.fastq.gz
        else
            #${TIME_COMMAND} ~{docker_dir}/seqkit scat --threads ${N_THREADS} -f --out-format fastq ./fastqs > tmp1.fastq.gz
            ${TIME_COMMAND} zcat ./fastqs/*.fastq.gz | gzip -c --fast > tmp1.fastq.gz
        fi
        rm -rf ./fastqs/
        df -h
        
        # 2. Subsampling
        function sample() {
            local INPUT_FASTQ_GZ=$1
            local INPUT_FASTQ_GZ_N_CHARS=$2
            local COVERAGE=$3
            
            TARGET_N_CHARS=$(( ${COVERAGE} * ${HAPLOID_GENOME_LENGTH_GB} ))
            if [ ${TARGET_N_CHARS} -ge ${INPUT_FASTQ_GZ_N_CHARS} ]; then
                echo "WARNING: the remote files contain just ${INPUT_FASTQ_GZ_N_CHARS} bps in total, but we need ${TARGET_N_CHARS} bps to achieve coverage ${COVERAGE}."
                cp ${INPUT_FASTQ_GZ} ~{sample_id}_${COVERAGE}.fastq.gz
            else
                PROPORTION=$( bc <<< "scale=2; ${TARGET_N_CHARS}/${INPUT_FASTQ_GZ_N_CHARS}" )
                ${TIME_COMMAND} ~{docker_dir}/seqkit sample --proportion ${PROPORTION} -o ~{sample_id}_${COVERAGE}.fastq.gz ${INPUT_FASTQ_GZ}
            fi
        }
        echo ~{coverages} | tr ',' '\n' > coverages.txt
        while read COVERAGE; do
            sample tmp1.fastq.gz ${TOTAL_N_CHARS} ${COVERAGE} &
        done < coverages.txt
        wait
        rm -f tmp1.fastq.gz
        
        # 3. Uploading
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{sample_id}_'*.fastq.gz' ~{remote_dir} && echo 0 || echo 1)
            if [[ ${TEST} -eq 1 ]]; then
                echo "Error uploading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
