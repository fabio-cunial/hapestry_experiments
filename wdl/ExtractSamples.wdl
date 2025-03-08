version 1.0


#
workflow ExtractSamples {
    input {
        File cohort_vcf_gz
        File cohort_tbi
        Int include_missing
    }
    parameter_meta {
        include_missing: "(0/1) Includes ./. calls in the output."
    }
    
    call ExtractImpl {
        input:
            cohort_vcf_gz = cohort_vcf_gz,
            cohort_tbi = cohort_tbi,
            include_missing = include_missing
    }
    
    output {
    	Array[File] extracted_vcf_gz = ExtractImpl.extracted_vcf_gz
    	Array[File] extracted_tbi = ExtractImpl.extracted_tbi
    }
}


#
task ExtractImpl {
    input {
        File cohort_vcf_gz
        File cohort_tbi
        Int include_missing
    }
    parameter_meta {
        include_missing: "(0/1) Includes ./. calls in the output."
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    Int disk_size_gb = 100 + ceil(size(cohort_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        FILTER_STRING_WITHOUT_MISSING='COUNT(GT="alt")>0'
        FILTER_STRING_WITH_MISSING='COUNT(GT="RR")=0'
        if [ ~{include_missing} -eq 1 ]; then
            FILTER_STRING=${FILTER_STRING_WITH_MISSING}
        else
            FILTER_STRING=${FILTER_STRING_WITHOUT_MISSING}
        fi
        
        
        function extractThread() {
            SAMPLES_FILE=$1
        
            while read SAMPLE; do
                ${TIME_COMMAND} bcftools view --samples ${SAMPLE} ~{cohort_vcf_gz} | bcftools filter --include "${FILTER_STRING}" --output-type z > sample_${SAMPLE}.vcf.gz
                tabix -f sample_${SAMPLE}.vcf.gz
            done < ${SAMPLES_FILE}
        }
        
        
        bcftools view --header-only ~{cohort_vcf_gz} > header.txt
        tail -n 1 header.txt | cut -f 10- | tr '\t' '\n' > samples.txt
        N_SAMPLES=$(wc -l < samples.txt)
        N_ROWS=$(( (${N_SAMPLES}+${N_THREADS}-1)/${N_THREADS} ))  # Ceiling
        split -d -l ${N_ROWS} samples.txt samples_
        for SAMPLES_FILE in $(find . -maxdepth 1 -name 'samples_*' | sort); do
            extractThread ${SAMPLES_FILE} &
        done
        wait
    >>>
    
    output {
    	Array[File] extracted_vcf_gz = glob(work_dir+"/sample_*.vcf.gz")
    	Array[File] extracted_tbi = glob(work_dir+"/sample_*.vcf.gz.tbi")
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 16
        memory: "16GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
