version 1.0

import "Resolve.wdl" as resolve


#
workflow PbsvIntersample {
    input {
        Array[File] svsig
        Array[String] regions
        Int min_sv_length
        File reference_fa
        Int n_cores_call = 32
        Int mem_gb_call = 128
        Int n_cores_concat = 32
        Int mem_gb_concat = 128
    }

    scatter(region in regions) {
        call CallRegion {
            input:
                region = region,
                svsig = svsig,
                min_sv_length = min_sv_length,
                reference_fa = reference_fa,
                n_cores = n_cores_call,
                mem_gb = mem_gb_call
        }
        call resolve.Resolve as res {
            input:
                sample_id = region,
                vcf_gz = CallRegion.vcf_gz,
                tbi = CallRegion.tbi,
                reference_fa = reference_fa
        }
    }
    call Concat {
        input:
            vcf_gz = res.resolved_vcf_gz,
            tbi = res.resolved_tbi,
            n_cores = n_cores_concat,
            mem_gb = mem_gb_concat
    }

    output {
         File output_vcf_gz = Concat.vcf_gz
         File output_tbi = Concat.tbi
    }
}


#
task CallRegion {
    input {
        String region
        Array[File] svsig
        Int min_sv_length
        File reference_fa
        Int n_cores
        Int mem_gb
    }
    parameter_meta {
        svsig: "All SVSIG files, from all regions (wasteful but easier to implement)."
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
        
        INPUT_FILES=~{sep=',' svsig}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        for INPUT_FILE in ${INPUT_FILES}; do
            mv ${INPUT_FILE} .
        done
        ${TIME_COMMAND} pbsv call \
            --num-threads ${N_THREADS} \
            --ccs \
            --min-sv-length ~{min_sv_length} \
            ~{reference_fa} *.~{region}.svsig.gz ~{region}.pbsv_joint.vcf
        bgzip ~{region}.pbsv_joint.vcf
        tabix -f ~{region}.pbsv_joint.vcf.gz
    >>>

    output {
        File vcf_gz = work_dir + "/" + region + ".pbsv_joint.vcf.gz"
        File tbi = work_dir + "/" + region + ".pbsv_joint.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


#
task Concat {
    input {
        Array[File] vcf_gz
        Array[File] tbi
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
        
        INPUT_FILES=~{sep=',' vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            echo ${INPUT_FILE} >> list.txt
        done
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --file-list list.txt --output-type z > pbsv_joint.vcf.gz
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
