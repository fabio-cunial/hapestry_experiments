version 1.0


#
workflow SnifflesIntersample {
    input {
        Array[File] input_snf
        Int ram_gb
    }

    call JointCalling {
        input:
            input_snf = input_snf,
            ram_gb = ram_gb
    }

    output {
         File output_vcf_gz = JointCalling.output_vcf_gz
         File output_tbi = JointCalling.output_tbi
    }
}


task JointCalling {
    input {
        Array[File] input_snf
        Int ram_gb
    }

    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    Int disk_size_gb = 10*ceil(size(input_snf, "GB")) + 100

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}

        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        
        INPUT_FILES=~{sep=',' input_snf}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        ${TIME_COMMAND} sniffles --threads ${N_THREADS} --input ${INPUT_FILES} --vcf sniffles_joint.vcf
        bgzip sniffles_joint.vcf
        tabix sniffles_joint.vcf.gz
        ls -laht
        tree
    >>>

    output {
        File output_vcf_gz = work_dir + "/sniffles_joint.vcf.gz"
        File output_tbi = work_dir + "/sniffles_joint.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 16
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}