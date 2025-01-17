version 1.0


#
workflow MapCCS {
    input {
        String sample_id
        File reference_fa
        File reference_fai
        File reads_fastq_gz
        Int n_cpus
        Int ram_size_gb
        Int disk_gb
    }
    parameter_meta {
    }
    
    call MapCCSImpl {
        input:
            sample_id = sample_id,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            reads_fastq_gz = reads_fastq_gz,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb,
            disk_gb = disk_gb
    }
    
    output {
        File output_bam = MapCCSImpl.output_bam
        File output_bai = MapCCSImpl.output_bai
    }
}


task MapCCSImpl {
    input {
        String sample_id
        File reference_fa
        File reference_fai
        File reads_fastq_gz
        Int n_cpus
        Int ram_size_gb
        Int disk_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = ceil(size(reads_fastq_gz, "GB")) + ceil(size(reference_fa, "GB")) + disk_gb
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Remark: pbmm2 automatically uses all cores.
        ${TIME_COMMAND} pbmm2 align --preset CCS --sort --sample ~{sample_id} ~{reference_fa} ~{reads_fastq_gz} out.bam
        ${TIME_COMMAND} samtools calmd -@ ${N_THREADS} --no-PG -b out.bam ~{reference_fa} > ~{sample_id}.bam
        ${TIME_COMMAND} samtools index -@ ${N_THREADS} ~{sample_id}.bam
    >>>
    
    output {
        File output_bam = work_dir + "/" + sample_id + ".bam"
        File output_bai = work_dir + "/" + sample_id + ".bam.bai"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
