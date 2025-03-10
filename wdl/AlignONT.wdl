version 1.0


# 
#
workflow AlignONT {
    input {
        String sample_id
        File reference_fa
        File reference_fai
        File reads_fastq_gz
        Boolean is_r10
        Int n_cpus
        Int ram_size_gb
        Int disk_gb
    }
    parameter_meta {
    }
    
    call AlignONTImpl {
        input:
            sample_id = sample_id,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            reads_fastq_gz = reads_fastq_gz,
            is_r10 = is_r10,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb,
            disk_gb = disk_gb
    }
    
    output {
        File output_bam = AlignONTImpl.output_bam
        File output_bai = AlignONTImpl.output_bai
    }
}


task AlignONTImpl {
    input {
        String sample_id
        File reference_fa
        File reference_fai
        File reads_fastq_gz
        Boolean is_r10
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
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        FLAGS=""
        if [[ ~{is_r10} == true ]]; then
            FLAGS="-x lr:hq"
        else
            FLAGS="-x map-ont"
        fi
        ${TIME_COMMAND} ~{docker_dir}/minimap2/minimap2 -t ${N_THREADS} ${FLAGS} -ayYL --MD --eqx --cs ~{reference_fa} ~{reads_fastq_gz} > out.sam
        ${TIME_COMMAND} samtools sort -@ ${N_THREADS} --no-PG -O BAM out.sam > ~{sample_id}.bam
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
