version 1.0

#
workflow AlignAssembly {
    input {
        String sample_id

        File hap1_fa
        File hap2_fa
        
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            sample_id = sample_id,
            hap1_fa = hap1_fa,
            hap2_fa = hap2_fa,
            reference_fa = reference_fa,
            reference_fai = reference_fai
    }
    
    output {
    	File hap1_bam = Impl.hap1_bam
        File hap1_bai = Impl.hap1_bai

        File hap2_bam = Impl.hap2_bam
        File hap2_bai = Impl.hap2_bai
    }
}


#
task Impl {
    input {
        String sample_id

        File hap1_fa
        File hap2_fa
        File reference_fa
        File reference_fai

        Int ram_gb = 64
        Int n_cpu = 16
        Int disk_size_gb = 50
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))


        ${TIME_COMMAND} minimap2 -t ${N_THREADS} -a -x asm10 ~{reference_fa} ~{hap1_fa} > hap1.sam
        ${TIME_COMMAND} samtools sort -@ ${N_THREADS} -o ~{sample_id}_hap1.bam hap1.sam
        ${TIME_COMMAND} samtools index -@ ${N_THREADS} ~{sample_id}_hap1.bam  
        rm -f hap1.sam

        ${TIME_COMMAND} minimap2 -t ${N_THREADS} -a -x asm10 ~{reference_fa} ~{hap2_fa} > hap2.sam
        ${TIME_COMMAND} samtools sort -@ ${N_THREADS} -o ~{sample_id}_hap2.bam hap2.sam
        ${TIME_COMMAND} samtools index -@ ${N_THREADS} ~{sample_id}_hap2.bam
        rm -f hap2.sam
        
        ls -laht
    >>>
    
    output {
    	File hap1_bam = sample_id + "_hap1.bam"
        File hap1_bai = sample_id + "_hap1.bam.bai"

        File hap2_bam = sample_id + "_hap2.bam"
        File hap2_bai = sample_id + "_hap2.bam.bai"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}