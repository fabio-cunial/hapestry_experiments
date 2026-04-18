version 1.0


workflow Vcf2Haplotypes {
    input {
        String sample_id

        File vcf_gz
        File tbi
        String region
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            sample_id = sample_id,
            vcf_gz = vcf_gz,
            tbi = tbi,
            region = region,
            reference_fa = reference_fa,
            reference_fai = reference_fai
    }
    
    output {
    	File hap1 = Impl.hap1
        File hap2 = Impl.hap2
    }
}


#
task Impl {
    input {
        String sample_id

        File vcf_gz
        File tbi
        String region
        File reference_fa
        File reference_fai

        Int ram_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(reference_fa,"GB") + size(vcf_gz,"GB"))
    String docker_dir = "/hapestry"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        bcftools view --output-type z ~{vcf_gz} ~{region} --output in.vcf.gz
        tabix -f in.vcf.gz
        ${TIME_COMMAND} bcftools consensus --fasta-ref  ~{reference_fa} --haplotype 1 in.vcf.gz --output ~{sample_id}_hap1.fa &
        ${TIME_COMMAND} bcftools consensus --fasta-ref  ~{reference_fa} --haplotype 2 in.vcf.gz --output ~{sample_id}_hap2.fa &
        wait
        ls -laht
    >>>
    
    output {
    	File hap1 = sample_id + "_hap1.fa"
        File hap2 = sample_id + "_hap2.fa"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 2
        memory: ram_gb + "GB"
        disks: "local-disk 50 HDD"
        preemptible: 0
    }
}