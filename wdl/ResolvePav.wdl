version 1.0

import "Resolve.wdl" as resolve


#
workflow ResolvePav {
    input {
        String sample_id
        File input_vcf_gz
        Int min_sv_length
        File reference_fa
    }
    parameter_meta {
    }
    
    call ResolvePavImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            min_sv_length = min_sv_length
    }
    call resolve.ResolveImpl as r {
        input:
            sample_id = sample_id,
            vcf_gz = ResolvePavImpl.sv_vcf_gz,
            tbi = ResolvePavImpl.sv_tbi,
            reference_fa = reference_fa
    }
    output {
        File resolved_vcf_gz = r.resolved_vcf_gz
        File resolved_tbi = r.resolved_tbi
    }
}


task ResolvePavImpl {
    input {
        File input_vcf_gz
        Int min_sv_length
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(input_vcf_gz,"GB"))
    Int ram_size_gb = 4
    Int compression_level = 1
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
        
        # - Normalizing multiallelic records
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type v ~{input_vcf_gz} > input.vcf
        
        # - Keeping only long calls
        ${TIME_COMMAND} java -cp ~{docker_dir} PAV2SVs input.vcf ~{min_sv_length} sv.vcf snp.vcf
        ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level ~{compression_level} sv.vcf
        tabix -f sv.vcf.gz
    >>>
    
    output {
        File sv_vcf_gz = work_dir + "/sv.vcf.gz"
        File sv_tbi = work_dir + "/sv.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 4
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
