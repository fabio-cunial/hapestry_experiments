version 1.0


# Keeps only long calls.
#
workflow FilterLength {
    input {
        String sample_id
        File vcf_gz
        File tbi
        Int min_sv_length
    }
    parameter_meta {
        vcf_gz: "Assumed to be the output of `Resolve.wdl`."
    }
    
    call FilterImpl {
        input:
            sample_id = sample_id,
            vcf_gz = vcf_gz,
            tbi = tbi,
            min_sv_length = min_sv_length
    }
    
    output {
    	File filtered_vcf_gz = FilterImpl.filtered_vcf_gz
    	File filtered_tbi = FilterImpl.filtered_tbi
    }
}


#
task FilterImpl {
    input {
        String sample_id
        File vcf_gz
        File tbi
        Int min_sv_length
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    Int ram_gb = 8
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "SVLEN>=~{min_sv_length} || SVLEN<=-~{min_sv_length}" --output-type z ~{vcf_gz} > ~{sample_id}_filtered.vcf.gz
        tabix -f ~{sample_id}_filtered.vcf.gz
    >>>
    
    output {
    	File filtered_vcf_gz = work_dir + "/" + sample_id + "_filtered.vcf.gz"
    	File filtered_tbi = work_dir + "/" + sample_id + "_filtered.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 2
        memory: ram_gb + "GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
