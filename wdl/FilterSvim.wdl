version 1.0


#
workflow FilterSvim {
    input {
        String sample_id
        File input_vcf_gz
        File input_tbi
        String filter_string = ""
    }
    parameter_meta {
        filter_string: "With single quotes. Examples: 'QUAL>=10', 'INFO/SUPPORT>=2'"
    }
    
    call FilterSvimImpl {
        input:
            sample_id = sample_id,
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi,
            filter_string = filter_string
    }
    
    output {
        File vcf_gz = FilterSvimImpl.vcf_gz
        File tbi = FilterSvimImpl.tbi
    }
}


task FilterSvimImpl {
    input {
        String sample_id
        File input_vcf_gz
        File input_tbi
        String filter_string
    }
    
    Int disk_size_gb = 10*ceil(size(input_vcf_gz,"GB"))
    Int ram_size_gb = 4
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

        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include ~{filter_string} --output-type z ~{input_vcf_gz} > ~{sample_id}_filtered.vcf.gz
        tabix -f filtered.vcf.gz
    >>>
    
    output {
        File vcf_gz = work_dir + "/" + sample_id + "_filtered.vcf.gz"
        File tbi = work_dir + "/" + sample_id + "_filtered.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 2
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
