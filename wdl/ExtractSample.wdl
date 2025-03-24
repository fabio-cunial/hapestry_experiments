version 1.0


# Extracts a specific sample from a multi-sample VCF and keeps only records that
# occur in that sample.
#
workflow ExtractSample {
    input {
        String sample_id
        File vcf_gz
        File tbi
    }
    parameter_meta {
    }
    
    call ExtractImpl {
        input:
            sample_id = sample_id,
            vcf_gz = vcf_gz,
            tbi = tbi
    }
    
    output {
    	File extracted_vcf_gz = ExtractImpl.extracted_vcf_gz
    	File extracted_tbi = ExtractImpl.extracted_tbi
    }
}


#
task ExtractImpl {
    input {
        String sample_id
        File vcf_gz
        File tbi
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    Int disk_size_gb = 100 + ceil(size(vcf_gz,"GB"))
    Int ram_gb = 8
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ~{sample_id} --output-type z ~{vcf_gz} > tmp.vcf.gz
        tabix -f tmp.vcf.gz
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'GT="alt"' --output-type z tmp.vcf.gz > ~{sample_id}.vcf.gz
        tabix -f ~{sample_id}.vcf.gz
    >>>
    
    output {
    	File extracted_vcf_gz = work_dir + "/" + sample_id + ".vcf.gz"
    	File extracted_tbi = work_dir + "/" + sample_id + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 2
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
