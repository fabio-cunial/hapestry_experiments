version 1.0


# Normalizes multiallelic sites and then filters SVs by length.
#
workflow FilterDipcall {
    input {
        File input_vcf_gz
        String min_sv_length
    }
    parameter_meta {
        min_sv_length: "Comma-separated. Example: 10,20,30,40,50"
    }
    
    call FilterDipcallImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            min_sv_length = min_sv_length
    }
    
    output {
        Array[File] vcf_gz = FilterDipcallImpl.vcf_gz
        Array[File] tbi = FilterDipcallImpl.tbi
    }
}


task FilterDipcallImpl {
    input {
        File input_vcf_gz
        String min_sv_length
    }
    parameter_meta {
        min_sv_length: "Comma-separated. Example: 10,20,30,40,50"
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
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        bcftools norm --multiallelics - --output-type v ~{input_vcf_gz} > input.vcf
        rm -f ~{input_vcf_gz}
        echo ~{min_sv_length} | tr ',' '\n' > lengths.txt
        while read LENGTH; do
            ${TIME_COMMAND} java -cp ~{docker_dir} FilterDipcall input.vcf ${LENGTH} dipcall_${LENGTH}.vcf
            ${TIME_COMMAND} bgzip dipcall_${LENGTH}.vcf
            tabix -f dipcall_${LENGTH}.vcf.gz
        done < lengths.txt
    >>>
    
    output {
        Array[File] vcf_gz = glob(work_dir+"dipcall_*.vcf.gz")
        Array[File] tbi = glob(work_dir+"dipcall_*.vcf.gz.tbi")
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 1
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
