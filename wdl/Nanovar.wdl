version 1.0


# Performance on a machine with 16 cores and 32 GB of RAM, HG002,
# min_sv_length=10:
#
# COVERAGE  TIME    %CPU    RAM
# 4x        
# 8x        
# 16x       
# 32x       
#
workflow Nanovar {
    input {
        String sample_id
        File input_bam
        File input_bai
        Int min_sv_length = 10
        File reference_fa
        Int n_cores = 16
        Int mem_gb = 32
    }

    call NanovarImpl {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            min_sv_length = min_sv_length,
            reference_fa = reference_fa,
            n_cores = n_cores,
            mem_gb = mem_gb
    }

    output {
        
    }
}


task NanovarImpl {
    input {
        String sample_id
        File input_bam
        File input_bai
        Int min_sv_length
        File reference_fa
        Int n_cores
        Int mem_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(input_bam, "GB")) + 10*ceil(size(reference_fa, "GB")) + 100
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    String nanovar_suffix = "nanovar"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        
        ${TIME_COMMAND} nanovar \
            --threads ${N_THREADS} \
            --data_type pacbio-ccs \
            --minlen ~{min_sv_length} \
            ~{input_bam} \
            ~{reference_fa} \
            ./working
        tree
        #bcftools sort --max-mem $(( ~{mem_gb} - 2 )) --output-type z ~{sample_id}.~{nanovar_suffix}.vcf > ~{sample_id}.~{nanovar_suffix}.vcf.gz
        #tabix -f ~{sample_id}.~{nanovar_suffix}.vcf.gz
    >>>

    output {
        
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
