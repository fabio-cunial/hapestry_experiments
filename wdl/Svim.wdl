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
workflow Svim {
    input {
        String sample_id
        File input_bam
        File input_bai
        Int min_sv_length = 10
        File reference_fa
        Int n_cores = 16
        Int mem_gb = 32
    }

    call SvimImpl {
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
        File vcf_gz = SvimImpl.vcf_gz
        File tbi = SvimImpl.tbi
    }
}


task SvimImpl {
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
    String svim_suffix = "svim"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        
        ${TIME_COMMAND} svim alignment \
            --min_sv_size ~{min_sv_length} \
            --tandem_duplications_as_insertions \
            --interspersed_duplications_as_insertions \
            --sample ~{sample_id} \
            ./working \
            ~{input_bam} \
            ~{reference_fa}
        tree
        bcftools sort --max-mem $(( ~{mem_gb} - 2 )) --output-type z ./working/variants.vcf > ~{sample_id}.~{svim_suffix}.vcf.gz
        tabix -f ~{sample_id}.~{svim_suffix}.vcf.gz
    >>>

    output {
        File vcf_gz = work_dir + "/" + sample_id + "." + svim_suffix + ".vcf.gz"
        File tbi = work_dir + "/" + sample_id + "." + svim_suffix + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
