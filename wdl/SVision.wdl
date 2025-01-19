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
workflow SVision {
    input {
        String sample_id
        File input_bam
        File input_bai
        Int min_sv_length = 10
        String model_file = "model_liteunet_256_8_16_32_32_32.pth"
        String access_file = "hg38.access.10M.bed"
        File reference_fa
        Int n_cores = 16
        Int mem_gb = 32
    }

    call SVisionImpl {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            min_sv_length = min_sv_length,
            model_file = model_file,
            access_file = access_file,
            reference_fa = reference_fa,
            n_cores = n_cores,
            mem_gb = mem_gb
    }

    output {
        File output_vcf_gz = SVisionImpl.vcf_gz
        File output_tbi = SVisionImpl.tbi
    }
}


task SVisionImpl {
    input {
        String sample_id
        File input_bam
        File input_bai
        Int min_sv_length
        String model_file
        String access_file
        File reference_fa
        Int n_cores
        Int mem_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 2*ceil(size(input_bam, "GB")) + ceil(size(reference_fa, "GB")) + 50
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    String svision_suffix = "svision_pro_v2.4.s5"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        MKL_INTERFACE_LAYER="LP64"
        export MKL_INTERFACE_LAYER
        
        source activate svision-pro-env
        ${TIME_COMMAND} SVision-pro \
            --device cpu --process_num ${N_THREADS} \
            --preset hifi \
            --detect_mode germline \
            --min_sv_size ~{min_sv_length} \
            --model_path ~{docker_dir}/svision/src/pre_process/~{model_file} \
            --access_path ~{docker_dir}/svision/src/pre_process/~{access_file} \
            --genome_path ~{reference_fa} \
            --sample_name ~{sample_id} \
            --target_path ~{input_bam} \
            --out_path .
        conda deactivate
        bcftools sort --max-mem $(( ~{mem_gb} - 2 )) --output-type z ~{sample_id}.~{svision_suffix}.vcf > ~{sample_id}.~{svision_suffix}.vcf.gz
        tabix -f ~{sample_id}.~{svision_suffix}.vcf.gz
    >>>

    output {
        File vcf_gz = work_dir + "/" + sample_id + "." + svision_suffix + ".vcf.gz"
        File tbi = work_dir + "/" + sample_id + "." + svision_suffix + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
