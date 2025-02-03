version 1.0


# Performance on a machine with 16 cores and 32 GB of RAM, HG002,
# min_sv_length=10:
#
# COVERAGE      TIME    %CPU    RAM
# 4x call       
#
# 8x call       
#
# 16x discover  
# 16x call      
#
# 32x discover  
# 32x call      
#
workflow Sawfish {
    input {
        String sample_id
        Boolean is_male
        File input_bam
        File input_bai
        Int min_sv_length = 10
        File reference_fa
        File expected_cn_male
        File expected_cn_female
        Int n_cores = 16
        Int mem_gb = 32
    }

    call SawfishImpl {
        input:
            sample_id = sample_id,
            is_male = is_male,
            input_bam = input_bam,
            input_bai = input_bai,
            min_sv_length = min_sv_length,
            reference_fa = reference_fa,
            expected_cn_male = expected_cn_male,
            expected_cn_female = expected_cn_female,
            n_cores = n_cores,
            mem_gb = mem_gb
    }

    output {
         File output_vcf_gz = SawfishImpl.vcf_gz
         File output_tbi = SawfishImpl.tbi
    }
}


task SawfishImpl {
    input {
        String sample_id
        Boolean is_male
        File input_bam
        File input_bai
        Int min_sv_length
        File reference_fa
        File expected_cn_male
        File expected_cn_female
        Int n_cores
        Int mem_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 2*ceil(size(input_bam, "GB")) + ceil(size(reference_fa, "GB")) + 50
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        MIN_INDEL_SIZE=$(( ~{min_sv_length} - 5 ))  # Arbitrary

        if [ ~{is_male} == true ]; then
            CN_FLAG="--expected-cn ~{expected_cn_male}"
        else
            CN_FLAG="--expected-cn ~{expected_cn_female}"
        fi
        export RUST_BACKTRACE="full"
        ${TIME_COMMAND} ~{docker_dir}/sawfish/bin/sawfish discover \
            --threads ${N_THREADS} \
            --ref ~{reference_fa} \
            --bam ~{input_bam} \
            --min-indel-size ${MIN_INDEL_SIZE} \
            --output-dir ./sawfish_discover_output
            ${CN_FLAG}
        ${TIME_COMMAND} ~{docker_dir}/sawfish/bin/sawfish joint-call \
            --threads ${N_THREADS} \
            --sample ./sawfish_discover_output \
            --output-dir ./sawfish_joint_call_dir
        mv ./sawfish_joint_call_dir/genotyped.sv.vcf.gz ./~{sample_id}.sawfish.vcf.gz
        tabix -f ./~{sample_id}.sawfish.vcf.gz
    >>>

    output {
        File vcf_gz = work_dir + "/" + sample_id + ".sawfish.vcf.gz"
        File tbi = work_dir + "/" + sample_id + ".sawfish.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
