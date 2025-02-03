version 1.0


# Performance on a machine with 16 cores and 32 GB of RAM, HG002,
# min_sv_length=10:
#
# COVERAGE  TIME    %CPU    RAM
# 4x        31s     1300%   800m
# 8x        1m40s   700%    850m
# 16x       9m14s   200%    900m
# 32x       11m     200%    980m
#
workflow Cutesv {
    input {
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        String cutesv_params_ccs = "--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --merge_ins_threshold 500 --merge_del_threshold 500"
        String cutesv_params_ont = "--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --merge_ins_threshold 500 --merge_del_threshold 500"
        Int is_ont
        Int n_cpu
        Int ram_size_gb
    }
    parameter_meta {
        cutesv_params_ccs: "Default values are taken from the command-line help."
        cutesv_params_ont: "Default values are taken from the command-line help."
        is_ont: "0/1"
    }
    
    call SingleSampleCalling {
        input:
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            cutesv_params_ccs = cutesv_params_ccs,
            cutesv_params_ont = cutesv_params_ont,
            is_ont = is_ont,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb
    }
    
    output {
        File vcf_gz = SingleSampleCalling.vcf_gz
        File tbi = SingleSampleCalling.tbi
    }
}


task SingleSampleCalling {
    input {
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        String cutesv_params_ccs
        String cutesv_params_ont
        Int is_ont
        Int n_cpu
        Int ram_size_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    String output_prefix = "cutesv"
    Int disk_size_gb = 200 + ceil(size(reference_fa,"GB")) + 10*ceil(size(alignments_bam,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        df -h
        
        if [ ~{is_ont} -eq 1 ]; then
            PARAMS=$(echo ~{cutesv_params_ont})
        else
            PARAMS=$(echo ~{cutesv_params_ccs})
        fi
        mkdir ./cutesv_tmp
        ${TIME_COMMAND} cuteSV --threads ${N_THREADS} \
            ${PARAMS} \
            ~{alignments_bam} \
            ~{reference_fa} \
            ~{output_prefix}.vcf \
            ./cutesv_tmp
        rm -rf ./cutesv_tmp
        bgzip ~{output_prefix}.vcf
        tabix -f ~{output_prefix}.vcf.gz
    >>>

    output {
        File vcf_gz = work_dir + "/" + output_prefix + ".vcf.gz"
        File tbi = work_dir + "/" + output_prefix + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
