version 1.0


# Performance on a machine with 16 cores and 32 GB of RAM, HG002,
# min_sv_length=10:
#
# COVERAGE  TIME    %CPU    RAM
# 4x        1m      900%    300m
# 8x        1m      1000%   600m
# 16x       50m     30%     400m       
# 32x       1h      30%     400m
#
workflow Sniffles {
    input {
        String sample_id
        File input_bam
        File input_bai
        Int min_sv_length = 10
        File reference_fa
        File tandems_bed
        Int n_cores = 16
        Int mem_gb = 32
    }

    call SnifflesImpl {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            min_sv_length = min_sv_length,
            reference_fa = reference_fa,
            tandems_bed = tandems_bed,
            n_cores = n_cores,
            mem_gb = mem_gb
    }

    output {
         File output_vcf_gz = SnifflesImpl.vcf_gz
         File output_tbi = SnifflesImpl.tbi
         File output_snf = SnifflesImpl.snf
    }
}


task SnifflesImpl {
    input {
        String sample_id
        File input_bam
        File input_bai
        Int min_sv_length
        File reference_fa
        File tandems_bed
        Int n_cores
        Int mem_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 2*ceil(size(input_bam, "GB")) + ceil(size(reference_fa, "GB")) + 50
    String docker_dir = "/hapestry_experiments"
    String work_dir = "/cromwell_root/hapestry_experiments"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"

        ${TIME_COMMAND} sniffles --threads ${N_THREADS} \
            --minsvlen ~{min_sv_length} \
            --input ~{input_bam} \
            --reference ~{reference_fa} \
            --tandem-repeats ~{tandems_bed} \
            --sample-id ~{sample_id} \
            --vcf ~{sample_id}.sniffles.vcf \
            --snf ~{sample_id}.sniffles.snf
        bgzip ~{sample_id}.sniffles.vcf
        tabix -f ~{sample_id}.sniffles.vcf.gz
    >>>

    output {
        File vcf_gz = work_dir + "/" + sample_id + ".sniffles.vcf.gz"
        File tbi = work_dir + "/" + sample_id + ".sniffles.vcf.gz.tbi"
        File snf = work_dir + "/" + sample_id + ".sniffles.snf"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
