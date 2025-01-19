version 1.0


# Performance on a machine with 16 cores and 32 GB of RAM, HG002,
# min_sv_length=10:
#
# COVERAGE      TIME    %CPU    RAM
# 4x discover   6m      100%    1G 
# 4x call       1h30m   200%    42G
# 8x discover   12m     100%    2.3G
# 8x call       2h30m   250%    44G
# 16x discover  30m     100%    4.5G
# 16x call      ---------> out of ram
# 32x discover  38m     100%    5.5G
# 32x call      ---------> out of ram
#
workflow Pbsv {
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

    call PbsvImpl {
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
         File output_vcf_gz = PbsvImpl.vcf_gz
         File output_tbi = PbsvImpl.tbi
         File output_svsig = PbsvImpl.svsig
    }
}


task PbsvImpl {
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
        MIN_SVSIG_LENGTH=$(( ~{min_sv_length} - 3 ))  # Arbitrary

        #${TIME_COMMAND} pbsv discover \
        #    --ccs \
        #    --sample ~{sample_id} \
        #    --min-svsig-length ${MIN_SVSIG_LENGTH} \
        #    --tandem-repeats ~{tandems_bed} \
        #    ~{input_bam} ~{sample_id}.svsig.gz
        for REGION in $(samtools view -H ~{input_bam} | grep '^@SQ' | cut -f2 | cut -d ':' -f2); do
            ${TIME_COMMAND} pbsv discover \
                --region ${REGION} \
                --ccs \
                --sample ~{sample_id} \
                --min-svsig-length ${MIN_SVSIG_LENGTH} \
                --tandem-repeats ~{tandems_bed} \
                ~{input_bam} ~{sample_id}.${REGION}.svsig.gz &
        done
        wait
        ${TIME_COMMAND} pbsv call \
            --num-threads ${N_THREADS} \
            --ccs \
            --min-sv-length ~{min_sv_length} \
            ~{reference_fa} *.svsig.gz ~{sample_id}.pbsv.vcf
        bgzip ~{sample_id}.pbsv.vcf
        tabix ~{sample_id}.pbsv.vcf.gz
    >>>

    output {
        File vcf_gz = work_dir + "/" + sample_id + ".pbsv.vcf.gz"
        File tbi = work_dir + "/" + sample_id + ".pbsv.vcf.gz.tbi"
        File svsig = work_dir + "/" + sample_id + ".svsig.gz"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
