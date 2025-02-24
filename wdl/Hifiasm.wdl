version 1.0

#
workflow Hifiasm {
    input {
        String sample_id
        File reads_fq
        Int n_cpu = 64
        Int ram_gb = 128
    }

    call Assemble {
        input:
            sample_id  = sample_id,
            reads_fq = reads_fq,
            n_cpu = n_cpu,
            ram_gb = ram_gb
    }

    output {
        File hap1_fa_gz = Assemble.hap1_fa_gz
        File hap2_fa_gz = Assemble.hap2_fa_gz
    }
}


# Performance on a VM with 64 cores and 128GB of RAM:
#
# COVERAGE  CPU     RAM     TIME
# 4x
# 8x        2100%   51.5G   1h
# 16x
# 32x       2300%   54G     3h
#
task Assemble {
    input {
        String sample_id
        File reads_fq
        Int n_cpu
        Int ram_gb
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    Int disk_size_gb = 10 * ceil(size(reads_fq, "GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"

        ${TIME_COMMAND} ~{docker_dir}/hifiasm/hifiasm -t ${N_THREADS} -o ~{sample_id} ~{reads_fq}
        ls -laht
        tree

        # Outputting
        for haplotype_gfa in ~{sample_id}.bp.hap*.p_ctg.gfa ; do
            filename=$(basename ${haplotype_gfa})
            haplotype="${filename%.*}"
            awk '/^S/{print ">"$2; print $3}' ${haplotype_gfa} | gzip -1 -c > ${haplotype}.fa.gz
        done
    >>>

    output {
        File hap1_fa_gz = work_dir + "/" + sample_id + ".bp.hap1.p_ctg.fa.gz"
        File hap2_fa_gz = work_dir + "/" + sample_id + ".bp.hap2.p_ctg.fa.gz"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
