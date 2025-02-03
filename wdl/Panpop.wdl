version 1.0


#
workflow Panpop {
    input {
        String sample_id
        File bcftools_merge_vcf_gz
        File bcftools_merge_tbi
        File reference_fa
        File reference_fai
        Int n_cpu
        Int ram_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    call PanpopImpl {
        input:
            sample_id = sample_id,
            bcftools_merge_vcf_gz = bcftools_merge_vcf_gz,
            bcftools_merge_tbi = bcftools_merge_tbi,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            n_cpu = n_cpu,
            ram_gb = ram_gb,
            disk_size_gb = disk_size_gb
    }
    output {
        File output_vcf_gz = PanpopImpl.output_vcf_gz
        File output_tbi = PanpopImpl.output_tbi
    }
}


#
task PanpopImpl {
    input {
        String sample_id
        File bcftools_merge_vcf_gz
        File bcftools_merge_tbi
        File reference_fa
        File reference_fai
        Int n_cpu
        Int ram_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        EFFECTIVE_MEM_GB=$(( ~{ram_gb} - 2 ))
        PANPOP_COMMAND="perl ~{docker_dir}/panpop-NC2024/bin/PART_run.pl"
        
        
        ${TIME_COMMAND} ${PANPOP_COMMAND} -t ${N_THREADS} --in_vcf ~{bcftools_merge_vcf_gz} -r ~{reference_fa} --tmpdir ./tmpdir1/ -o ./output1/
        ls -laht ./tmpdir1/
        ${TIME_COMMAND} ${PANPOP_COMMAND} -t ${N_THREADS} --in_vcf ./output1/3.final.vcf.gz -r ~{reference_fa} --tmpdir ./tmpdir2/ -o ./output2/ -not_first_merge
        ls -laht ./tmpdir2/
        tree -a
        mv ./tmpdir2/3.final.vcf.gz ~{sample_id}.panpop.vcf.gz
        tabix -f ~{sample_id}.panpop.vcf.gz
    >>>
    output {
        File output_vcf_gz = work_dir + "/" + sample_id + ".panpop.vcf.gz"
        File output_tbi = work_dir + "/" + sample_id + ".panpop.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
