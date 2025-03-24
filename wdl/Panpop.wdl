version 1.0


#
workflow Panpop {
    input {
        String sample_id
        File bcftools_merge_vcf_gz
        File bcftools_merge_tbi
        File reference_fa
        File reference_fai
        Int max_sv_length = 50000
        Int n_cpu
        Int ram_gb
        Int disk_size_gb
    }
    parameter_meta {
        max_sv_length: "SVs that are too long make even an intra-sample merge too slow. 0=do not remove long SVs."
    }
    
    call PanpopImpl {
        input:
            sample_id = sample_id,
            bcftools_merge_vcf_gz = bcftools_merge_vcf_gz,
            bcftools_merge_tbi = bcftools_merge_tbi,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            max_sv_length = max_sv_length,
            n_cpu = n_cpu,
            ram_gb = ram_gb,
            disk_size_gb = disk_size_gb
    }
    output {
        File output_vcf_gz = PanpopImpl.output_vcf_gz
        File output_tbi = PanpopImpl.output_tbi
    }
}


# Performance on a VM with 16 cores and 16GB of RAM:
#
# COVERAGE  CPU     RAM     TIME
# 4x        600%    3G      2m
# 8x        1200%   3G      3m
# 16x       1200%   9G      3m
# 32x       1300%   9G      5m
#
task PanpopImpl {
    input {
        String sample_id
        File bcftools_merge_vcf_gz
        File bcftools_merge_tbi
        File reference_fa
        File reference_fai
        Int max_sv_length
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
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        EFFECTIVE_MEM_GB=$(( ~{ram_gb} - 2 ))
        PANPOP_COMMAND="perl ~{docker_dir}/panpop-NC2024/bin/PART_run.pl"
        
        # - Removing calls that are too long
        if [ ~{max_sv_length} -eq 0 ]; then
            mv ~{bcftools_merge_vcf_gz} tmp1.vcf.gz
            mv ~{bcftools_merge_tbi} tmp1.vcf.gz.tbi
        else
            bcftools filter --include "INFO/SVLEN>=-~{max_sv_length} && INFO/SVLEN<=~{max_sv_length}" --output-type z ~{bcftools_merge_vcf_gz} > tmp1.vcf.gz
            tabix -f tmp1.vcf.gz
        fi
        
        # - Running panpop
        source activate panpop
        cpanm MCE::Channel Tie::CharArray
        ${TIME_COMMAND} ${PANPOP_COMMAND} -t ${N_THREADS} --in_vcf tmp1.vcf.gz -r ~{reference_fa} --tmpdir ./tmpdir1/ -o ./output1/
        ls -laht ./tmpdir1/
        rm -f tmp1.vcf.gz*
        ${TIME_COMMAND} ${PANPOP_COMMAND} -t ${N_THREADS} --in_vcf ./output1/3.final.vcf.gz -r ~{reference_fa} --tmpdir ./tmpdir2/ -o ./output2/ -not_first_merge
        ls -laht ./tmpdir2/
        tree -a
        conda deactivate
        
        # - Removing multiallelic records, which might have been created by
        #   panpop.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z ./output2/3.final.vcf.gz > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        
        # - Outputting
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp2.vcf.gz > ~{sample_id}.panpop.vcf.gz
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
