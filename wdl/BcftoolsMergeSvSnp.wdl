version 1.0


# Merges the raw joint SNP GVCF emitted by `LRJointCallGVCFs.wdl` with the
# joint normalized SV VCF created by `BcftoolsMergeIntersample.wdl`.
#
workflow BcftoolsMergeSvSnp {
    input {
        File sv_vcf_gz
        File sv_tbi
        File snp_vcf_gz
        File snp_tbi
        Int n_cpu = 8
        Int ram_gb = 200
    }
    parameter_meta {
    }
    
    call Merge {
        input:
            sv_vcf_gz = sv_vcf_gz,
            sv_tbi = sv_tbi,
            snp_vcf_gz = snp_vcf_gz,
            snp_tbi = snp_tbi,
            n_cpu = n_cpu,
            ram_gb = ram_gb
    }
    
    output {
        File output_vcf_gz = Merge.output_vcf_gz
        File output_tbi = Merge.output_tbi
    }
}


#
task Merge {
    input {
        File sv_vcf_gz
        File sv_tbi
        File snp_vcf_gz
        File snp_tbi
        Int n_cpu
        Int ram_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 2*ceil(size(sv_vcf_gz, "GB")) + 2*ceil(size(snp_vcf_gz, "GB")) + 100
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        REGIONS="chr1"
        for i in $(seq 2 22) X Y M; do
            REGIONS="${REGIONS},chr${i}"
        done
        
        # - Restricting the SNP VCF to standard chromosomes
        # - Removing multiallelic records from the SNP VCF
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --rm-dup exact --multiallelics - --do-not-normalize --output-type z --regions $(echo ${REGIONS}) ~{snp_vcf_gz} > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        rm -f ~{snp_vcf_gz}*
        ls -laht; df -h
        
        # - Merging
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --rm-dups exact --allow-overlaps --output-type z tmp1.vcf.gz > merged.vcf.gz
        tabix -f merged.vcf.gz
        ls -laht; df -h
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/merged.vcf.gz"
        File output_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
