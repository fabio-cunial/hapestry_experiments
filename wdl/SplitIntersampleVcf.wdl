version 1.0


#
task SplitIntersampleVcf {
    input {
        File vcf_gz
        File vcf_gz_tbi
        Int min_sv_length
        String coverage_id
        String caller_id
        
        String remote_outdir
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 2
    }

    Int disk_size_gb = 1 + ceil(3 * (size(vcf_gz, "GiB")))
    String docker_dir = "/hapestry"

    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        
        # Splitting
        bcftools view --header-only ~{vcf_gz} | tail -n 1 | tr '\t' '\n' | tail -n +10 > samples.txt
        ${TIME_COMMAND} bcftools +split --samples-file samples.txt --output-type z --output . ~{vcf_gz}
        
        # Keeping only present records
        for FILE in $(ls *.vcf.gz); do
            SAMPLE_ID=$(basename ${FILE} .vcf.gz)
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'GT="alt"' --output-type z ${FILE} --output ${SAMPLE_ID}_extracted.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_extracted.vcf.gz
            rm -f ${FILE}*
        done
        ls -laht
        gcloud storage mv '*_extracted.vcf.gz*' ~{remote_outdir}/~{min_sv_length}bp/~{coverage_id}/~{caller_id}/
    >>>

    output {
    }

    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
