version 1.0


#
workflow TruvariEvaluation {
    input {
        String sample_id
        String remote_input_dir
        String coverage_id
        String caller_id
        File sample_dipcall_vcf_gz
        File sample_dipcall_tbi
        
        Int min_sv_length
        String? region = "chr1"
        File confident_bed
        
        File reference_fa
        File reference_fai

        String docker_image = "fcunial/truvari_refine"
    }
    parameter_meta {
        sample_dipcall_vcf_gz: "Assumed to contain SNPs and short INDELs, but no multiallelics."
        confident_bed: "Might not be the dipcall confident BED of this sample, since hapestry was run on a fixed confident BED. This BED should be a subset of the BED that was input to hapestry."
    }
    
    call Impl { 
        input:
            sample_id = sample_id,
            remote_input_dir = remote_input_dir,
            coverage_id = coverage_id,
            caller_id = caller_id,
            sample_dipcall_vcf_gz = sample_dipcall_vcf_gz,
            sample_dipcall_tbi = sample_dipcall_tbi,
            min_sv_length = min_sv_length,
            region = region,
            confident_bed = confident_bed,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            docker_image = docker_image
    }
}


# Performance on an 8-core, 32GB VM, chr1 only, truvari+kanpig 10bp, 4x:
#
# COMMAND           CPU         RAM         TIME        COST
# 
#
task Impl {
    input {
        String sample_id
        String remote_input_dir
        String coverage_id
        String caller_id
        File sample_dipcall_vcf_gz
        File sample_dipcall_tbi
        
        Int min_sv_length
        String? region
        File confident_bed
        
        File reference_fa
        File reference_fai

        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 32
        Int disk_size_gb = 50
    }

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        TIME_COMMAND="/usr/bin/time --verbose"
        
        # Preparing the input VCFs
        ${TIME_COMMAND} gcloud storage cp ~{remote_input_dir}/~{min_sv_length}bp/~{coverage_id}/~{caller_id}/~{sample_id}_extracted.vcf.'gz*' .
        mv ~{sample_dipcall_vcf_gz} dipcall.vcf.gz
        mv ~{sample_dipcall_tbi} dipcall.vcf.gz.tbi
        mv ~{reference_fa} reference.fa
        mv ~{reference_fai} reference.fa.fai
        if ~{defined(region)}
        then
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --output-type z ~{sample_id}_extracted.vcf.gz ~{region} --output query.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t query.vcf.gz
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --output-type z dipcall.vcf.gz ~{region} --output truth.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t truth.vcf.gz
        else
            mv ~{sample_id}_extracted.vcf.gz query.vcf.gz
            mv ~{sample_id}_extracted.vcf.gz.tbi query.vcf.gz.tbi
            mv dipcall.vcf.gz truth.vcf.gz
            mv dipcall.vcf.gz.tbi truth.vcf.gz.tbi
        fi
        
        # Benchmarking
        ${TIME_COMMAND} truvari bench --sizemin ~{min_sv_length} --sizefilt ~{min_sv_length} --sizemax 10000 --includebed ~{confident_bed} --base truth.vcf.gz --comp query.vcf.gz --reference reference.fa --refine --output truvari_output/
        mv truvari_output/summary.json .
        mv truvari_output/candidate.refine.bed .
    >>>

    output {
        File summary_json = "summary.json"
        File refine_bed = "candidate.refine.bed"
    }

    runtime {
        docker: docker_image
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: ram_size_gb + " GiB"
        cpu: n_cpu
        preemptible: 0
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
