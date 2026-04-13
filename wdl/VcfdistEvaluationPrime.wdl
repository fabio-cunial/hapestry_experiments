version 1.0


#
workflow VcfdistEvaluationPrime {
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
        Int vcfdist_mode
        
        File reference_fa
        File reference_fai

        String docker_image = "fcunial/hapestry_experiments"
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
            vcfdist_mode = vcfdist_mode,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            docker_image = docker_image
    }
}


# Performance on an 8-core, 32GB VM, chr1 only, truvari+kanpig 10bp, 4x:
#
# COMMAND           CPU         RAM         TIME        COST
# vcfdist           200%        20G         4h          $2
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
        Int vcfdist_mode
        
        File reference_fa
        File reference_fai

        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 24
        Int disk_size_gb = 50
        Int preemptible_number = 0
        String time_command = "/usr/bin/time --verbose"
    }

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        
        # Preparing the input VCFs
        ~{time_command} gcloud storage cp ~{remote_input_dir}/~{min_sv_length}bp/~{coverage_id}/~{caller_id}/~{sample_id}_extracted.vcf.'gz*' .
        mv ~{sample_dipcall_vcf_gz} dipcall.vcf.gz
        mv ~{sample_dipcall_tbi} dipcall.vcf.gz.tbi
        mv ~{reference_fa} reference.fa
        mv ~{reference_fai} reference.fa.fai
        if ~{defined(region)}
        then
            ~{time_command} bcftools view --threads ${N_THREADS} --output-type z ~{sample_id}_extracted.vcf.gz ~{region} --output query.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t query.vcf.gz
            ~{time_command} bcftools view --threads ${N_THREADS} --output-type z dipcall.vcf.gz ~{region} --output truth.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t truth.vcf.gz
        else
            mv ~{sample_id}_extracted.vcf.gz query.vcf.gz
            mv ~{sample_id}_extracted.vcf.gz.tbi query.vcf.gz.tbi
            mv dipcall.vcf.gz truth.vcf.gz
            mv dipcall.vcf.gz.tbi truth.vcf.gz.tbi
        fi
        
        # See https://github.com/TimD1/vcfdist/wiki/02-Parameters-and-Usage
        if [ ~{vcfdist_mode} -eq 0 ]; then
            # Default
            ~{time_command} vcfdist query.vcf.gz truth.vcf.gz reference.fa \
                --max-threads ${N_THREADS} --max-ram ${EFFECTIVE_RAM_GB} --verbosity 1 \
                --realign-query --realign-truth \
                --bed ~{confident_bed} \
                --distance \
                --prefix ~{sample_id}_
        elif [ ~{vcfdist_mode} -eq 1 ]; then        
            # Default plus --sv-threshold
            ~{time_command} vcfdist query.vcf.gz truth.vcf.gz reference.fa \
                --sv-threshold ~{min_sv_length} \
                \
                --max-threads ${N_THREADS} --max-ram ${EFFECTIVE_RAM_GB} --verbosity 1 \
                --realign-query --realign-truth \
                --bed ~{confident_bed} \
                --distance \
                --prefix ~{sample_id}_
        elif [ ~{vcfdist_mode} -eq 2 ]; then
            # Default plus --largest-variant
            # Remark: `--max-supercluster-size` has to be >= `--largest-variant
            # + 2`
            ~{time_command} vcfdist query.vcf.gz truth.vcf.gz reference.fa \
                --largest-variant 10000 \
                --max-supercluster-size 10002 \
                \        
                --max-threads ${N_THREADS} --max-ram ${EFFECTIVE_RAM_GB} --verbosity 1 \
                --realign-query --realign-truth \
                --bed ~{confident_bed} \
                --distance \
                --prefix ~{sample_id}_
        elif [ ~{vcfdist_mode} -eq 3 ]; then
            # Default plus --sv-threshold plus --largest-variant
            ~{time_command} vcfdist query.vcf.gz truth.vcf.gz reference.fa \
                --sv-threshold ~{min_sv_length} \
                \
                --largest-variant 10000 \
                --max-supercluster-size 10002 \
                \
                --max-threads ${N_THREADS} --max-ram ${EFFECTIVE_RAM_GB} --verbosity 1 \
                --realign-query --realign-truth \
                --bed ~{confident_bed} \
                --distance \
                --prefix ~{sample_id}_
        fi
        ls -laht 1>&2
        df -h 1>&2
        ${TIME_COMMAND} tar -czf ~{min_sv_length}bp_~{coverage_id}_~{caller_id}_~{sample_id}_vcfdist.tar.gz ~{sample_id}_*
    >>>

    output {
        File out_tar_gz = min_sv_length + "bp_" + coverage_id + "_" + caller_id + "_" + sample_id + "_vcfdist.tar.gz"
    }

    runtime {
        docker: docker_image
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: ram_size_gb + " GiB"
        cpu: n_cpu
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
