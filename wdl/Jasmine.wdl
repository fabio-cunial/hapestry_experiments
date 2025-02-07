version 1.0

# In intra-sample jasmine we use:
#
# --allow_intrasample min_seq_id=0.9
#
# In particular, `--allow_intrasample` is necessary, since the bcftools merge
# VCF in input contains just one sample column and variants from different
# callers. `min_seq_id=0.9` mimics what is done in a production pipeline like
# AoU, but in jasmine it only works with INS.
#
# When we do inter-sample jasmine on the bcftools merge of all raw calls from
# all callers, we use:
#
# --allow_intrasample min_seq_id=0.7
#
# Once again, `--allow_intrasample` is needed, since variants from the same
# sample might come from different callers, and `min_seq_id=0.7` mimics the AoU
# setting above, even though this is a different input VCF.
#
workflow Jasmine {
    input {
        String sample_id
        File bcftools_merge_vcf_gz
        File bcftools_merge_tbi
        String jasmine_params = " "
        Int max_sv_length = 0
        Int n_cpu
        Int ram_gb
        Int disk_size_gb
    }
    parameter_meta {
        max_sv_length: "SVs that are too long make working on the inter-sample bcftools merge too slow (50000 works well in this case). 0=do not remove long SVs."
    }
    
    call JasmineImpl {
        input:
            sample_id = sample_id,
            bcftools_merge_vcf_gz = bcftools_merge_vcf_gz,
            bcftools_merge_tbi = bcftools_merge_tbi,
            jasmine_params = jasmine_params,
            max_sv_length = max_sv_length,
            n_cpu = n_cpu,
            ram_gb = ram_gb,
            disk_size_gb = disk_size_gb
    }
    output {
        File output_vcf_gz = JasmineImpl.output_vcf_gz
        File output_tbi = JasmineImpl.output_tbi
    }
}


# Performance on a VM with 16 cores and 16GB of RAM:
#
# COVERAGE  CPU     RAM     TIME
# 4x        100%    2.5G    9m
# 8x        150%    1.9G    2m
# 16x       150%    1.8G    2m
# 32x       150%    2.2G    2m
#
task JasmineImpl {
    input {
        String sample_id
        File bcftools_merge_vcf_gz
        File bcftools_merge_tbi
        String jasmine_params
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
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        EFFECTIVE_MEM_GB=$(( ~{ram_gb} - 2 ))
        JAVA_PATH="/usr/bin/java"  # Using the latest JRE
        
        # - Removing calls that are too long. Without this, processing the
        #   bcftools merge of all samples takes days.
        if [ ~{max_sv_length} -eq 0 ]; then
            gunzip -c ~{bcftools_merge_vcf_gz} > tmp1.vcf
        else
            bcftools filter --include "INFO/SVLEN>=-~{max_sv_length} && INFO/SVLEN<=~{max_sv_length}" --output-type v ~{bcftools_merge_vcf_gz} > tmp1.vcf
        fi
        rm -f ~{bcftools_merge_vcf_gz}
        
        # - Making all DELs and INVs symbolic to speed up Jasmine. Without this,
        #   in degenerate cases it takes hours to intra-sample merge on a
        #   reasonably-sized cloud VM, even with SSD.
        python ~{docker_dir}/symbolic_jasmine.py tmp1.vcf > input.vcf
        rm -f tmp1.vcf
        echo "input.vcf" > list.txt
        
        # - Using `--output_genotypes` on a bcftools merge VCF with only one
        #   sample leads to a NullPointerException:
        #   at AddGenotypes.addGenotypes(AddGenotypes.java:152).
        ${TIME_COMMAND} ${JAVA_PATH} -jar /opt/conda/bin/jasmine.jar -Xms${EFFECTIVE_MEM_GB}G -Xmx${EFFECTIVE_MEM_GB}G threads=${N_THREADS} ~{jasmine_params} file_list=list.txt out_file=tmp2.vcf
        
        # - Removing a suffix of the INFO field added by Jasmine, since it makes
        #   bcftools sort crash.
        bcftools view --header-only --no-version tmp2.vcf > tmp3.vcf
        N_ROWS=$(wc -l < tmp3.vcf)
        tail -n +$(( ${N_ROWS} + 1 )) tmp2.vcf > body.txt
        rm -f tmp2.vcf
        ${TIME_COMMAND} cat body.txt | awk '{ \
            pattern="ALLVARS_EXT"; \
            \
            printf("%s",$1); \
            for (i=2; i<8; i++) printf("\t%s",$i); \
            i=match($8,pattern); \
            if (i!=0) printf("\t%s",substr($8,1,i-1)); \
            else printf("\t%s",$8); \
            for (i=9; i<=NF; i++) printf("\t%s",$i); \
            printf("\n"); \
        }' >> tmp3.vcf
        rm -f body.txt
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp3.vcf > ~{sample_id}.jasmine.vcf.gz
        tabix -f ~{sample_id}.jasmine.vcf.gz
    >>>
    output {
        File output_vcf_gz = work_dir + "/" + sample_id + ".jasmine.vcf.gz"
        File output_tbi = work_dir + "/" + sample_id + ".jasmine.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
