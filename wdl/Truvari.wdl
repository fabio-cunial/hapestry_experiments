version 1.0

# In intra-sample truvari we use:
#
# --sizemin 0 --sizemax 1000000 --pctseq 0.90 --pctsize 0.90
#
# In particular, we keep `--gt` to default (off), since this allows collapsing
# variants from the same sample. This is necessary, since the bcftools merge VCF
# in input contains just one sample column and variants from different callers.
# `--pctseq` and `--pctsize` are as suggested by Adam in AoU.
#
# When we do inter-sample truvari on the VCFs emitted by the intra-sample
# truvari above, we use:
#
# --sizemin 0 --sizemax 1000000 --keep common --gt all
#
# In particular, `--gt all` disables merging calls from the same sample, since
# we assume that intra-sample merging has already been performed. `--pctseq`
# and `--pctsize` are left to default (0.7), as suggested by Adam in AoU.
#
# When we do inter-sample truvari on the bcftools merge of all raw calls from
# all callers, we use:
#
# --sizemin 0 --sizemax 1000000 --keep common
#
# Once again, we keep `--gt` to default (off) to collapse variants from the
# same sample, which might come from different callers. `--pctseq` and
# `--pctsize` are left to default (0.7), to mimic the AoU setting above, even
# though this is a different input VCF.
#
workflow Truvari {
    input {
        String sample_id
        File bcftools_merge_vcf_gz
        File bcftools_merge_tbi
        String truvari_args = " "
        Int ram_gb
    }
    parameter_meta {
    }
    
    call TruvariImpl {
        input:
            sample_id = sample_id,
            bcftools_merge_vcf_gz = bcftools_merge_vcf_gz,
            bcftools_merge_tbi = bcftools_merge_tbi,
            truvari_args = truvari_args,
            ram_gb = ram_gb
    }
    
    output {
    	File vcf_gz = TruvariImpl.vcf_gz
    	File tbi = TruvariImpl.tbi
    }
}


#
task TruvariImpl {
    input {
        String sample_id
        File bcftools_merge_vcf_gz
        File bcftools_merge_tbi
        String truvari_args
        Int ram_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(bcftools_merge_vcf_gz,"GB")) + 100
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_MEM_GB=$(( ~{ram_gb} - 2 ))

        ${TIME_COMMAND} truvari collapse --input ~{bcftools_merge_vcf_gz} ~{truvari_args} --output tmp.vcf
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp.vcf > ~{sample_id}.truvari_collapsed.vcf.gz
        tabix -f ~{sample_id}.truvari_collapsed.vcf.gz
    >>>
    
    output {
    	File vcf_gz = work_dir + "/" + sample_id + ".truvari_collapsed.vcf.gz"
        File tbi = work_dir + "/" + sample_id + ".truvari_collapsed.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 1
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
