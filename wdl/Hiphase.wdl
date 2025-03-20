version 1.0


#
workflow Hiphase {
    input {
        String samplename
        File bam
        File bai
        File unphased_snp_vcf
        File unphased_snp_tbi
        File unphased_sv_vcf
        File unphased_sv_tbi
        File ref_fasta
        File ref_fasta_fai
        Int n_cpu
        Int ram_gb
        String extra_args
    }
    parameter_meta {
    }
    
    call HiphaseImpl {
        input:
            samplename = samplename,
            bam = bam,
            bai = bai,
            unphased_snp_vcf = unphased_snp_vcf,
            unphased_snp_tbi = unphased_snp_tbi,
            unphased_sv_vcf = unphased_sv_vcf,
            unphased_sv_tbi = unphased_sv_tbi,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            n_cpu = n_cpu,
            ram_gb = ram_gb,
            extra_args = extra_args
    }
    
    output {
        File phased_snp_vcf = HiphaseImpl.phased_snp_vcf
        File phased_snp_vcf_tbi = HiphaseImpl.phased_snp_vcf_tbi
        File phased_sv_vcf = HiphaseImpl.phased_sv_vcf
        File phased_sv_vcf_tbi = HiphaseImpl.phased_sv_vcf_tbi
        File merged_vcf = HiphaseImpl.merged_vcf
        File merged_vcf_tbi = HiphaseImpl.merged_vcf_tbi
        File haplotag_file = HiphaseImpl.haplotag_file
        File stats_file = HiphaseImpl.stats_file
        File blocks_file = HiphaseImpl.blocks_file
        File summary_file = HiphaseImpl.summary_file
    }
}


#
task HiphaseImpl {
    input {
        String samplename
        File bam
        File bai
        File unphased_snp_vcf
        File unphased_snp_tbi
        File unphased_sv_vcf
        File unphased_sv_tbi
        File ref_fasta
        File ref_fasta_fai
        Int n_cpu
        Int ram_gb
        String extra_args
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    Int disk_size_gb = 10*ceil(size(unphased_snp_vcf,"GB") + size(unphased_sv_vcf,"GB") + size(ref_fasta,"GB") + size(bam,"GB")) + 100
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} ~{docker_dir}/hiphase --threads ${N_THREADS} --global-realignment-cputime 300 \
            --reference ~{ref_fasta} \
            --bam ~{bam} \
            --vcf ~{unphased_snp_vcf} \
            --output-vcf ~{samplename}_phased_snp.vcf.gz \
            --vcf ~{unphased_sv_vcf} \
            --output-vcf ~{samplename}_phased_sv.vcf.gz \
            --haplotag-file ~{samplename}_phased_sv_haplotag.tsv \
            --stats-file ~{samplename}.stats.csv \
            --blocks-file ~{samplename}.blocks.tsv \
            --summary-file ~{samplename}.summary.tsv \
            ~{extra_args}
            
        # Sorting and merging
        ${TIME_COMMAND} bcftools sort ~{samplename}_phased_snp.vcf.gz -O z -o ~{samplename}_phased_snp.sorted.vcf.gz
        tabix -p vcf ~{samplename}_phased_snp.sorted.vcf.gz
        ${TIME_COMMAND} bcftools sort ~{samplename}_phased_sv.vcf.gz -O z -o ~{samplename}_phased_sv.sorted.vcf.gz
        tabix -p vcf ~{samplename}_phased_sv.sorted.vcf.gz
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS}--allow-overlaps --remove-duplicates --output-type z ~{samplename}_phased_snp.sorted.vcf.gz ~{samplename}_phased_sv.sorted.vcf.gz > ~{samplename}_merged.vcf.gz
        tabix -f ~{samplename}_merged.vcf.gz
    >>>
    
    output {
        File phased_snp_vcf =     work_dir + "/" + samplename + "_phased_snp.sorted.vcf.gz"
        File phased_snp_vcf_tbi = work_dir + "/" + samplename + "_phased_snp.sorted.vcf.gz.tbi"
        File phased_sv_vcf =      work_dir + "/" + samplename + "_phased_sv.sorted.vcf.gz"
        File phased_sv_vcf_tbi =  work_dir + "/" + samplename + "_phased_sv.sorted.vcf.gz.tbi"
        File merged_vcf =         work_dir + "/" + samplename + "_merged.vcf.gz"
        File merged_vcf_tbi =     work_dir + "/" + samplename + "_merged.vcf.gz.tbi"
        File haplotag_file =      work_dir + "/" + samplename + "_phased_sv_haplotag.tsv"
        File stats_file =         work_dir + "/" + samplename + ".stats.csv"
        File blocks_file =        work_dir + "/" + samplename + ".blocks.tsv"
        File summary_file =       work_dir + "/" + samplename + ".summary.tsv"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
