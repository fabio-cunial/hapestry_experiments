version 1.0


# Builds a BED file that contains all and only the intervals in the PBSV tandem
# repeat track that overlap with more than one interval in the TRGT database
# (this is to enrich for intervals that contain sections with different periods,
# since the TRGT database is very strict in defining the boundaries of
# intervals), and that also overlap with at least one SV call from a given
# callset.
#
workflow GetCompositeTRs {
    input {
        File tandem_repeat_track_bed
        File trgt_repeat_catalog
        File bcftools_merge_vcf_gz
        File bcftools_merge_tbi
        
        Int ram_gb = 32
        Int disk_size_gb = 50
    }
    parameter_meta {
    }
    
    call GetCompositeTRsImpl {
        input:
            tandem_repeat_track_bed = tandem_repeat_track_bed,
            trgt_repeat_catalog = trgt_repeat_catalog,
            bcftools_merge_vcf_gz = bcftools_merge_vcf_gz,
            bcftools_merge_tbi = bcftools_merge_tbi,
            ram_gb = ram_gb,
            disk_size_gb = disk_size_gb
    }
    output {
        File composite_trs_bed = GetCompositeTRsImpl.composite_trs_bed
    }
}

 
#
task GetCompositeTRsImpl {
    input {
        File tandem_repeat_track_bed
        File trgt_repeat_catalog
        File bcftools_merge_vcf_gz
        File bcftools_merge_tbi
        
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
        
        ${TIME_COMMAND} bedtools intersect -c -a ~{tandem_repeat_track_bed} -b ~{trgt_repeat_catalog} | awk '{ if ($4>1) print $0; }' > tmp1.bed
        ${TIME_COMMAND} bedtools intersect -u -a tmp1.bed -b ~{bcftools_merge_vcf_gz} > composite_trs.bed
    >>>
    output {
        File composite_trs_bed = work_dir + "/composite_trs.bed"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 1
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
