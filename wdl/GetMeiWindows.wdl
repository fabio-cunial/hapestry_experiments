version 1.0


# Builds BED files that contain all and only the hapestry windows of a given
# VCF that contain at least two MEI calls of the same type (this is to enrich
# for windows that likely contain similar INS at different locations).
#
workflow GetMeiWindows {
    input {
        File bcftools_merge_vcf_gz
        File bcftools_merge_tbi
        
        File tandems_bed
        File ref_fasta
        
        Int interval_max_length = 20000
        Int min_sv_length = 50
        Int flank_length = 200
    }
    parameter_meta {
    }
    
    call FindWindows {
        input:
            bcftools_merge_vcf_gz = bcftools_merge_vcf_gz,
            bcftools_merge_tbi = bcftools_merge_tbi,
            tandems_bed = tandems_bed,
            ref_fasta = ref_fasta,
            interval_max_length = interval_max_length,
            min_sv_length = min_sv_length,
            flank_length = flank_length
    }
    call GetMeis {
        input:
            windows_bed = FindWindows.windows_bed,
            bcftools_merge_vcf_gz = bcftools_merge_vcf_gz,
            bcftools_merge_tbi = bcftools_merge_tbi
    }
    
    output {
        File alu_bed = GetMeis.alu_bed
        File l1_bed = GetMeis.l1_bed
        File sva_bed = GetMeis.sva_bed
        File all_bed = GetMeis.all_bed
    }
}


# Just a call to hapestry's `find_windows` executable.
#
task FindWindows {
    input {
        File bcftools_merge_vcf_gz
        File bcftools_merge_tbi
        
        File tandems_bed
        File ref_fasta
        
        Int interval_max_length = 20000
        Int min_sv_length = 50
        Int flank_length = 200
        
        Int n_cpu = 2
        Int ram_gb = 8
        Int disk_size_gb = 50
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
        
        mv ~{bcftools_merge_vcf_gz} ./bcftools_merge.vcf.gz
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} --decompress bcftools_merge.vcf.gz
        mkdir ./out
        ${TIME_COMMAND} ~{docker_dir}/sv_merge/build/find_windows --n_chunks 1 --vcf bcftools_merge.vcf --tandems ~{tandems_bed} --interval_max_length ~{interval_max_length} --min_sv_length ~{min_sv_length} --flank_length ~{flank_length} --ref_fasta ~{ref_fasta} --output_dir ./out
    >>>
    output {
        File windows_bed = work_dir + "/out/run/windows_0_unflanked.bed"
    }
    runtime {
        docker: "fcunial/hapestry:merge"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


# Keeps the intervals of `windows_bed` that contain at least two calls annotated
# with the same SVANN value by PBSV (ALU, L1 or SVA).
#
task GetMeis {
    input {
        File windows_bed
        File bcftools_merge_vcf_gz
        File bcftools_merge_tbi
        
        Int n_cpu = 2
        Int ram_gb = 8
        Int disk_size_gb = 50
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
        
        
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'INFO/SVANN=="ALU"' --output-type z ~{bcftools_merge_vcf_gz} > alu.vcf.gz
        tabix -f alu.vcf.gz
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'INFO/SVANN=="L1"' --output-type z ~{bcftools_merge_vcf_gz} > l1.vcf.gz
        tabix -f l1.vcf.gz
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'INFO/SVANN=="SVA"' --output-type z ~{bcftools_merge_vcf_gz} > sva.vcf.gz
        tabix -f sva.vcf.gz
        date
        bedtools intersect -c -a ~{windows_bed} -b alu.vcf.gz | awk '{ if ($4>1) print $0; }' | sort -k1,1 -k2,2n > alu.bed
        date
        bedtools intersect -c -a ~{windows_bed} -b l1.vcf.gz | awk '{ if ($4>1) print $0; }' | sort -k1,1 -k2,2n > l1.bed
        date
        bedtools intersect -c -a ~{windows_bed} -b sva.vcf.gz | awk '{ if ($4>1) print $0; }' | sort -k1,1 -k2,2n > sva.bed
        date
        cat alu.bed l1.bed sva.bed | sort -k1,1 -k2,2n > tmp.bed
        ${TIME_COMMAND} bedtools merge -c 4 -o sum -i tmp.bed > all.bed
    >>>
    output {
        File alu_bed = work_dir + "/alu.bed"
        File l1_bed = work_dir + "/l1.bed"
        File sva_bed = work_dir + "/sva.bed"
        File all_bed = work_dir + "/all.bed"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
