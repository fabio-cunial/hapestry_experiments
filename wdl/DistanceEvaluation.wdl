version 1.0

#
workflow DistanceEvaluation {
    input {
        String min_sv_length
        String coverage
        String sample_id
        
        File flanked_windows_bed
        File confident_bed
        
        File assembly_bam_csv
        File hapestry_bam_csv
        File kanpig_bam_csv
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            min_sv_length = min_sv_length,
            coverage = coverage,
            sample_id = sample_id,
            flanked_windows_bed = flanked_windows_bed,
            confident_bed = confident_bed,
            assembly_bam_csv = assembly_bam_csv,
            hapestry_bam_csv = hapestry_bam_csv,
            kanpig_bam_csv = kanpig_bam_csv
    }
    
    output {
        File alignment_distances_hapestry = Impl.alignment_distances_hapestry
        File alignment_distances_kanpig = Impl.alignment_distances_kanpig
        File alignment_distances_hapestry_kanpig_csv = Impl.alignment_distances_hapestry_kanpig_csv
    }
}


# Performance on a VM with 8 cores and 16 GB of RAM (chr1 only):
#
# TOOL                                      CPU%     RAM      TIME
# extract_reads_from_windows
# distance_evaluation_align_windows.py
# sort
# join
#
task Impl {
    input {
        String min_sv_length
        String coverage
        String sample_id
        
        File flanked_windows_bed
        File confident_bed
        
        File assembly_bam_csv
        File hapestry_bam_csv
        File kanpig_bam_csv
        
        Int ram_gb = 16
        Int n_cpu = 8
        Int disk_size_gb = 50
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        FETCH_MAX_LENGTH="100000"  # Arbitrary
        
        # Downloading windows
        ${TIME_COMMAND} /hapestry/sv_merge/build/extract_reads_from_windows --n_threads ${N_THREADS} --bam_csv ~{assembly_bam_csv} --output_dir ./windows_assembly/ --windows ~{flanked_windows_bed} --flank_length 0 --fetch_max_length ${FETCH_MAX_LENGTH} --require_spanning
        ${TIME_COMMAND} /hapestry/sv_merge/build/extract_reads_from_windows --n_threads ${N_THREADS} --bam_csv ~{hapestry_bam_csv} --output_dir ./windows_hapestry/ --windows ~{flanked_windows_bed} --flank_length 0 --fetch_max_length ${FETCH_MAX_LENGTH} --require_spanning
        ${TIME_COMMAND} /hapestry/sv_merge/build/extract_reads_from_windows --n_threads ${N_THREADS} --bam_csv ~{kanpig_bam_csv}   --output_dir ./windows_kanpig/   --windows ~{flanked_windows_bed} --flank_length 0 --fetch_max_length ${FETCH_MAX_LENGTH} --require_spanning
        
        # Aligning windows
        ${TIME_COMMAND} python distance_evaluation_align_windows.py ./windows_assembly/ ./windows_hapestry/ ./alignment_distances_hapestry.csv --confident-bed ~{confident_bed}
        ${TIME_COMMAND} python distance_evaluation_align_windows.py ./windows_assembly/ ./windows_kanpig/   ./alignment_distances_kanpig.csv   --confident-bed ~{confident_bed}
        ${TIME_COMMAND} sort -t , -k 1,1 alignment_distances_hapestry.csv > alignment_distances_hapestry_sorted.csv
        ${TIME_COMMAND} sort -t , -k 1,1 alignment_distances_kanpig.csv > alignment_distances_kanpig_sorted.csv
        ${TIME_COMMAND} join -t ',' -1 1 -2 1 alignment_distances_hapestry_sorted.csv alignment_distances_kanpig_sorted.csv > join.csv
        ls -laht
        
        # Outputting
        mv join.csv ~{min_sv_length}bp_~{coverage}_~{sample_id}_alignment_distances_hapestry_kanpig.csv
        mv alignment_distances_hapestry_sorted.csv ~{min_sv_length}bp_~{coverage}_~{sample_id}_alignment_distances_hapestry.csv
        mv alignment_distances_kanpig_sorted.csv ~{min_sv_length}bp_~{coverage}_~{sample_id}_alignment_distances_kanpig.csv
    >>>
    
    output {
        File alignment_distances_hapestry = min_sv_length + "bp_" + coverage + "_" + sample_id + "_alignment_distances_hapestry.csv"
        File alignment_distances_kanpig = min_sv_length + "bp_" + coverage + "_" + sample_id + "_alignment_distances_kanpig.csv"
        File alignment_distances_hapestry_kanpig_csv = min_sv_length + "bp_" + coverage + "_" + sample_id + "_alignment_distances_hapestry_kanpig.csv"
    }
    runtime {
        docker: "fcunial/hapestry_distance_evaluation"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
