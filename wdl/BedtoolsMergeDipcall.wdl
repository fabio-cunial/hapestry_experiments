version 1.0


#
workflow BedtoolsMergeDipcall {
    input {
        Array[File] sample_bed
    }
    parameter_meta {
    }
    
    call BedtoolsMergeDipcallImpl {
        input:
            sample_bed = sample_bed
    }
    
    output {
        File union_bed = BedtoolsMergeDipcallImpl.union_bed
        File intersection_bed = BedtoolsMergeDipcallImpl.intersection_bed
        File stats_bed = BedtoolsMergeDipcallImpl.stats_bed
    }
}


task BedtoolsMergeDipcallImpl {
    input {
        Array[File] sample_bed
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    Int disk_size_gb = 10*ceil(size(sample_bed, "GB"))
    Int n_files = length(sample_bed)
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        INPUT_FILES=~{sep=' ' sample_bed}
        ${TIME_COMMAND} bedtools multiinter -header -i ${INPUT_FILES} > intersection_stats.bed
        awk '{ if ($4==~{n_files}) printf("%s\t%s\t%s\n",$1,$2,$3); }' intersection_stats.bed > intersection.bed
        ${TIME_COMMAND} bedtools merge -header -i ${INPUT_FILES} > union.bed
        ls -laht
    >>>
    
    output {
        File union_bed = work_dir + "/union.bed"
        File intersection_bed = work_dir + "/intersection.bed"
        File stats_bed = work_dir + "/intersection_stats.bed"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 2
        memory: "4GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
