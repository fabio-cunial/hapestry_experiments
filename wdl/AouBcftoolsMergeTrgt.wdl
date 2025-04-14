version 1.0


# - Cleans the raw, single-sample PBSV and TRGT calls in AoU. After this, let a
#   TRGT call be a TRGT record whose GT is different from 0/0 and ./.. Let a
#   PBSV call be a PBSV record with any GT.
# - For each sample, keeps only PBSV calls that overlap with the complement of
#   its TRGT calls. In particular, if a sample has no TRGT call in a TRGT
#   region, but has PBSV calls in that region, the PBSV calls are kept.
# - Merges all the resulting single-sample VCFs.
#
workflow AouBcftoolsMergeTrgt {
    input {
        Array[String] sample_ids
        Array[File] pbsv_vcf_gz
        Array[File] pbsv_tbi
        Array[File] trgt_vcf_gz
        Array[File] trgt_tbi
        
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    
    scatter (i in range(length(pbsv_vcf_gz))) {
        call CleanPbsv {
            input:
                sample_id = sample_ids[i],
                input_vcf_gz = pbsv_vcf_gz[i],
                input_tbi = pbsv_tbi[i],
                reference_fa = reference_fa,
                reference_fai = reference_fai
        }
        call CleanTrgt {
            input:
                sample_id = sample_ids[i],
                input_vcf_gz = trgt_vcf_gz[i],
                input_tbi = trgt_tbi[i]
        }
        call CombinePbsvTrgt {
            input:
                pbsv_vcf_gz = CleanPbsv.output_vcf_gz,
                pbsv_tbi = CleanPbsv.output_tbi,
                trgt_vcf_gz = CleanTrgt.output_vcf_gz,
                trgt_tbi = CleanTrgt.output_tbi,
                reference_fai = reference_fai
        }
    }
    call Merge {
        input:
            input_vcf_gz = CombinePbsvTrgt.output_vcf_gz,
            input_tbi = CombinePbsvTrgt.output_tbi
    }
    
    output {
        File output_vcf_gz = Merge.output_vcf_gz
        File output_tbi = Merge.output_tbi
    }
}


# Remark: all calls are kept, regardless of their GT.
#
# Remark: $bcftools merge$ collapses into the same record every record with the
# same $CHR,POS,REF,ALT$, disregarding the INFO field and in particular 
# differences in SVLEN and END. This may delete information for symbolic
# ALTs. Our script makes sure that only symbolic records with the same SVLEN and
# END are collapsed into the same record.
#
# Remark: we do not consider STRAND in the above.
#
# Remark: symbolic ALTs in the output VCF are not necessarily identical to the
# symbolic ALTs in the input. The original genotypes are preserved.
#
# Performance on each AoU 8x sample:
#
# COMMAND           RUNTIME     N_CPUS      MAX_RSS
# bcftools sort     4m          1           230M
# bcftools norm     60s         2           230M
# awk               30s         1           156M
# rm-dup            2s          1           16M
# bcftools merge    3m          2           250M
# bcftools norm     3m          2           230M
# bgzip             50s         2           9M
#
task CleanPbsv {
    input {
        String sample_id
        File input_vcf_gz
        File input_tbi
        File reference_fa
        File reference_fai
        
        Int ram_gb = 4
    }
    parameter_meta {
    }
    
    Int compression_level = 1
    Int disk_size_gb = 10*ceil(size(input_vcf_gz, "GB")) + 10
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_gb} - 2 ))
        
        function cleanVCF() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_VCF=$2
            
            # - Ensuring that the input file is sorted
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${INPUT_VCF_GZ} > tmp0.vcf.gz
            tabix -f tmp0.vcf.gz
            
            # - Removing multiallelic records.
            # - Fixing wrong REF values (which may occur e.g. in sniffles).
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z tmp0.vcf.gz > tmp1.vcf.gz
            tabix -f tmp1.vcf.gz
            rm -f tmp0.vcf.gz*
            
            # - Storing SVLEN and END in symbolic ALTs
            bcftools view --header-only tmp1.vcf.gz > tmp2.vcf
            ${TIME_COMMAND} bcftools view --no-header tmp1.vcf.gz | awk '{ \
                tag="artificial"; \
                if ($5=="<DEL>" || $5=="<INS>" || $5=="<INV>" || $5=="<DUP>" || $5=="<CNV>") { \
                    svtype=substr($5,2,3); \
                    end=""; \
                    svlen=""; \
                    n=split($8,A,";"); \
                    for (i=1; i<=n; i++) { \
                        if (substr(A[i],1,4)=="END=") end=substr(A[i],5); \
                        else if (substr(A[i],1,6)=="SVLEN=") { \
                            if (substr(A[i],7,1)=="-") svlen=substr(A[i],8); \
                            else svlen=substr(A[i],7); \
                        } \
                    } \
                    $5="<" svtype ":" tag ":" (length(end)==0?"?":end) ":" (length(svlen)==0?"?":svlen) ">" \
                }; \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            }' >> tmp2.vcf
            rm -f tmp1.vcf.gz*
            ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level ~{compression_level} tmp2.vcf
            tabix -f tmp2.vcf.gz
            
            # - Removing identical records
            # See <https://github.com/samtools/bcftools/issues/1089>.
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --rm-dup exact --output-type v tmp2.vcf.gz > ${OUTPUT_VCF}
            rm -f tmp2.vcf.gz*
            ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level ~{compression_level} ${OUTPUT_VCF}
            tabix -f ${OUTPUT_VCF}.gz
        }
        
        # - Cleaning the single-caller VCF
        cleanVCF ~{input_vcf_gz} cleaned.vcf
        
        # - Making sure the sample name is correct
        echo ~{sample_id} > samples.txt
        ${TIME_COMMAND} bcftools reheader --threads ${N_THREADS} --samples samples.txt --output out.vcf.gz cleaned.vcf.gz
        tabix -f out.vcf.gz
    >>>

    output {
        File output_vcf_gz = work_dir + "/out.vcf.gz"
        File output_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 2
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


# Removes records with 0/0 or ./. GT.
#
task CleanTrgt {
    input {
        String sample_id
        File input_vcf_gz
        File input_tbi
        
        Int ram_gb = 4
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(input_vcf_gz, "GB")) + 10
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_gb} - 2 ))
        
        # - Ensuring that the input file is sorted
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ~{input_vcf_gz} > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        
        # - Removing multiallelic records, which might be created by TRGT.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp1.vcf.gz > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        rm -f tmp1.vcf.gz*
        
        # - Discarding records with missing or ref GT
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --exclude 'GT="mis" || GT="0/0"' --output-type z tmp2.vcf.gz > tmp3.vcf.gz
        tabix -f tmp3.vcf.gz
        rm -f tmp2.vcf.gz
        
        # - Making sure the sample name is correct
        echo ~{sample_id} > samples.txt
        ${TIME_COMMAND} bcftools reheader --threads ${N_THREADS} --samples samples.txt --output out.vcf.gz tmp3.vcf.gz
        tabix -f out.vcf.gz
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/out.vcf.gz"
        File output_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 2
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


#
task CombinePbsvTrgt {
    input {
        File pbsv_vcf_gz
        File pbsv_tbi
        File trgt_vcf_gz
        File trgt_tbi
        
        File reference_fai
        
        Int n_cpu = 4
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 100
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # - Removing SVs inside TRGT intervals
        ${TIME_COMMAND} bedtools complement -i ~{trgt_vcf_gz} -g ~{reference_fai} > not_trgt.bed
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --regions-file not_trgt.bed --output-type z ~{pbsv_vcf_gz} > pbsv_cleaned.vcf.gz
        tabix -f pbsv_cleaned.vcf.gz
        
        # - Combining SV and TRGT calls
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --output-type z pbsv_cleaned.vcf.gz ~{trgt_vcf_gz} > out.vcf.gz
        tabix -f out.vcf.gz
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/out.vcf.gz"
        File output_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


# Performance on 1074 AoU 8x samples:
#
# COMMAND           RUNTIME     N_CPUS      MAX_RSS
# bcftools merge    
#
task Merge {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
        
        Int ram_gb = 256
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil(size(input_vcf_gz, "GB")) + 100
    Int compression_level = 1
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
        
        INPUT_FILES=~{sep=',' input_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            echo ${INPUT_FILE} >> list.txt
        done
        
        # - Merging all samples
        # $--info-rules -$ disables default rules, and it is used just to avoid
        # the following error:
        # Only fixed-length vectors are supported with -i sum:AC
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --info-rules - --file-list list.txt --output-type z > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        ls -laht
        
        # - Removing multiallelic records, if any are generated during the
        # merge. This is just an extra safeguard and might be dropped if
        # $bcftools merge --merge none$ always behaves correctly.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp1.vcf.gz > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        rm -f tmp1.vcf.gz*
        ls -laht
        
        # - Restoring symbolic ALTs to their original states
        bcftools view --header-only tmp2.vcf.gz > merged.vcf
        ${TIME_COMMAND} bcftools view --no-header tmp2.vcf.gz | awk '{ \
            tag="artificial"; \
            if (substr($5,6,length(tag))==tag) $5=substr($5,1,4) ">"; \
            printf("%s",$1); \
            for (i=2; i<=NF; i++) printf("\t%s",$i); \
            printf("\n"); \
        }' >> merged.vcf
        rm -f tmp2.vcf.gz*
        ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level ~{compression_level} merged.vcf
        tabix -f merged.vcf.gz
        ls -laht
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/merged.vcf.gz"
        File output_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 4
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}









#
task Concat {
    input {
        File sv_merge_vcf_gz
        File sv_merge_tbi
        File trgt_merge_vcf_gz
        File trgt_merge_tbi
        
        File samples_file
        File reference_fai
        
        Int n_cpu = 4
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil( size(sv_merge_vcf_gz, "GB") + size(trgt_merge_vcf_gz, "GB") ) + 100
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
        
        # Ensuring that samples have the same order in both files
        ${TIME_COMMAND} bcftools view --samples-file ~{samples_file} --output-type z ~{sv_merge_vcf_gz} > sv.vcf.gz
        tabix -f sv.vcf.gz
        rm -f ~{sv_merge_vcf_gz}
        ${TIME_COMMAND} bcftools view --samples-file ~{samples_file} --output-type z ~{trgt_merge_vcf_gz} > trgt.vcf.gz
        tabix -f trgt.vcf.gz
        rm -f ~{trgt_merge_vcf_gz}
        
        # Removing SVs inside TRGT intervals
        ${TIME_COMMAND} bedtools complement -i trgt.vcf.gz -g ~{reference_fai} > not_trgt.bed
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --regions-file not_trgt.bed --output-type z sv.vcf.gz > sv_cleaned.vcf.gz
        tabix -f sv_cleaned.vcf.gz
        rm -f sv.vcf.gz*
        
        # Combining SV and TRGT calls
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --output-type z sv_cleaned.vcf.gz trgt.vcf.gz > out.vcf.gz
        tabix -f out.vcf.gz
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/out.vcf.gz"
        File output_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}