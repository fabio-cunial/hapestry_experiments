version 1.0


# Re-genotypes an inter-sample VCF with kanpig, and keeps only records that
# occur in at least one sample after re-genotyping.
#
workflow KanpigGenotypeAndFilter {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        Array[Boolean] is_male
        Array[File] alignments_bam
        Array[File] alignments_bai
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
    }
    parameter_meta {
    }
    
    call RemoveSamples {
        input:
            intersample_vcf_gz = intersample_vcf_gz,
            intersample_tbi = intersample_tbi
    }
    scatter(i in range(length(alignments_bam))) {
        call Kanpig {
            input:
                is_male = is_male[i],
                input_vcf_gz = RemoveSamples.cleaned_vcf_gz,
                input_tbi = RemoveSamples.cleaned_tbi,
                alignments_bam = alignments_bam[i],
                alignments_bai = alignments_bai[i],
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                ploidy_bed_female = ploidy_bed_female,
                ploidy_bed_male = ploidy_bed_male
        }
    }
    call Merge {
        input:
            regenotyped_vcfs = Kanpig.regenotyped_kanpig
    }
    output {
        File output_vcf_gz = Merge.output_vcf_gz
        File output_tbi = Merge.output_tbi
    }
}


# Replaces all sample columns with a single sample column where all calls have
# 0/1 GT. Assigns a unique integer ID to every call, moving the original ID to
# the INFO field (this is necessary for kanpig, which otherwise complains
# about duplicated IDs).
#
task RemoveSamples {
    input {
        File intersample_vcf_gz
        File intersample_tbi
    }
    parameter_meta {
    }
    
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
        
        
        bcftools view --header-only ~{intersample_vcf_gz} > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > cleaned.vcf
        echo '##INFO=<ID=ORIGINAL_ID,Number=1,Type=String,Description="Original ID">' >> cleaned.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> cleaned.vcf
        bcftools view --no-header ~{intersample_vcf_gz} | awk 'BEGIN { i=0; } { gsub(/;/,"_",$3); printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s;ORIGINAL_ID=%s\tGT\t0/1\n",$1,$2,++i,$4,$5,$6,$7,$8,$3); }' >> cleaned.vcf
        rm -f ~{intersample_vcf_gz}
        bgzip -@ ${N_THREADS} cleaned.vcf
        tabix -f cleaned.vcf.gz
    >>>

    output {
        File cleaned_vcf_gz = work_dir + "/cleaned.vcf.gz"
        File cleaned_tbi = work_dir + "/cleaned.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 1
        memory: "32GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}


# Parameters are tuned for a multi-sample VCF.
#
task Kanpig {
    input {
        Boolean is_male
        File input_vcf_gz
        File input_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        Int n_cpu = 4
        Int ram_size_gb = 32
        File ploidy_bed_female
        File ploidy_bed_male
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    String output_prefix = "kanpig_regenotyped"
    String kanpig_params_multisample =  "--sizemin 50 --sizemax 10000 --neighdist  500 --gpenalty 0.04 --hapsim 0.97"
    Int disk_size_gb = 200 + ceil(size(reference_fa,"GB")) + 100*ceil(size(input_vcf_gz,"GB")) + 2*ceil(size(alignments_bam,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_MEM_GB=~{ram_size_gb}
        EFFECTIVE_MEM_GB=$(( ${EFFECTIVE_MEM_GB} - 4 ))
        df -h
        
        if [ ~{is_male} == true ]; then
            PLOIDY_BED=$(echo ~{ploidy_bed_male})
        else
            PLOIDY_BED=$(echo ~{ploidy_bed_female})
        fi
        export RUST_BACKTRACE="full"
        ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ~{kanpig_params_multisample} --reference ~{reference_fa} --input ~{input_vcf_gz} --reads ~{alignments_bam} --out tmp1.vcf.gz
        bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp1.vcf.gz > ~{output_prefix}.vcf.gz
        tabix -f ~{output_prefix}.vcf.gz
    >>>

    output {
        File regenotyped_kanpig = work_dir + "/" + output_prefix + ".vcf.gz"
        File regenotyped_kanpig_tbi = work_dir + "/" + output_prefix + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


# Assume that we have built an inter-sample-merged VCF, that we have
# subsequently kept just one column of it (setting it to all 0/1), and that we
# have re-genotyped this single-column VCF using the BAM of every original
# sample. The program combines all such re-genotyped single-sample VCFs into
# an inter-sample-merged VCF with updated genotypes.
#
# Remark: every file in $regenotyped_vcfs$ is assumed to have exactly the
# same set of calls in the same order. This is usually the case if the files
# are the result of re-genotyping the same inter-sample VCF with different BAMs.
# If this is not the case the program gives a wrong output, and it should be
# replaced with $bcftools merge --merge none$.
#
task Merge {
    input {
        Array[File] regenotyped_vcfs
        Int n_cpus = 8
    }
    parameter_meta {
    }
    
    String work_dir = "/cromwell_root/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        function pasteThread() {
            local THREAD_ID=$1
            
            OUTPUT_FILE="columns_${THREAD_ID}.txt"; touch ${OUTPUT_FILE}
            FIELDS_FILE="fields_${THREAD_ID}.txt"; touch ${FIELDS_FILE}
            TMP_PREFIX="tmp_${THREAD_ID}"
            i="0";
            while read ADDRESS; do
                i=$(( $i + 1 ))
                # Adding the new sample to the set of columns
                cp ${ADDRESS} ${TMP_PREFIX}.vcf.gz
                bcftools view --header-only ${TMP_PREFIX}.vcf.gz > ${TMP_PREFIX}.txt
                N_ROWS=$(wc -l < ${TMP_PREFIX}.txt)
                tail -n 1 ${TMP_PREFIX}.txt | cut -f 10 > sample_${THREAD_ID}.txt
                if [ $i = "1" ]; then
                    mv sample_${THREAD_ID}.txt ${FIELDS_FILE}
                else
                    paste ${FIELDS_FILE} sample_${THREAD_ID}.txt > ${FIELDS_FILE}.prime
                    mv ${FIELDS_FILE}.prime ${FIELDS_FILE}
                fi
                echo "Current fields of thread ${THREAD_ID}:"; cat ${FIELDS_FILE}
                # Adding the new column to the body
                bcftools view --no-header ${TMP_PREFIX}.vcf.gz | cut -f 10 > ${TMP_PREFIX}.txt
                if [ $i = "1" ]; then
                    mv ${TMP_PREFIX}.txt ${OUTPUT_FILE}
                else
                    ${TIME_COMMAND} paste ${OUTPUT_FILE} ${TMP_PREFIX}.txt > ${OUTPUT_FILE}.prime
                    mv ${OUTPUT_FILE}.prime ${OUTPUT_FILE}
                fi
            done < list_${THREAD_ID}
        }
        
        
        # Main program
        INPUT_FILES=~{sep=',' regenotyped_vcfs}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            echo ${INPUT_FILE} >> list.txt
        done
        
        # Initializing the inter-sample VCF with the first file
        ADDRESS=$(head -n 1 list.txt)
        mv ${ADDRESS} first.vcf.gz
        bcftools view --header-only first.vcf.gz > tmp.txt
        N_ROWS=$(wc -l < tmp.txt)
        head -n $(( ${N_ROWS} - 1 )) tmp.txt > header.txt
        tail -n 1 tmp.txt | cut -f 1,2,3,4,5,6,7,8,9 > fields.txt
        bcftools view --no-header first.vcf.gz | cut -f 1,2,3,4,5,6,7,8,9 > calls.txt
        rm -f first.vcf.gz
        
        # Appending all the remaining files
        N_ROWS=$(wc -l < list.txt)
        N_ROWS=$(( ${N_ROWS} / ${N_THREADS} ))
        split -d -l ${N_ROWS} list.txt list_
        COLUMNS_FILES=""; FIELDS_FILES=""
        for LIST_FILE in $(find . -maxdepth 1 -name 'list_*' | sort); do
            ID=${LIST_FILE#./list_}
            pasteThread ${ID} &
            COLUMNS_FILES="${COLUMNS_FILES} columns_${ID}.txt"
            FIELDS_FILES="${FIELDS_FILES} fields_${ID}.txt"
        done
        wait
        paste fields.txt ${FIELDS_FILES} > fields_all.txt
        paste calls.txt ${COLUMNS_FILES} > body.txt
        cat header.txt fields_all.txt body.txt > merged.vcf
        rm -f *.txt
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} merged.vcf
        tabix -f merged.vcf.gz
        
        # Keeping only the records that occur at least once in some sample
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="0/1" || GT="0|1" || GT="1/0" || GT="1|0" || GT="1/1" || GT="1|1")>0' --output-type z merged.vcf.gz > filtered.vcf.gz
        tabix -f filtered.vcf.gz
    >>>

    output {
        File output_vcf_gz = work_dir + "/filtered.vcf.gz"
        File output_tbi = work_dir + "/filtered.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cpus
        disks: "local-disk 128 HDD"
        preemptible: 0
    }
}
