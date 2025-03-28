version 1.0


# Re-genotypes an inter-sample VCF with kanpig, and keeps only records that
# occur in at least one sample after re-genotyping.
#
workflow KanpigGenotypeAndFilter {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        Array[String] sample_id
        Array[Boolean] is_male
        Array[File] alignments_bam
        Array[File] alignments_bai
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        Int min_sv_length = 10
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
                sample_id = sample_id[i],
                is_male = is_male[i],
                input_vcf_gz = RemoveSamples.cleaned_vcf_gz,
                input_tbi = RemoveSamples.cleaned_tbi,
                alignments_bam = alignments_bam[i],
                alignments_bai = alignments_bai[i],
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                ploidy_bed_female = ploidy_bed_female,
                ploidy_bed_male = ploidy_bed_male,
                min_sv_length = min_sv_length
        }
    }
    call Merge {
        input:
            intersample_vcf_gz = RemoveSamples.cleaned_vcf_gz,
            format_column = Kanpig.regenotyped_format[0],
            gt_columns = Kanpig.regenotyped_gt
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
        bgzip -@ ${N_THREADS} --compress-level 1 cleaned.vcf
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
        disks: "local-disk 256 SSD"
        preemptible: 0
    }
}


# Parameters are tuned for a multi-sample VCF.
#
task Kanpig {
    input {
        String sample_id
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
        Int min_sv_length
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    String output_prefix = "kanpig_regenotyped"
    String kanpig_params_multisample =  "--sizemin ~{min_sv_length} --sizemax 10000 --neighdist 500 --gpenalty 0.04 --hapsim 0.97"
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
        bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp1.vcf.gz > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        rm -f tmp1.vcf.gz*
        
        # Outputting just the FORMAT and GT columns
        bcftools view --threads ${N_THREADS} --no-header tmp2.vcf.gz | cut -f 9 >> format.txt
        echo ~{sample_id} > gt.txt
        bcftools view --threads ${N_THREADS} --no-header tmp2.vcf.gz | cut -f 10 >> gt.txt
    >>>

    output {
        File regenotyped_gt = work_dir + "/gt.txt"
        File regenotyped_format = work_dir + "/format.txt"
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
        File intersample_vcf_gz
        File format_column
        Array[File] gt_columns
        Int n_cpus = 8
    }
    parameter_meta {
        intersample_vcf_gz: "The input to kanpig."
        format_column: "A FORMAT column from the output of kanpig (header NOT included)."
        gt_columns: "Just the GT column from the output of kanpig. The first line must be the sample ID."
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
                head -n 1 ${ADDRESS} > sample_${THREAD_ID}.txt
                if [ $i = "1" ]; then
                    mv sample_${THREAD_ID}.txt ${FIELDS_FILE}
                else
                    paste ${FIELDS_FILE} sample_${THREAD_ID}.txt > ${FIELDS_FILE}.prime
                    mv ${FIELDS_FILE}.prime ${FIELDS_FILE}
                fi
                echo "Current fields of thread ${THREAD_ID}:"; cat ${FIELDS_FILE}
                # Adding the new column to the body
                tail -n +2 ${ADDRESS} > ${TMP_PREFIX}.txt
                if [ $i = "1" ]; then
                    mv ${TMP_PREFIX}.txt ${OUTPUT_FILE}
                else
                    ${TIME_COMMAND} paste ${OUTPUT_FILE} ${TMP_PREFIX}.txt > ${OUTPUT_FILE}.prime
                    mv ${OUTPUT_FILE}.prime ${OUTPUT_FILE}
                fi
            done < list_${THREAD_ID}
        }

        
        # Main program
        INPUT_FILES=~{sep=',' gt_columns}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            echo ${INPUT_FILE} >> list.txt
        done
        
        # Initializing the inter-sample VCF
        bcftools view --header-only ~{intersample_vcf_gz} > tmp.txt
        N_ROWS=$(wc -l < tmp.txt)
        head -n $(( ${N_ROWS} - 1 )) tmp.txt > header.txt
        echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Kanplug genotype\">" >> header.txt
        echo -e "##FORMAT=<ID=FT,Number=1,Type=Integer,Description=\"Kanpig filter\">" >> header.txt
        echo -e "##FORMAT=<ID=SQ,Number=1,Type=Integer,Description=\"Phred scaled quality of sample being non-ref at this variant\">" >> header.txt
        echo -e "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred scaled quality of genotype\">" >> header.txt
        echo -e "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Local phase group of entries\">" >> header.txt
        echo -e "##FORMAT=<ID=NE,Number=1,Type=Integer,Description=\"Neighborhood id of variants evaluated together for short-range phasing\">" >> header.txt
        echo -e "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Coverage over region\">" >> header.txt
        echo -e "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Coverage for reference and alternate alleles\">" >> header.txt
        echo -e "##FORMAT=<ID=KS,Number=.,Type=Integer,Description=\"Kanpig score\">" >> header.txt
        tail -n 1 tmp.txt | cut -f 1,2,3,4,5,6,7,8,9 > fields.txt
        bcftools view --no-header ~{intersample_vcf_gz} | cut -f 1,2,3,4,5,6,7,8 > calls.txt
        ${TIME_COMMAND} paste calls.txt ~{format_column} > calls2.txt
        rm -f calls.txt
        mv calls2.txt calls.txt
        
        # Appending the sample columns
        N_ROWS=$(wc -l < list.txt)
        N_ROWS=$(( (${N_ROWS}+${N_THREADS}-1)/${N_THREADS} ))  # Ceiling
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
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} --compress-level 1 merged.vcf
        tabix -f merged.vcf.gz
        
        # Keeping all and only the records that occur at least once in some
        # sample
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
        memory: "32GB"
        disks: "local-disk 100 SSD"
        preemptible: 0
    }
}
