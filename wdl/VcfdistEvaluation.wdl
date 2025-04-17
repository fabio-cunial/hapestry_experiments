version 1.0
import "clean_vcf.wdl" as clean_vcf


# task to average any Array[Float]
task ComputeAverage {
    input {
        Array[Float] x
    }

    command <<<
python <<CODE
arr = ["~{sep="\", \"" x}"]
floats = [float(x) for x in arr]
mean = sum(floats) / len(floats)
with open('mean.txt', 'w') as file:
    file.write(f'{mean}\n')
CODE
    >>>

    runtime {
        cpu: 1
        memory: "8 GiB"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "python:3.9.20-slim-bullseye"
    }

    output {
        Float y = read_float("mean.txt")
    }
}


workflow VcfdistEvaluation {
    input {
        Array[String] samples
        File truth_vcf
        File? truth_vcf_tbi
        File eval_vcf
        File? eval_vcf_tbi
        String region
        File reference_fasta
        File reference_fasta_fai

        File vcfdist_bed_file
        String? vcfdist_extra_args
    }

    scatter (sample in samples) {

        call SubsetSampleFromVcf as SubsetSampleFromVcfEval { input:
            vcf = eval_vcf,
            vcf_tbi = eval_vcf_tbi,
            sample = sample,
            region = region,
            reference_fasta_fai = reference_fasta_fai
        }

        # Clean the single-sample VCF (don't do this in first step because VCF too big)
        call clean_vcf.CleanVcfAlleles as clean { input:
            vcf_gz = SubsetSampleFromVcfEval.single_sample_vcf,
            vcf_tbi = SubsetSampleFromVcfEval.single_sample_vcf_tbi,
            ref_fasta = reference_fasta
        }

        call SubsetSampleFromVcf as SubsetSampleFromVcfTruth { input:
            vcf = truth_vcf,
            vcf_tbi = truth_vcf_tbi,
            sample = sample,
            region = region,
            reference_fasta_fai = reference_fasta_fai
        }

        call Vcfdist as vcfdist { input:
            sample = sample,
            eval_vcf = clean.output_vcf_gz,
            truth_vcf = SubsetSampleFromVcfTruth.single_sample_vcf,
            bed_file = vcfdist_bed_file,
            reference_fasta = reference_fasta,
            extra_args = vcfdist_extra_args
        }
    }

    # Compute averages for all categories
    call ComputeAverage as avg_sv_p { input: x = vcfdist.SV_PREC }
    call ComputeAverage as avg_sv_r { input: x = vcfdist.SV_RECALL }
    call ComputeAverage as avg_sv_f { input: x = vcfdist.SV_F1_SCORE }

    call ComputeAverage as avg_snp_p { input: x = vcfdist.SNP_PREC }
    call ComputeAverage as avg_snp_r { input: x = vcfdist.SNP_RECALL }
    call ComputeAverage as avg_snp_f { input: x = vcfdist.SNP_F1_SCORE }

    call ComputeAverage as avg_indel_p { input: x = vcfdist.INDEL_PREC }
    call ComputeAverage as avg_indel_r { input: x = vcfdist.INDEL_RECALL }
    call ComputeAverage as avg_indel_f { input: x = vcfdist.INDEL_F1_SCORE }

    call ComputeAverage as avg_all_p { input: x = vcfdist.ALL_PREC }
    call ComputeAverage as avg_all_r { input: x = vcfdist.ALL_RECALL }
    call ComputeAverage as avg_all_f { input: x = vcfdist.ALL_F1_SCORE }

    output {
        # per-sample
        Array[File] vcfdist_summary_vcf = vcfdist.summary_vcf
        Array[File] vcfdist_precision_recall_summary_tsv = vcfdist.precision_recall_summary_tsv
        Array[File] vcfdist_superclusters_tsv = vcfdist.superclusters_tsv

        # SV metrics
        Float SV_PREC_avg = avg_sv_p.y
        Float SV_RECALL_avg = avg_sv_r.y
        Float SV_F1_SCORE_avg = avg_sv_f.y

        # SNP metrics
        Float SNP_PREC_avg = avg_snp_p.y
        Float SNP_RECALL_avg = avg_snp_r.y
        Float SNP_F1_SCORE_avg = avg_snp_f.y

        # INDEL metrics
        Float INDEL_PREC_avg = avg_indel_p.y
        Float INDEL_RECALL_avg = avg_indel_r.y
        Float INDEL_F1_SCORE_avg = avg_indel_f.y

        # ALL metrics
        Float ALL_PREC_avg = avg_all_p.y
        Float ALL_RECALL_avg = avg_all_r.y
        Float ALL_F1_SCORE_avg = avg_all_f.y
    }
}


task SubsetSampleFromVcf {
    input {
        File vcf
        File? vcf_tbi
        String sample
        String region
        File reference_fasta_fai
    }

    Int disk_size = 1 + ceil(2 * (size(vcf, "GiB")))

    command <<<
        set -euxo pipefail

        if ~{!defined(vcf_tbi)}; then
            bcftools index --threads 4 ~{vcf}
        fi

        bcftools view ~{vcf} \
            -r ~{region} \
            --threads 4 \
            -Oz -o region.vcf.gz

        bcftools index -t --threads 4 region.vcf.gz

        bcftools view region.vcf.gz \
            --threads 4 \
            -s ~{sample} \
            -r ~{region} \
            -Oz | bcftools +fill-tags -- -t AF | bcftools view -i 'AF>0' -Oz -o ~{sample}.subset.g.vcf.gz

        bcftools reheader ~{sample}.subset.g.vcf.gz \
            --fai ~{reference_fasta_fai} \
            -o ~{sample}.subset.reheadered.g.vcf.gz

        bcftools index --threads 4 -t ~{sample}.subset.reheadered.g.vcf.gz
    >>>

    output {
        File single_sample_vcf = "~{sample}.subset.reheadered.g.vcf.gz"
        File single_sample_vcf_tbi = "~{sample}.subset.reheadered.g.vcf.gz.tbi"
    }

    runtime {
        cpu: 16
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " SSD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "staphb/bcftools:1.20"
    }
}


task Vcfdist {
    input {
        String sample
        File eval_vcf
        File truth_vcf
        File bed_file
        File reference_fasta
        String? extra_args
        Int verbosity = 1

        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 16
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        vcfdist \
            ~{eval_vcf} \
            ~{truth_vcf} \
            ~{reference_fasta} \
            -b ~{bed_file} \
            -v ~{verbosity} \
            ~{extra_args}

        extract_metrics() {
            local VAR_TYPE="$1"
            # Extract the row where VAR_TYPE matches the provided variant and THRESHOLD = 'BEST'
            row=$(awk -F'\t' -v var="$VAR_TYPE" '
                NR==1 {for (i=1; i<=NF; i++) header[i]=$i}
                $1==var && $2=="BEST" {
                    for (i=1; i<=NF; i++) print header[i], $i
                }' precision-recall-summary.tsv)

            # Extract PREC, RECALL, and F1_SCORE values
            PREC=$(echo "$row" | awk '/PREC/ {print $2}')
            RECALL=$(echo "$row" | awk '/RECALL/ {print $2}')
            F1_SCORE=$(echo "$row" | awk '/F1_SCORE/ {print $2}')

            # Print the extracted values
            echo "${VAR_TYPE}_PREC: $PREC"
            echo "${VAR_TYPE}_RECALL: $RECALL"
            echo "${VAR_TYPE}_F1_SCORE: $F1_SCORE"

            # Write the values to files
            echo "$PREC" > "${VAR_TYPE}_PREC"
            echo "$RECALL" > "${VAR_TYPE}_RECALL"
            echo "$F1_SCORE" > "${VAR_TYPE}_F1_SCORE"
        }

        extract_metrics "SV"
        extract_metrics "INDEL"
        extract_metrics "SNP"
        extract_metrics "ALL"

        for tsv in $(ls *.tsv); do mv $tsv ~{sample}.$tsv; done
        mv summary.vcf ~{sample}.summary.vcf
    >>>

    output {
        File summary_vcf = "~{sample}.summary.vcf"
        File precision_recall_summary_tsv = "~{sample}.precision-recall-summary.tsv"
        File precision_recall_tsv = "~{sample}.precision-recall.tsv"
        File query_tsv = "~{sample}.query.tsv"
        File truth_tsv = "~{sample}.truth.tsv"
        File phasing_summary_tsv = "~{sample}.phasing-summary.tsv"
        File switchflips_tsv = "~{sample}.switchflips.tsv"
        File superclusters_tsv = "~{sample}.superclusters.tsv"
        File phase_blocks_tsv = "~{sample}.phase-blocks.tsv"

        # SV metrics
        Float SV_PREC = read_float("SV_PREC")
        Float SV_RECALL = read_float("SV_RECALL")
        Float SV_F1_SCORE = read_float("SV_F1_SCORE")

        # SNP metrics
        Float SNP_PREC = read_float("SNP_PREC")
        Float SNP_RECALL = read_float("SNP_RECALL")
        Float SNP_F1_SCORE = read_float("SNP_F1_SCORE")

        # INDEL metrics
        Float INDEL_PREC = read_float("INDEL_PREC")
        Float INDEL_RECALL = read_float("INDEL_RECALL")
        Float INDEL_F1_SCORE = read_float("INDEL_F1_SCORE")

        # ALL metrics
        Float ALL_PREC = read_float("ALL_PREC")
        Float ALL_RECALL = read_float("ALL_RECALL")
        Float ALL_F1_SCORE = read_float("ALL_F1_SCORE")
    }

    runtime {
        docker: "timd1/vcfdist:v2.5.3"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}
