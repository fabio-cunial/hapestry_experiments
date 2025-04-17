version 1.0


task SubsetNSamples {
    input {
        File vcf_gz
        File vcf_gz_tbi
        Int n
        Boolean use_last
    }

    Int disk_size = 1 + ceil(2 * (size(vcf_gz, "GiB")))

    command <<<
        set -euxo pipefail

        COMMAND="head"
        if ~{use_last}; then
            COMMAND="tail"
        fi

        # Extract the list of sample names
        SAMPLES=$(bcftools query -l ~{vcf_gz} | $COMMAND -n ~{n} | paste -sd, -)

        # subset to first/last n samples and filter on AF > 0
        bcftools view -Oz -s "$SAMPLES" ~{vcf_gz} | bcftools +fill-tags -- -t AF | bcftools view -i 'AF>0' -Oz -o n_samples.vcf.gz

        bcftools index -t --threads 4 n_samples.vcf.gz
    >>>

    output {
        File n_samples_vcf_gz = "n_samples.vcf.gz"
        File n_samples_vcf_gz_tbi = "n_samples.vcf.gz.tbi"
    }

    runtime {
        cpu: 4
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " SSD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "staphb/bcftools:1.20"
    }
}

workflow SubsetNSamplesFromVcf {
    input {
        Int n
        File vcf_gz
        File vcf_gz_tbi
        Boolean use_last = false
    }

    call SubsetNSamples as subset {
        input:
            vcf_gz = vcf_gz,
            vcf_gz_tbi = vcf_gz_tbi,
            n = n,
            use_last = use_last
    }

    output {
        File n_samples_vcf_gz = subset.n_samples_vcf_gz
        File n_samples_vcf_gz_tbi = subset.n_samples_vcf_gz_tbi
    }

    parameter_meta {
        n: "Number of samples to include in the subset"
        vcf_gz: "Input VCF file in bgzipped format"
        vcf_gz_tbi: "Index file (.tbi) for the input VCF"
        use_last: "If true, take the last n samples instead of the first n"
    }
}

