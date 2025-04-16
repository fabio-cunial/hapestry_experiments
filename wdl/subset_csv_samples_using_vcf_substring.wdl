version 1.0


task SubsetSamplesFromCsv {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File csv
    }

    Int disk_size = 1 + ceil(2 * (size(vcf_gz, "GiB")))

  command <<<
    # Extract sample names from VCF
    bcftools query -l ~{vcf_gz} > vcf_samples.txt

    python3 - <<EOF
import csv

# Read VCF sample names
with open("vcf_samples.txt") as f:
    vcf_samples = f.read().splitlines()

# Open input and output CSVs
with open("~{csv}", newline='') as infile, open("output.csv", "w", newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)
    for row in reader:
        if any(sample in row[0] for sample in vcf_samples):
            writer.writerow(row)
EOF
  >>>
    output {
        File subsampled_csv = "output.csv"
    }

    runtime {
        cpu: 4
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " SSD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "fcunial/hapestry_experiments"
    }
}


workflow SubsetSamplesFromCsvUsingVcfSubstring {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File csv
    }

    call SubsetSamplesFromCsv as subset { input:
        vcf_gz = vcf_gz,
        vcf_gz_tbi = vcf_gz_tbi,
        csv = csv
    }

    output {
        File subsampled_csv = subset.subsampled_csv
    }
}
