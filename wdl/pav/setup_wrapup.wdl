version 1.0

task tar_asm {
  input {
    File ref
    File hapOne
    File hapTwo
    String sample
    String threads
    String mem_gb
  }
  command <<<
    set -eux
    if cmp --silent ~{hapOne} ~{hapTwo} ; then
      echo "Haplotype FASTAs are identical" && exit 1
    fi
    mkdir -p asm/~{sample}
    cp ~{ref} asm/ref.fa
    samtools faidx asm/ref.fa
    cp ~{hapOne} asm/~{sample}/h1.fa.gz
    cp ~{hapTwo} asm/~{sample}/h2.fa.gz
    tar zcvf asm.tgz asm/
  >>>
  output {
    File asm_tar = "asm.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 30 + " SSD"
      bootDiskSizeGb: 10
      preemptible:    1
      maxRetries:     0
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1"
  }
}

task call_final_bed {
  input {
    File pav_conf
    File pav_asm
    File invBed
    File insBed
    File delBed
    File snvBed
    String threads
    String mem_gb
    String sample
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{invBed}
    tar zxvf ~{snvBed}
    tar zxvf ~{insBed}
    tar zxvf ~{delBed}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} results/~{sample}/bed/snv_snv.bed.gz results/~{sample}/bed/indel_ins.bed.gz results/~{sample}/bed/indel_del.bed.gz results/~{sample}/bed/sv_ins.bed.gz results/~{sample}/bed/sv_del.bed.gz results/~{sample}/bed/sv_inv.bed.gz results/~{sample}/bed/fa/indel_ins.fa.gz results/~{sample}/bed/fa/indel_del.fa.gz results/~{sample}/bed/fa/sv_ins.fa.gz results/~{sample}/bed/fa/sv_del.fa.gz results/~{sample}/bed/fa/sv_inv.fa.gz
    tar zcvf final_bed.tgz results/~{sample}/bed/snv_snv.bed.gz results/~{sample}/bed/indel_ins.bed.gz results/~{sample}/bed/indel_del.bed.gz results/~{sample}/bed/sv_ins.bed.gz results/~{sample}/bed/sv_del.bed.gz results/~{sample}/bed/sv_inv.bed.gz results/~{sample}/bed/fa/indel_ins.fa.gz results/~{sample}/bed/fa/indel_del.fa.gz results/~{sample}/bed/fa/sv_ins.fa.gz results/~{sample}/bed/fa/sv_del.fa.gz results/~{sample}/bed/fa/sv_inv.fa.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File bed = "final_bed.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 30 + " SSD"
      bootDiskSizeGb: 10
      preemptible:    1
      maxRetries:     0
      docker:         "fcunial/hapestry_experiments"
  }
}

task data_ref_contig_table{
  input {
    File pav_conf
    File pav_asm
    String threads
    String mem_gb
    String sample
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s snakemake -s pav/Snakefile --cores ~{threads} data/ref/contig_info.tsv.gz
    tar zcvf contig_info.tgz data/ref/contig_info.tsv.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File contigInfo = 'contig_info.tgz'
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 30 + " SSD"
      bootDiskSizeGb: 10
      preemptible:    1
      maxRetries:     0
      docker:         "fcunial/hapestry_experiments"
  }
}

task write_vcf {
  input {
    File pav_conf
    File pav_asm
    File contigInfo
    File finalBedOut
    String threads
    String mem_gb
    String sample
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{contigInfo}
    tar zxvf ~{finalBedOut}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} pav_~{sample}.vcf.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File vcf = "pav_~{sample}.vcf.gz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 30 + " SSD"
      bootDiskSizeGb: 10
      preemptible:    1
      maxRetries:     0
      docker:         "fcunial/hapestry_experiments"
  }
}

task FilterChromosomes {
  meta {
        description: "For down the contigs listed in the FAI to just those targetted."
  }
  input {
    File refFai
    File targetChromsList
  }

  command <<<
    set -eux
    grep -wf ~{targetChromsList} ~{refFai} > "filtered.fai"
  >>>

  output {
    Array[Array[String]] targetChromes = read_tsv("filtered.fai")
  }

  runtime {
    disks: "local-disk 10 HDD"
    docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    preemptible:    1
    maxRetries:     0
  }
}

task CollapseStrings {
    meta {
        description: "For collapsing an array of strings into a long single-space-delimited string"
    }
    input {
        Array[String] whatever
    }

    command <<<
        echo ~{sep=' ' whatever}
    >>>

    output {
        String collapsed = read_string(stdout())
    }

    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
        preemptible:    1
    maxRetries:     0
    }
}

task IndexVcf {

    meta {
        description: "Indexing vcf.gz. Note: do NOT use remote index as that's buggy."
    }

    input {
        File vcf
    }

    String prefix = basename(vcf, ".vcf.gz")

    command <<<
      set -eux
      cp ~{vcf} ~{prefix}.vcf.gz && \
        tabix -p vcf ~{prefix}.vcf.gz && \
        find ./ -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g'
    >>>

    output {
        File tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    runtime {
        cpu:          8
        memory:       "32 GiB"
        disks:        "local-disk 375 LOCAL"
        preemptible:  1
        maxRetries:   0
        docker:       "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
    }
}

task FinalizeToFile {
    input {
        File file
        String outdir
        String? name
    }

    parameter_meta {
        file: {
            localization_optional: true
        }
        outdir: "directory to which files should be uploaded"
        name:   "name to set for uploaded file"
    }

    String gcs_output_dir = sub(outdir, "/+$", "")
    String gcs_output_file = gcs_output_dir + "/" + select_first([name, basename(file)])

    command <<<
        set -euxo pipefail

        gsutil -m cp "~{file}" "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    runtime {
        cpu:    1
        memory: "4 GiB"
        disks:  "local-disk " + 30 + " SSD"
        preemptible: 2
        maxRetries:  0
        docker: "fcunial/hapestry_experiments"
    }
}
