# Code
All the WDLs in this workspace are available [on GitHub](https://github.com/fabio-cunial/hapestry_experiments).


# Tools
See [the dockerfile](https://github.com/fabio-cunial/hapestry_experiments/blob/main/docker/Dockerfile) for more details.

* bcftools 1.21
* cutesv 2.1.1
* deepvariant 1.6.0
* dipcall 0.3
* hifiasm 0.24.0
* jasmine 1.1.5
* panpop NC2024
* [pav 1.1.2](https://dockstore.org/workflows/github.com/broadinstitute/pav-wdl/pav:sh_more_resources_pete_t2thotfix?tab=info), the version used in AoU (v2.4.6 does not work with the AoU WDLs)
* pbsv 2.10.0
* sniffles 2.5.3
* trgt 1.5.1
* truvari 5.0.0



# Bucket address
Bucket address: `gs://fc-a90ab401-9c4b-43d1-b891-f0410c667ff2`

# Reference genome
* GRCh38 without ALT contigs, as suggested [in this post](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) (see the corresponding workspace variable).
* As input to SV callers, we used the tandem repeat track [distributed with pbsv](https://github.com/PacificBiosciences/pbsv/tree/master/annotations), since we observed it captures most of the TRs in the reference.
* For TRGT we used the [tandem repeat catalog + variation clusters v1.0](https://github.com/broadinstitute/tandem-repeat-catalog/releases/download/v1.0/variation_clusters_and_isolated_TRs_v1.hg38.TRGT.bed.gz), since it was the most extensive at the time of the study.


# HPRC Y1
The data is organized in the following tables.

## Table `hprc_y1_sample`
Contains basic information about the samples, including the ground truth used in all experiments. Columns are organized in the following sections:

### Section 1 - General information on samples and remote addresses of source files
This information is copied from the [official HPRC Y1 AnVIL workspace](https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_HPRC/data) (tables `sample` and `assembly_sample`). 

* Column `01_souce_ccs` is a copy of AnVIL's column `hifi`of table `sample`, but for every cell and for every distinct filename (disregarding extension) we kept just one file with that name (preferring FASTQs to BAMs).

### Section 2 - Truth data: assemblies and dipcall

* The remote assembly addresses are copied from AnVIL's columns `mat_fasta`, `pat_fasta` of table `assembly_sample`.
* We created a fork of dipcall ([dipcall_asm10.wdl](https://github.com/rlorigro/dipcall_asm10)) that runs minimap2 with the `-x asm10` assembly-to-reference preset, instead of the default `-x asm5`. This is because we observed that the default setting misses several SVs that are supported by read alignments.
* The raw dipcall output (column `dipcall_vcf`) contains all variants, including SNPs and short INDELs, and some of its records are multiallelic. We transform multiallelic records into biallelic using `bcftools norm`, and we keep only variants of length at least `L` for `L` in {10,20,30,40,50} with a custom script (column `dipcall_normed_filtered_L`, workflow [FilterDipcall.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/FilterDipcall.wdl)). 
	* Column `dipcall_normed_filtered` is the result of applying FilterDipcall.wdl with `L`=0.
* For use in IGV, we separately aligned every assembly haplotype to the reference, using [wdl-minimap2](https://github.com/jmonlong/wdl-minimap2/tree/main) with parameters `kmerSize=19, secondary=false, softClipOnly=false, eqx=true, indexSplitBp=8g, preset=asm10, otherArgs=-L, batchSize=2g` (see e.g. [this run](https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/callset_integration/job_history/3c8ce068-12e3-4fb9-8b5f-8dcd02286461) in the [callset_integration workspace](https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/callset_integration); these alignments were created in that workspace and are just pointed from the tables in this workspace). Apart from asm10, such parameters were chosen to speed up computation.

## Table `hprc_y1_Cx`

Contains subsampled reads at coverage `C`, and the corresponding alignments and calls. Cotains also intra-sample merges of multiple callsets. Columns are organized in the following sections:

### Section 1 - Subsampled reads and corresponding alignments
* Reads are subsampled with [Subsample.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/Subsample.wdl). Briefly, given the list of remote source files of a sample, the list is randomly permuted. For each file in the permuted order, the file is appended to a growing FASTQ until the desired number of bases (measured with `seqkit stats`) is reached or exceeded. Finally, for each coverage, the FASTQ is randomly sampled using `seqkit sample`. If the concatenation of all remote files does not achieve a given coverage, that coverage's FASTQ is set to the concatenation of all files and a warning is issued. The subsampled FASTQ at the `i`-th coverage is not necessarily a superset of the subsampled FASTQ at the `i-1`-th coverage.
* CCS reads are mapped to the reference with pbmm2 using default parameters (workflow [MapCCS.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/MapCCS.wdl)).

### Section 2 - SV callers
We consider a production run to consist only of callers pbsv, sniffles, cutesv. Except otherwise noted, every SV merging tool is fed all and only such callers in input. Pav is used in a separate, ad hoc set of experiments that study the effect of FPs.

* All alignment-based callers are run with `min_sv_length`=10 and all other parameters left to default values. 
* Columns `02_C_resolved_vcf` for every caller `C`contain a normalized version of its raw output, produced by [Resolve.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/Resolve.wdl).
* Pav and hifiasm were run with default parameters (workflows [Hifiasm.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/Hifiasm.wdl), and [pav.wdl](https://dockstore.org/workflows/github.com/broadinstitute/pav-wdl/pav:sh_more_resources_pete_t2thotfix?tab=info) from AoU production), and then only calls >=10bp were kept (workflow [ResolvePav.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/ResolvePav.wdl)).

### Section 3 - Intra-sample merge of SV callers

* Column `03_bcftools_merge_vcf` contains a trivial merge of all the `02_C_resolved_vcf` columns, for all `C` ([BcftoolsMergeIntrasample.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/BcftoolsMergeIntrasample.wdl)). This VCF contains exactly one sample column.
* Column `03_truvari_vcf` contains the truvari collapse of column `03_bcftools_merge_vcf`, with all parameters set to default except the following: `--sizemin 0 --sizemax 1000000 --pctseq 0.90 --pctsize 0.90`. See [Truvari.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/Truvari.wdl) for an explanation of these parameters.
* Column `03_jasmine_vcf` contains the jasmine collapse of column `03_bcftools_merge_vcf`, with all parameters set to default except the following: `--allow_intrasample min_seq_id=0.9`. See [Jasmine.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/Jasmine.wdl) for an explanation of these parameters and of some cleaning steps. 
	* **Remark: all the DEL and INV calls in these files are symbolic: this was needed to make jasmine finish in a reasonable amount of time.**
	* **Remark: jasmine does not allow to set min SV length on the command line, so its behavior with calls <50bp is unclear.**
* Column `03_panpop_vcf` contains the panpop collapse of column `03_bcftools_merge_vcf` (workflow [Panpop.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/Panpop.wdl)).
	* **Remark: only calls <=50kb are kept, to make panpop finish in a reasonable amount of time.** 
	* **Remark: panpop has no parameters that can be set from the command line.**

### Section 4 - Non-SV calls
* Column `04_dv_g_vcf` contains deepvariant GVCFs computed with workflow [CallVariantsReadBased.wdl](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/CallVariantsReadBased:sh_update_wgs_callers?tab=versions) (branch `sh_update_wgs_callers`).
* Column `04_trgt_vcf` contains the raw TRGT VCFs computed with workflow [Trgt.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/Trgt.wdl). Column `04_trgt_resolved_vcf` contains a version of `04_trgt_vcf` where records with missing or ref GT are discarded, and multiallelic records are split into multiple records. 

### Section 5 - Intra-sample merge of SV callers with PAV

Like Section 3, but with PAV added.

### Section 6 - Intra-sample merge of pbsv, sniffles, PAV

Like Section 3, but with pbsv, sniffles, PAV.

* Column `06_jasmine_vcf` contains the jasmine collapse of column `06_bcftools_merge_vcf`, with all parameters set to default except the following: `--allow_intrasample min_seq_id=0.9`. See [Jasmine.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/Jasmine.wdl) for an explanation of these parameters and of some cleaning steps. 
	* **Remark: only calls <=50kb are kept, to make jasmine finish in a reasonable amount of time.**
	* **DEL and INV calls in these files are NOT forced to be symbolic.**
	* **Remark: jasmine does not allow to set min SV length on the command line, so its behavior with calls <50bp is unclear.**


## Table `hprc_y1_cohort`
Cohort-level VCFs that aggregate all samples at each coverage.

### Section 1 - Naive merge of all the dipcall calls
* Column `01_dipcall_L_vcf`  for each `L`contains the bcftools merge of all the HPRC Y1 dipcall VCFs in column `hprc_y1_sample.02_dipcall_normed_filtered_L_vcf` (workflow [BcftoolsMergeDipcall.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/BcftoolsMergeDipcall.wdl)). 
* Column `01_dipcall_vcf` contains the bcftools merge  of all the VCFs in column `hprc_y1_sample.02_dipcall_normed_filtered_vcf`.

### Section 2 - Joint SV callers
* Column `02_sniffles_joint_vcf` contains the joint VCF created by sniffles when run on the SNF files of all samples at a specific coverage with default parameters (workflow [SnifflesIntersample.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/SnifflesIntersample.wdl)). The raw VCF is then cleaned with [Resolve.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/Resolve.wdl).
* Column `02_pbsv_joint_vcf` contains the joint VCF created by pbsv when run on the SVSIG files of all samples at a specific coverage with default parameters (workflow [PbsvIntersample.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/PbsvIntersample.wdl)). The raw VCF is then cleaned with [Resolve.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/Resolve.wdl).

### Section 3 - SV merging tools
* Column `03_bcftools_merge_vcf` contains the bcftools merge of all the VCFs in column `hprc_y1_Cx.03_bcftools_merge_vcf` for a specific coverage `C` (workflow [BcftoolsMergeIntersample.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/BcftoolsMergeIntersample.wdl)). These calls are the naive union of all callers and all samples VCFs at coverage `C`, obained by just doing an inter-sample bcftools merge of all the intra-sample bcftools merge VCFs.
* Truvari:
	* Column `03_truvari_mode1_vcf` represents how truvari collapse would be run in a production pipeline similar to AoU. It corresponds to the pipeline: Intra-sample bcftools merge -> `truvari collapse --sizemin 0 --sizemax 1000000 --pctseq 0.90 --pctsize 0.90` -> Inter-sample bcftools merge -> `truvari collapse --sizemin 0 --sizemax 1000000 --keep common --gt all`. See [Truvari.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/Truvari.wdl) for an explanation of these parameters.
	* Column `03_truvari_mode2_vcf` tests truvari collapse on the set of all raw calls from all callers stored in column `03_bcftools_merge_vcf`. It corresponds to the pipeline: Intra-sample bcftools merge -> Inter-sample bcftools merge -> `truvari collapse --sizemin 0 --sizemax 1000000 --keep common`. See [Truvari.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/Truvari.wdl) for an explanation of these parameters.
* Jasmine:
	* Column `03_jasmine_mode1_vcf` represents how jasmine would be run in a production pipeline similar to AoU. It corresponds to the pipeline: Intra-sample bcftools merge -> `jasmine --allow_intrasample min_seq_id=0.9` -> `jasmine min_seq_id=0.7`. See [JasmineIntersample1.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/JasmineIntersample1.wdl) for an explanation of these parameters. 
		* **Remark: this inter-sample VCF has only one sample column, since using --output_genotypes makes jasmine crash.**
		* **Remark: jasmine does not allow to set min SV length on the command line, so its behavior with calls <50bp is unclear.**
	* Column  `03_jasmine_mode2_vcf` tests jasmine on the set of all raw calls from all callers stored in column `03_bcftools_merge_vcf`. It corresponds to the pipeline: Intra-sample bcftools merge -> Inter-sample bcftools merge -> `jasmine --allow_intrasample min_seq_id=0.7`. See [Jasmine.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/Jasmine.wdl) for an explanation of these parameters. 
		* **Remark: all the DEL and INV calls in these files are symbolic, and only calls <=50kb are kept: this was needed to make jasmine finish in a reasonable amount of time.**
		* **Remark: jasmine does not allow to set min SV length on the command line, so its behavior with calls <50bp is unclear.**
* Panpop:
	* Column `03_panpop_mode1_vcf` represents how panpop would be run in a production pipeline similar to AoU. It corresponds to the pipeline: Intra-sample bcftools merge -> panpop -> Inter-sample bcftools merge -> panpop. Note that panpop has no parameters that can be set from the command line.
	* Column `03_panpop_mode2_vcf` tests panpop on the set of all raw calls from all callers stored in column `03_bcftools_merge_vcf`. It corresponds to the pipeline: Intra-sample bcftools merge -> Inter-sample bcftools merge -> panpop. 
		* **Remark: only calls <=50kbp are  kept, to make panpop finish in a reasonable amount of time.**
		* **Remark: panpop has no parameters that can be set from the command line.**

### Section 4 - Non-SV calls
* Column `04_dv_vcf` contains joint-called SNPs and small indels created by running workflow [LRJointCallGVCFs.wdl](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/LRJointCallGVCFs:sh_a_bit_refactor?tab=info) on the GVCFs in column `hprc_y1_Cx.04_dv_g_vcf`. 
	* **Remark: this file contains multiallelic records.**
* Column `04_trgt_merge_vcf` contains the result of running `trgt merge` on the raw TRGT VCFs in column `hprc_y1_Cx.04_trgt_vcf`. Only calls in the output that are non-ref in some sample are kept, and multiallelic records are converted to biallelic (workflow [TrgtMerge.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/TrgtMerge.wdl)).
* Column `04_bcftols_merge_trgt_vcf` contains the bcftools merge of the resolved TRGT VCFs in column `hprc_y1_Cx.04_trgt_resolved_vcf` (workflow [BcftoolsMergeIntersample.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/BcftoolsMergeIntersample.wdl)).

### Section 5 - SV and non-SV calls
* Column `05_bcftools_merge_sv_dv_vcf` contains the bcftools concat of a non-multiallelic version of `04_dv_vcf` and `03_bcftools_merge_vcf` (workflow [BcftoolsMergeSvSnp.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/BcftoolsMergeSvSnp.wdl)).
* Column `05_bcftools_merge_sv_trgt1_vcf` contains the bcftools concat of `04_trgt_merge_vcf` and `03_bcftools_merge_vcf`  (workflow [BcftoolsMergeSvTrgt.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/BcftoolsMergeSvTrgt.wdl)).
* Column `05_bcftools_merge_sv_trgt2_vcf` contains the bcftools concat of `04_bcftols_merge_trgt_vcf` and `03_bcftools_merge_vcf`  (workflow [BcftoolsMergeSvTrgt.wdl](https://github.com/fabio-cunial/hapestry_experiments/blob/main/wdl/BcftoolsMergeSvTrgt.wdl)).

### Section 6 - SV merging tools with PAV

Same as Section 3, but with PAV added.

### Section 7 - Pbsv, sniffles, pav.

Same as Section 3, but with pbsv, sniffles, PAV.