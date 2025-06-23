To get the tandem repeat track windows that intersect with at least two TRGT intervals:

bedtools intersect -c -a human_GRCh38_no_alt_analysis_set.trf.bed -b variation_clusters_and_isolated_TRs_v1.hg38.TRGT.bed.gz | awk '{ if ($4>1) print $0; }' > human_GRCh38_no_alt_analysis_set.trf.two_target_intervals.bed

Then, to select only windows with at least a call in the 8x bcftools merge VCF at 50bp:

bedtools intersect -u -a human_GRCh38_no_alt_analysis_set.trf.two_target_intervals.bed -b merged.vcf.gz > 
