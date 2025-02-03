version 1.0


# Ensures that the VCF from an SV caller is in a standard format.
#
workflow Resolve {
    input {
        String sample_id
        File vcf_gz
        File tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    call ResolveImpl {
        input:
            sample_id = sample_id,
            vcf_gz = vcf_gz,
            tbi = tbi,
            reference_fa = reference_fa
    }
    
    output {
    	File resolved_vcf_gz = ResolveImpl.resolved_vcf_gz
    	File resolved_tbi = ResolveImpl.resolved_tbi
    }
}


#
task ResolveImpl {
    input {
        String sample_id
        File vcf_gz
        File tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    Int ram_gb = 8
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_gb} - 2 ))
        REGIONS=""
        for i in $(seq 1 22) X Y M; do
            REGIONS="${REGIONS} chr${i}"
        done
        
        # - Restricting to standard chromosomes
        # - Removing BND calls
        # - Ensuring that the input file is sorted
        ${TIME_COMMAND} bcftools view --include "SVTYPE != 'BND'" ~{vcf_gz} $(echo ${REGIONS}) | bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z > tmp0.vcf.gz
        tabix -f tmp0.vcf.gz
        
        # - Normalizing multiallelic records
        # - Fixing wrong REF values (which may occur e.g. in sniffles).
        ${TIME_COMMAND} bcftools norm --multiallelics -any --check-ref s --do-not-normalize --fasta-ref ~{reference_fa} --output-type z tmp0.vcf.gz > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        rm -f tmp0.vcf.gz*
        
        # - Removing identical records
        #   See <https://github.com/samtools/bcftools/issues/1089>.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --rm-dup exact --output-type z tmp1.vcf.gz > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        rm -f tmp1.vcf.gz*
        
        # - Removing INFO/STRAND, to avoid the following error when merging
        #   e.g. sniffles and cutesv VCFs with bcftools merge:
        #   [W::bcf_hdr_merge] Trying to combine "STRAND" tag definitions of
        #   different lengths
        #   Error occurred while processing INFO tag "STRAND" at chr1:...
        bcftools annotate --remove INFO/STRAND --output-type z tmp2.vcf.gz > tmp3.vcf.gz
        tabix -f tmp3.vcf.gz
        rm -f tmp2.vcf.gz*
        
        # - Additional fixes
        ${TIME_COMMAND} python ~{docker_dir}/resolve_light.py tmp3.vcf.gz ~{reference_fa} > tmp4.vcf
        rm -f tmp3.vcf.gz*
        
        bgzip -c tmp4.vcf > ~{sample_id}_resolved.vcf.gz
        tabix -f ~{sample_id}_resolved.vcf.gz
    >>>
    
    output {
    	File resolved_vcf_gz = "~{work_dir}/~{sample_id}_resolved.vcf.gz"
    	File resolved_tbi = "~{work_dir}/~{sample_id}_resolved.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: 1
        memory: ram_gb + "GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
