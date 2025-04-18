nextflow.enable.dsl=2

process SNP_BATCHES {
    input:
    path bim

    output:
    path "snp_batch_ranges.txt", emit: batch_file

    script:
    """
    wc -l < ${bim} | awk '{
        batch_size = 500000;
        batch_num = 1;
        for (start=1; start <= \$1; start+=batch_size) {
            end = (start + batch_size - 1 > \$1) ? \$1 : start + batch_size - 1;
            print start "," end "," batch_num;
            batch_num++;
        }
    }' > snp_batch_ranges.txt
    """
}

process RUN_BATCH {
    input:
    tuple val(start), val(end), val(batch_num), path(pheno_file), path(r_script)

    output:
    file "CRA_AKClusters_batch${batch_num}.tsv"

    script:
    """
    source /etc/profile.d/modules.sh
    module purge
    module load R/4.3.2
    export TMPDIR="/d/tmp/repsh"
    Rscript ${r_script} \\
        ${params.plink_prefix} \\
        ${pheno_file} \\
        ${start} \\
        ${end} \\
        ${batch_num}
    """
}

workflow {
    // Define channels here from params
    def bim_ch   = Channel.value(params.bim_file)
    def pheno_ch = Channel.value(params.pheno_file)
    def rscript_ch = Channel.value(params.r_script)

        // Get the output file from the SNP_BATCHES process
    def batch_file_ch = SNP_BATCHES(bim_ch)

    batch_file_ch
        .splitText()
        .map { line ->
        	def (start, end, batch_num) = line.trim().split(',')
        	[start, end, batch_num.toInteger()]
    	}
	.combine(pheno_ch)
        .combine(rscript_ch)
        | RUN_BATCH
}

