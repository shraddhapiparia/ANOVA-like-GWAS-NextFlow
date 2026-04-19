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
    path "gwas_results_batch${batch_num}.csv"

    script:
    """
    Rscript ${r_script} \
        ${params.plink_prefix} \
        ${pheno_file} \
        ${start} \
        ${end} \
        ${batch_num}
    """
}

workflow {
    def bim_ch     = Channel.value(params.bim_file)
    def pheno_ch   = Channel.value(params.pheno_file)
    def rscript_ch = Channel.value(params.r_script)

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
