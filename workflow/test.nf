// Smoke-test workflow: runs snp_anova.R on the tiny example dataset.
// Does NOT use module load or cluster-specific settings.
// Requires: R with data.table and snpStats installed and Rscript on PATH.
// Run from the workflow/ directory:
//   nextflow run test.nf -params-file params-test.yaml

nextflow.enable.dsl=2

process SMOKE_TEST {
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
    SMOKE_TEST(
        Channel.value([
            1,
            params.end_snp,
            1,
            file(params.pheno_file),
            file(params.r_script)
        ])
    )
    | view { "Output: $it" }
}
