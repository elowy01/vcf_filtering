/* 
 * Workflow to normalize a VCF file
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

nextflow.enable.dsl=2

include { ALLELIC_PRIMITIVES; RUN_VT_UNIQ } from '../nf_modules/normalization.nf' 
include { SPLIT_MULTIALLELIC; SELECT_VARIANTS; RUN_BCFTOOLS_SORT } from '../nf_modules/bcftools.nf'
include { SAVE_FILE } from '../nf_modules/utils.nf'

log.info 'Starting the normalization for the desired variant type.....'

workflow NORM_VT  {
    main:
        SPLIT_MULTIALLELIC( params.vcf, params.threads )
        ALLELIC_PRIMITIVES(SPLIT_MULTIALLELIC.out)
        SELECT_VARIANTS(ALLELIC_PRIMITIVES.out, params.vt, params.threads)
        RUN_BCFTOOLS_SORT(SELECT_VARIANTS.out)
        RUN_VT_UNIQ(RUN_BCFTOOLS_SORT.out)
}

log.info 'normalization completed.....'
