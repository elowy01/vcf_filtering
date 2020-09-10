/* 
 * normalize variants (snps/indels) in a VCF file
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

nextflow.enable.dsl=2

include { NORM_SNPS } from '../nf_modules/subworkflows/normalization_snps.nf' 
include { SAVE_FILE } from '../nf_modules/processes/utils.nf'

// params defaults
params.help = false
params.threads = 1
params.outdir = './results'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to normalize a VCF file'
    log.info '--------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow normalization.nf --vcf VCF --vt snps --threads 5 --outprefix out'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--vcf VCF    Path to the VCF file that will be normalized.'
    log.info '  --vt  VARIANT_TYPE   Type of variant that will be normalized. Poss1ible values are 'snps'/'indels'/'both'.'
    log.info '  --threads INT Number of threads used in the different BCFTools processes. Default=1.'
    log.info '  --outdir OUTDIR Name for directory for saving the normalized output VCF.'
    log.info ''
    exit 1
}

log.info 'Starting the normalization.....'

workflow  {
    main:
        SPLIT_MULTIALLELIC( params.vcf, params.threads )
        ALLELIC_PRIMITIVES(SPLIT_MULTIALLELIC.out)
        SELECT_VARIANTS(ALLELIC_PRIMITIVES.out, params.vt, params.threads)
        RUN_BCFTOOLS_SORT(SELECT_VARIANTS.out)
        RUN_VT_UNIQ(RUN_BCFTOOLS_SORT.out)
        SAVE_FILE(RUN_VT_UNIQ.out)
}

