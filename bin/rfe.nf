/* 
 * Train ML model apply in the variant filtering
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

nextflow.enable.dsl=2
include { ALLELIC_PRIMITIVES; RUN_VT_UNIQ } from "${params.NF_MODULES}/processes/normalization.nf"
include { SAVE_FILE} from "${params.NF_MODULES}/processes/utils.nf"
include { SPLIT_MULTIALLELIC; SELECT_VARIANTS; RUN_BCFTOOLS_SORT; EXC_NON_VTS; INTERSECTION_CALLSET} from "${params.NF_MODULES}/processes/bcftools.nf"
include { BCFT_QUERY as BCFT_QUERY_TP} from "${params.NF_MODULES}/processes/bcftools.nf"
include { BCFT_QUERY as BCFT_QUERY_FP} from "${params.NF_MODULES}/processes/bcftools.nf"
include { RFE } from "../nf_modules/processes/filter_modules.nf"

// params defaults
params.help = false
params.threads = 1
params.outdir = './results'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to train a ML model'
    log.info '----------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow train.nf --vcf VCF --vt snps --outprefix out'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--vcf VCF    Path to the VCF file TODO.'
    log.info '  --vt  VARIANT_TYPE   Type of variant that will be normalized. Poss1ible values are \'snps\'/\'indels\'/\'both\'.'
    log.info '  --true_cs VCF  Path to the VCF file containing the gold-standard sites.'
    log.info '  --no_feats INT  Number of features (variant annotations) to select using RFE.'
    log.info '  --outdir OUTDIR   Name for directory for saving the output files of this pipeline. Default: results/'
    log.info '  --outprefix OUTPREFIX Output filename.'
    log.info ''
    exit 1
}

log.info 'Starting the process.....'

// Parameter validation

vcfFile = file(params.vcf)
if( !vcfFile.exists() ) {
  exit 1, "The specified input VCF file does not exist: ${params.vcf}"
}

true_vcf = file(params.true_cs)
if( !true_vcf.exists() ) {
  exit 1, "The specified gold-standard VCF file does not exist: ${params.true_cs}"
}
true_vcf_ix = file("${true_vcf}.tbi")
if( !true_vcf_ix.exists() ) {
  exit 1, "I need an index for: ${params.true_cs}"
}
if(!params.outprefix)
{ exit 1, "You need to provide the --outprefix parameter" }


workflow  {
    main:
        SPLIT_MULTIALLELIC( vcfFile, params.threads )
        ALLELIC_PRIMITIVES(SPLIT_MULTIALLELIC.out)
        SELECT_VARIANTS(ALLELIC_PRIMITIVES.out, params.vt, params.threads)
        RUN_BCFTOOLS_SORT(SELECT_VARIANTS.out)
        RUN_VT_UNIQ(RUN_BCFTOOLS_SORT.out)
        EXC_NON_VTS(RUN_VT_UNIQ.out, params.threads)
        INTERSECTION_CALLSET(EXC_NON_VTS.out, params.vt, true_vcf, true_vcf_ix)
        BCFT_QUERY_TP(INTERSECTION_CALLSET.out.tp_vcf, '%CHROM\t%POS\t%INFO\n')
        BCFT_QUERY_FP(INTERSECTION_CALLSET.out.fp_vcf, '%CHROM\t%POS\t%INFO\n')
        RFE(BCFT_QUERY_TP.out, BCFT_QUERY_FP.out, params.no_feats)
        SAVE_FILE(RFE.out.sel_feats, params.outdir, params.outprefix, 'move')
}
