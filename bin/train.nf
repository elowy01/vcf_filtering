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
include { ALLELIC_PRIMITIVES; RUN_VT_UNIQ } from "../nf_modules/processes/normalization.nf"
include { SAVE_FILE as SAVE_MODEL} from "../nf_modules/processes/utils.nf"
include { SAVE_FILE as SAVE_MODEL_SCORE} from "../nf_modules/processes/utils.nf"
include { SPLIT_MULTIALLELIC; SELECT_VARIANTS; RUN_BCFTOOLS_SORT; EXC_NON_VTS; INTERSECTION_CALLSET} from "../nf_modules/processes/bcftools.nf"
include { BCFT_QUERY as BCFT_QUERY_TP} from "../nf_modules/processes/bcftools.nf"
include { BCFT_QUERY as BCFT_QUERY_FP} from "../nf_modules/processes/bcftools.nf"
include { TRAIN_MODULE } from "../nf_modules/processes/filter_modules.nf"

// params defaults
params.help = false
params.threads = 1
params.outdir = './results'
params.tmpdir = '/tmp/'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to train a ML model'
    log.info '----------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow train.nf --vcf VCF --vt snps --annotations \'%CHROM\t%POS\t%INFO/DP\t%INFO/FS\n\' --outprefix out'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--vcf VCF    Path to the VCF file TODO.'
    log.info '  --vt  VARIANT_TYPE   Type of variant that will be normalized. Poss1ible values are \'snps\'/\'indels\'/\'both\'.'
    log.info '  --true_cs VCF  Path to the VCF file containing the gold-standard sites.'
    log.info '  --annotations ANNOTATION_STRING	String containing the annotations to use for training the model, for example:'
    log.info '	%CHROM\t%POS\t%INFO/DP\t%INFO/RPB\t%INFO/MQB\t%INFO/BQB\t%INFO/MQSB\t%INFO/SGB\t%INFO/MQ0F\t%INFO/ICB\t%INFO/HOB\t%INFO/MQ\n.'
    log.info '  --outprefix OUTPREFIX Output filename.'
    log.info '  --outdir OUTDIR   Name for directory for saving the output files of this pipeline. Default: results/'
    log.info '  --tmpdir TMPDIR   Name of temp directory used to store the temporary files. Default: /tmp/'
    log.info '  --threads INT   Number of CPUs to use. Default: 1.'
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
        RUN_BCFTOOLS_SORT(SELECT_VARIANTS.out, params.tmpdir)
        RUN_VT_UNIQ(RUN_BCFTOOLS_SORT.out)
        EXC_NON_VTS(RUN_VT_UNIQ.out, params.threads)
        INTERSECTION_CALLSET(EXC_NON_VTS.out, params.vt, true_vcf, true_vcf_ix)
        BCFT_QUERY_TP(INTERSECTION_CALLSET.out.tp_vcf, params.annotations)
        BCFT_QUERY_FP(INTERSECTION_CALLSET.out.fp_vcf, params.annotations)
        TRAIN_MODULE(BCFT_QUERY_TP.out, BCFT_QUERY_FP.out)
        SAVE_MODEL(TRAIN_MODULE.out.trained_model, params.outdir, "${params.outprefix}.sav", 'copy')
        SAVE_MODEL_SCORE(TRAIN_MODULE.out.trained_model_score, params.outdir, "${params.outprefix}.score", 'copy')
}
