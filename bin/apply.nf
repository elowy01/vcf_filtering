/* 
 * Apply the fitted ML model on a VCF for filtering
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

nextflow.enable.dsl=2
include { GET_HEADER; SELECT_REGION; REHEADER; SPLIT_MULTIALLELIC; RUN_BCFTOOLS_SORT} from "${params.NF_MODULES}/processes/bcftools.nf"
include { BCFT_ANNOTATE } from "${params.NF_MODULES}/processes/bcftools.nf"
include { EXC_NON_VTS; DROP_GTPS; BCFT_QUERY } from "${params.NF_MODULES}/processes/bcftools.nf"
include { SELECT_VARIANTS as SELECT_SNPS} from "${params.NF_MODULES}/processes/bcftools.nf"
include { SELECT_VARIANTS as SELECT_INDELS} from "${params.NF_MODULES}/processes/bcftools.nf"
include { ALLELIC_PRIMITIVES; RUN_VT_UNIQ } from "${params.NF_MODULES}/processes/normalization.nf"
include { RUN_TABIX } from "${params.NF_MODULES}/processes/utils.nf"
include { MODIFY_HEADER; APPLY_MODEL; COMPRESS_PREDICTIONS } from "../nf_modules/processes/filter_modules.nf"
include { SAVE_FILE as SAVE_VCF} from "${params.NF_MODULES}/processes/utils.nf"
include { SAVE_FILE as SAVE_VCF_IX} from "${params.NF_MODULES}/processes/utils.nf"
 
// params defaults
params.help = false
params.threads = 1
params.outdir = './results'
params.tmpdir = '/tmp/'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to apply a certain fitted Logistic Regression model to filter a certain VCF file'
    log.info '-----------------------------------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow apply.nf --vcf VCF --model MODEL.sav --cutoff 0.95 --vt snps'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--vcf VCF    Path to the VCF file that will be filtered.'
    log.info '  --model FILE Path to serialized ML fitted model.'
    log.info '  --cutoff FLOAT Cutoff value used in the filtering. i.e.: 0.95.'
    log.info '  --annotations ANNOTATION_STRING String containing the annotations to filter, for example:'
    log.info '    %CHROM\t%POS\t%INFO/DP\t%INFO/RPB\t%INFO/MQB\t%INFO/BQB\t%INFO/MQSB\t%INFO/SGB\t%INFO/MQ0F\t%INFO/ICB\t%INFO/HOB\t%INFO/MQ\n.'
    log.info '  --chr chr1   Chromosome to be analyzed'
    log.info '  --vt  VARIANT_TYPE   Type of variant to filter. Poss1ible values are 'snps'/'indels'.'
    log.info '  --outdir OUTDIR   Name for directory for saving the output files of this pipeline. Default: results/'
    log.info '  --outprefix OUTPREFIX Output filename.'
    log.info ''
    exit 1
}

log.info 'Starting the process.....'

vcfFile = file(params.vcf)
if( !vcfFile.exists() ) {
  exit 1, "The specified input VCF file does not exist: ${params.vcf}"
}

workflow NORMALIZATION {
  /*
  This subworkflow will do the following:
  1) Normalize the VCF
  2) Split the normalized VCF into 2 VCF files (with SNPs and INDELs respectively)
  3) Remove the genotypes from either the SNP or INDEL vcf depending on params.vt
  */

  take: vcf
  main:
    SPLIT_MULTIALLELIC(vcf, params.threads )
    ALLELIC_PRIMITIVES(SPLIT_MULTIALLELIC.out)
    RUN_BCFTOOLS_SORT(ALLELIC_PRIMITIVES.out, params.tmpdir)
    SELECT_SNPS(RUN_BCFTOOLS_SORT.out, 'snps', params.threads)
    SELECT_INDELS(RUN_BCFTOOLS_SORT.out, 'indels', params.threads)
    DROP_GTPS(RUN_BCFTOOLS_SORT.out, params.vt, params.threads)
    RUN_VT_UNIQ(DROP_GTPS.out)
  emit:
    norm_no_gtps = RUN_VT_UNIQ.out
    norm_snps = SELECT_SNPS.out
    norm_indels = SELECT_INDELS.out
}

workflow APPLY_MODEL_WF {
  /*
  This subworkflow is used to apply the fitted model on the VCF
  */
  take: annotations
  take: vcf
  main:
    APPLY_MODEL(annotations, params.model, params.cutoff )
    COMPRESS_PREDICTIONS(APPLY_MODEL.out)
    BCFT_ANNOTATE(vcf,
                  COMPRESS_PREDICTIONS.out.predictions,
                  COMPRESS_PREDICTIONS.out.predictions_tbi,
                  'CHROM,POS,FILTER,prob_TP',
                  params.threads)
    RUN_BCFTOOLS_SORT(BCFT_ANNOTATE.out, params.tmpdir)
    RUN_TABIX(RUN_BCFTOOLS_SORT.out)
  emit:
    apply_vcf = RUN_BCFTOOLS_SORT.out
    apply_vcf_ix = RUN_TABIX.out
}

workflow  {
  main:
    RUN_TABIX(vcfFile)
    GET_HEADER(vcfFile)
    MODIFY_HEADER(GET_HEADER.out)
    SELECT_REGION(vcfFile, RUN_TABIX.out, params.chr)
    REHEADER(MODIFY_HEADER.out, SELECT_REGION.out)
    NORMALIZATION(REHEADER.out)
    EXC_NON_VTS(NORMALIZATION.out.norm_no_gtps, params.threads)
    BCFT_QUERY(EXC_NON_VTS.out, params.annotations)
    if ("${params.vt}"=="snps") {
      APPLY_MODEL_WF(BCFT_QUERY.out,NORMALIZATION.out.norm_snps)
    } else if ("${params.vt}"=="indels") {
      APPLY_MODEL_WF(BCFT_QUERY.out,NORMALIZATION.out.norm_indels)
    }
    SAVE_VCF(APPLY_MODEL_WF.out.apply_vcf, params.outdir,"${params.outprefix}.vcf.gz", 'copy')
    SAVE_VCF_IX(APPLY_MODEL_WF.out.apply_vcf_ix, params.outdir,"${params.outprefix}.vcf.gz.tbi", 'copy')
}