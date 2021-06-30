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
 include { GET_HEADER; SELECT_REGION; REHEADER; SPLIT_MULTIALLELIC} from "${params.NF_MODULES}/processes/bcftools.nf"
 include { MODIFY_HEADER } from "../nf_modules/processes/filter_modules.nf"
 

// params defaults
params.help = false
params.threads = 1
params.outdir = './results'

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
    log.info '  --cutoff FLOAT/FLOATs Cutoff value used in the filtering. It also accepts comma-separated list of floats: 0.95,0.96'
    log.info '  --annotations ANNOTATION_STRING String containing the annotations to filter, for example:'
    log.info '    %CHROM\t%POS\t%INFO/DP\t%INFO/RPB\t%INFO/MQB\t%INFO/BQB\t%INFO/MQSB\t%INFO/SGB\t%INFO/MQ0F\t%INFO/ICB\t%INFO/HOB\t%INFO/MQ\n.'
    log.info '  --chr chr1   Chromosome to be analyzed'
    log.info '  --vt  VARIANT_TYPE   Type of variant to filter. Poss1ible values are 'snps'/'indels'.'
    log.info ''
    exit 1
}

log.info 'Starting the process.....'

vcfFile = file(params.vcf)
if( !vcfFile.exists() ) {
  exit 1, "The specified input VCF file does not exist: ${params.vcf}"
}

workflow  {
    main:
        GET_HEADER(vcfFile)
        MODIFY_HEADER(GET_HEADER.out)
        SELECT_REGION(vcfFile, params.chr)
        REHEADER(MODIFY_HEADER.out, SELECT_REGION.out)
        SPLIT_MULTIALLELIC(REHEADER.out, params.threads )
}