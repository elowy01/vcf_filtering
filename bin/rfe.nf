/* 
 * Pipeline to perform Recursive Feature Elimination (RFE) to select the set of 
 * variant annotations that are more relevant to classify the genomic sites as 
 * variants / non-variants.
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

nextflow.enable.dsl=2
include { ALLELIC_PRIMITIVES; RUN_VT_UNIQ } from "../nf_modules/processes/normalization.nf"
include { SAVE_FILE} from "../nf_modules/processes/utils.nf"
include { SPLIT_MULTIALLELIC; SELECT_VARIANTS; RUN_BCFTOOLS_SORT; EXC_NON_VTS; INTERSECTION_CALLSET} from "../nf_modules/processes/bcftools.nf"
include { BCFT_QUERY as BCFT_QUERY_TP} from "../nf_modules/processes/bcftools.nf"
include { BCFT_QUERY as BCFT_QUERY_FP} from "../nf_modules/processes/bcftools.nf"
include { RFE } from "../nf_modules/processes/filter_modules.nf"

// params defaults
params.help = false
params.threads = 1
params.outdir = './results'
params.tmpdir = '/tmp/'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to perform Recursive Feature Elimination (RFE) on a VCF file'
    log.info '---------------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow rfe.nf --vcf VCF --vt snps --true_cs TRUE_VCF --no_feats 4 --outfile out_select_feats.txt'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--vcf VCF    Path to the VCF file.'
    log.info '  --vt  VARIANT_TYPE   Type of variant used in the analysis. Poss1ible values are \'snps\'/\'indels\'/\'both\'.'
    log.info '  --true_cs VCF  Path to the gold-standard VCF file.'
    log.info '  --no_feats INT  Number of features (variant annotations) to be selected using RFE.'
    log.info '  --outdir OUTDIR   Name of the directory used for saving the output files of this pipeline. Default: results/'
    log.info '  --outfile OUTFILE Output filename.'
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
if(!params.outfile)
{ exit 1, "You need to provide the --outfile parameter" }

workflow  {
    main:
        SPLIT_MULTIALLELIC( vcfFile, params.threads )
        ALLELIC_PRIMITIVES(SPLIT_MULTIALLELIC.out)
        SELECT_VARIANTS(ALLELIC_PRIMITIVES.out, params.vt, params.threads)
        RUN_BCFTOOLS_SORT(SELECT_VARIANTS.out, params.tmpdir)
        RUN_VT_UNIQ(RUN_BCFTOOLS_SORT.out)
        EXC_NON_VTS(RUN_VT_UNIQ.out, params.threads)
        INTERSECTION_CALLSET(EXC_NON_VTS.out, params.vt, true_vcf, true_vcf_ix)
        BCFT_QUERY_TP(INTERSECTION_CALLSET.out.tp_vcf, '%CHROM\t%POS\t%INFO\n')
        BCFT_QUERY_FP(INTERSECTION_CALLSET.out.fp_vcf, '%CHROM\t%POS\t%INFO\n')
        RFE(BCFT_QUERY_TP.out, BCFT_QUERY_FP.out, params.no_feats)
        SAVE_FILE(RFE.out.sel_feats, params.outdir, params.outfile, 'move')
}
