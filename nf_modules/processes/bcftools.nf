/* 
 *  Module containing processes for executing the different
 *  bcftools commands on a VCF
 *
 * This module relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

process SPLIT_MULTIALLELIC {
    /*
    This process will split the multiallelic variants by using BCFTools

    Returns
    -------
    Path to splitted VCF
    */
    input:
    path vcf
    val(threads)

    output:
    path 'out.bcftoolsnorm.vcf.gz'

    """
    bcftools norm -m -any ${vcf} -o out.bcftoolsnorm.vcf.gz -Oz --threads ${threads}
    """
}

process GET_HEADER {
    /*
    Process is used to get the header from a VCF file

    Returns
    -------
    Path to a text file with a header
    */
    executor 'local'

    input:
    path vcf

    output:
    path 'header.txt'

    """
    bcftools view -h ${vcf} > header.txt
    """
}

process REHEADER {
    /*
    Process to replace the header of a VCF file

    Parameters
    ----------
    header.txt : file with header
    vcf : vcf that will be altered

    Output
    ------
    gzipped file
    */
    executor 'local'

    input:
    path(header)
    path(vcf)

    output:
    path 'reheaded.vcf.gz'

    """
    bcftools reheader -h ${header} -o 'reheaded.vcf.gz' ${vcf}
    """
}

process SELECT_VARIANTS {
    /*
    Process to select the desired variants type (snps/indels)

    Parameters
    ----------
    vcf : input vcf
    vt : type of variant to select
    threads : Number of cpus to use
    */
    input:
    path(vcf)
    val(vt)
    val(threads)

    output:
    path "out.${vt}.vcf.gz"

    """
    bcftools view $vcf -v ${vt} -o out.${vt}.vcf.gz -Oz --threads ${threads}
    """
}

process RUN_BCFTOOLS_SORT {
    /*
    Process to run bcftools sort

    Parameters
    ----------
    vcf : vcf file to be sorted
    tmpdir : tmp dir keep the 'bcftools sort' intermediate files
    */
    input:
    path(vcf)
    val(tmpdir)

    output:
    path "out.sort.vcf.gz"

    """
	mkdir -p ${tmpdir}
	bcftools sort -T ${tmpdir} ${vcf} -o out.sort.vcf.gz -Oz
    """
}

process EXC_NON_VTS {
    /*
    This process will select the variants on a VCF

    Parameters
    ----------
    vcf : input vcf
    threads : Number of cpus to use

    Output
    -------
    Path to VCF containing only variants
    */
    input:
    path(vcf)
    val(threads)

    output:
    file "out.onlyvariants.vcf.gz"

    """
    bcftools view -c1 ${vcf} -o out.onlyvariants.vcf.gz --threads ${threads} -Oz
    """
}

process INTERSECTION_CALLSET {
    /*
    Process to find the intersection between a call set and the Gold
    standard call set

    Parameters
    ----------
    vcf : path to VCF file used for training
    vt : type of variant ('snps','indels')
    true_cs : path to gold standard call set
    true_cs_ix : path to index for 'true_cs'

    Output
    ------
    3 VCFs containing the False Positives in FP.vcf.gz, the False Negatives in FN.vcf.gz
    and the True Positives in TP.vcf.gz
    */
    input:
    path(vcf)
    val(vt)
    path(true_cs)
    path(true_cs_ix)

    output:
    path 'FP.vcf.gz', emit: fp_vcf
    path 'FN.vcf.gz', emit: fn_vcf 
    path 'TP.vcf.gz', emit: tp_vcf

    """
    tabix ${vcf}
    bcftools isec -c ${vt}  -p 'dir/' ${vcf} ${true_cs}
    bgzip -c dir/0000.vcf > FP.vcf.gz
    bgzip -c dir/0001.vcf > FN.vcf.gz
    bgzip -c dir/0002.vcf > TP.vcf.gz
    """
}

process BCFT_QUERY {
    /*
    Process to run 'bcftools query' on a VCF file

    Parameters
    ----------
    vcf : input vcf
    annotations : string with annotations to query

    Output
    ------
    .tsv file compressed with gzip containing the annotations
    that have been 
    */
    input:
    path(vcf)
    val(annotations)

    output:
    path("${vcf.baseName}.tsv.gz")

    """
    bcftools query -H -f '${annotations}' ${vcf} | bgzip -c > ${vcf.baseName}.tsv.gz
    """
}

process SELECT_REGION {
    /*
    Process to fetch a certain region from the VCF using 'bcftools view'

    Parameters
    ----------
    vcf : path to VCF file
    vcf_ix : path to tabix index for this VCF file
    region : region to fetch. For example: chr20
    threads : Number of cpus to use

    Output
    ------
    gzipped VCF file with a particular region
    */
    input:
    path(vcf)
    path(vcf_ix)
    val(region)
    val(threads)

    output:
    path("sub.${region}.vcf.gz")

    """
    bcftools view -r ${region} ${vcf} -o sub.${region}.vcf.gz ${threads} -Oz
    """
}

process DROP_GTPS {
    /*
    Process to drop the genotype information from a VCF
    Parameters
    ----------
    vcf : path to the VCF file
    vt : select/exclude comma-separated list of variant types: snps,indels,mnps,ref,bnd,other
    threads : Number of cpus to use

    Output
    ------
    gzipped VCF file without genotypes
    */
    input:
    path(vcf)
    val(vt)
    val(threads)

    output:
    path("out.${vt}.noGTPS.vcf.gz")

    """
    bcftools view -G -v ${vt} ${vcf} -o out.${vt}.noGTPS.vcf.gz -Oz --threads ${threads}
    """
}

process BCFT_ANNOTATE {
    /*
    Process to run 'bcftools annotate' on a VCF file

    Parameters
    ----------
    vcf : vcf file that will be reannotated
    tsv_f : tsv file with new annotations
    tsv_f_ix : tabix index for tsv_f
    columns : list of columns to use
    threads : Number of cpus to use

    Output
    ------
    vcf.gz : Reannotated vcf file 
    */
    input:
    path(vcf)
    path(tsv_f)
    path(tsv_f_ix)
    val(columns)
    val(threads)

    output:
    path("reannotated.vcf.gz")

    """
    bcftools annotate -a ${tsv_f} ${vcf} -c ${columns} -o reannotated.vcf.gz --threads ${threads} -Oz
    """
}