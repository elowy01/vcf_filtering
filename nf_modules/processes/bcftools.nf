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

    executor 'local'

    input:
    path vcf
    val(threads)


    output:
    path 'out.bcftoolsnorm.vcf.gz'

    """
    bcftools norm -m -any $vcf -o out.bcftoolsnorm.vcf.gz -Oz --threads $threads
    """
}

process SELECT_VARIANTS {
    /*
    Process to select the desired variants type (snps/indels)
    */

    executor 'local'

    input:
    path(vcf)
    val(vt)
    val(threads)

    output:
    path "out.${vt}.vcf.gz"

    """
    bcftools view $vcf -v $vt -o out.${vt}.vcf.gz -Oz --threads $threads
    """
}

process RUN_BCFTOOLS_SORT {
    /*
    Process to run bcftools sort
    */

    executor 'local'

    input:
    path(vcf)

    output:
    path "out.sort.vcf.gz"

    """
	mkdir tmpdir
	bcftools sort -T tmpdir/ ${vcf} -o out.sort.vcf.gz -Oz
    """
}