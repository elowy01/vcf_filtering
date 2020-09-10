/* 
 *  Module containing processes for normalizing a VCF
 *
 * This module relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

process ALLELIC_PRIMITIVES {
        /*
        Process to run vcflib vcfallelicprimitives to decompose of MNPs

        Returns
        -------
        Path to decomposed VCF
        */

        executor 'local'

        input:
	path vcf

        output:
        path "out.decomp.vcf.gz"

        """
        tabix -f $vcf
        vcfallelicprimitives -k -g $vcf |bgzip -c > out.decomp.vcf.gz
        """
}

process RUN_VT_UNIQ {
        /*
        Process to run vt uniq
        */

        executor 'local'

        input:
        path vcf

        output:
        file "out.norm.vcf.gz"

        """
        vt uniq $vcf | bgzip -c > out.norm.vcf.gz
        """
}