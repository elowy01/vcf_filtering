/* 
 *  Module to normalize a VCF
 *
 * This module relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

 process split_multiallelic {
        /*
        This process will split the multiallelic variants by using BCFTools

        Returns
        -------
        Path to splitted VCF
        */

        memory '10 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

        output:
        file "out.splitted.vcf.gz" into out_splitted

        when:
        !params.region

        """
        bcftools norm -m -any ${params.vcf} -o out.splitted.vcf.gz -Oz --threads ${params.threads}
        """
}
