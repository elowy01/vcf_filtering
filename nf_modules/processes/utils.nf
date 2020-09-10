process SAVE_FILE {
    /*
    This process will split the multiallelic variants by using BCFTools

    Returns
    -------
    Path to splitted VCF
    */
    
    publishDir "results/", mode: 'copy', overwrite: true

    executor 'local'

    input:
	path vcf

    output:
    path "out.norm.vcf.gz"

    """
    ls $vcf
    """
}