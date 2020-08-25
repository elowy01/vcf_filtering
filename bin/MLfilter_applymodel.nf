/* 
 * Workflow to filter a VCF
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false
params.split_multiallelics = false
params.threads = 1
params.queue = 'production-rh74'
params.executor = 'local'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to apply a certain fitted Logistic Regression model to filter a certain VCF file'
    log.info '-----------------------------------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow MLfilter_applymodel.nf --vcf VCF --model MODEL.sav --cutoff 0.95 --threads 5 --vt snps'
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
    log.info '  --threads INT Number of threads used in the different BCFTools processes. Default=1.'
    log.info ''
    exit 1
}

log.info 'Starting the analysis.....'

//Apply a fitted model obtained after running MLfilter_trainmodel.nf

process get_header {
        /*
        Process to get the header of the unfiltered VCF
        */

        memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        output:
        file 'header.txt' into header

        """
        bcftools view -h ${params.vcf} > header.txt
        """
}

process modify_header {
        /*
        Process to modify the header of the unfiltered VCF
        */

        memory '500 MB'
        executor 'local'
        queue "${params.queue}"
        cpus 1

        input:
        file header

        output:
        file 'newheader.txt' into newheader

        """
        #!/usr/bin/env python

        from VCF.VcfUtils import VcfUtils

        vcf_object=VcfUtils(vcf='${params.vcf}')

        vcf_object.add_to_header(header_f='${header}', outfilename='newheader1.txt',
                                 line_ann='##FILTER=<ID=MLFILT,Description="Binary classifier filter">')
        vcf_object.add_to_header(header_f='newheader1.txt', outfilename='newheader.txt',
                                 line_ann='##INFO=<ID=prob_TP,Number=1,Type=Float,Description="Probability of being a True positive">')
        """
}

process splitVCF {
        /*
        This process will select a single chromosome from the VCF.
	It will also drop the genotypes from the splitted VCF
        */

        memory '5 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

        output:
        file "unfilt.${params.chr}.vcf.gz" into unfilt_vcf_chr

        """
        bcftools view -r ${params.chr} ${params.vcf} -o unfilt.${params.chr}.vcf.gz --threads ${params.threads} -Oz
        """
}

process replace_header {
        /*
        Process to replace header in the unfiltered VCF
        */

        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
        file newheader
        file unfilt_vcf_chr

        output:
        file 'unfilt_reheaded.vcf.gz' into unfilt_vcf_chr_reheaded

        """
        bcftools reheader -h ${newheader} -o 'unfilt_reheaded.vcf.gz' ${unfilt_vcf_chr}
        """
}

// normalize vcf

process split_multiallelic {
        /*
        This process will split the multiallelic variants by using BCFTools
        Returns
        -------
        Path to splitted VCF
        */

        memory { 5.GB * task.attempt }
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

	errorStrategy 'retry'
        maxRetries 5

	input:
	file unfilt_vcf_chr_reheaded from unfilt_vcf_chr_reheaded

        output:
        file "out.splitted.vcf.gz" into out_splitted

        """
        bcftools norm -m -any ${unfilt_vcf_chr_reheaded} -o out.splitted.vcf.gz -Oz --threads ${params.threads}
        """
}

process allelic_primitives {
        /*
        Process to run vcflib vcfallelicprimitives to decompose of MNPs

        Returns
        -------
        Path to decomposed VCF
        */

	memory { 8.GB * task.attempt }
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        errorStrategy 'retry'
        maxRetries 5

        input:
	file out_splitted from out_splitted

        output:
        file "out.splitted.decomp.vcf.gz" into out_decomp

        """
        tabix -f ${out_splitted}
        vcfallelicprimitives -k -g ${out_splitted} |bgzip -c > out.splitted.decomp.vcf.gz
        """
}

process sort_out_decomp {
	/*
	Process to bcftools sort the 'out_decomp'
	*/

	memory '9 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

	input:
	file out_decomp

	output:
	file "out.splitted.decomp.sorted.vcf.gz" into out_decomp_sorted
	
	"""
	mkdir -p tmpdir
        bcftools sort -T tmpdir/ ${out_decomp} -o out.splitted.decomp.sorted.vcf.gz -Oz
	"""
}

process select_variants {
        /*
        Process to select the desired variants type (snps/indels)
        */

        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

        input:
        file out_decomp_sorted

        output:
	set file("out.snps.gtps.vcf.gz"), file("out.indels.gtps.vcf.gz") into gtps_vts
        file "out.${params.vt}.vcf.gz" into no_gtps_vts

        """
	bcftools view -v snps ${out_decomp_sorted} -o out.snps.gtps.vcf.gz --threads ${params.threads} -Oz
        bcftools view -v indels ${out_decomp_sorted} -o out.indels.gtps.vcf.gz --threads ${params.threads} -Oz
	bcftools view -G -v ${params.vt} ${out_decomp_sorted} -o out.${params.vt}.vcf.gz -O z --threads ${params.threads}
        """	 
}

process run_vt_uniq {
        /*
        Process to run vt uniq
        Returns
        -------
        Path to final normalized file
        */

        memory '9 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
        file no_gtps_vts

        output:
        file "out.normalized.vcf.gz" into out_uniq

        """
        vt uniq ${no_gtps_vts} | bgzip -c > out.normalized.vcf.gz
        """
}

process excludeNonVariants {
        /*
        This process will select the variants on the unfiltered vcf (all chros) for
        the particular type defined by 'params.vt'
        Returns
        -------
        Path to a site-VCF file containing just the variants on a particular chromosome
        */

        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

        input:
        file out_uniq

        output:
        file "out.onlyvariants.vcf.gz" into out_onlyvariants

        """
        bcftools view -c1 ${out_uniq} -o out.onlyvariants.vcf.gz --threads ${params.threads} -Oz
        """
}

// start the application of the model

process get_variant_annotations {
	/*
	Process to get the variant annotations for the selected ${params.vt} from the unfiltered VCF file
	*/
	
	memory '2 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

	input:
	file out_onlyvariants

	output:
	file 'unfilt_annotations.vt.tsv.gz' into unfilt_annotations

	"""
	tabix ${out_onlyvariants}
	bcftools view -r ${params.chr} -v ${params.vt} ${out_onlyvariants} -o out.onlyvariants.vt.vcf.gz -Oz --threads ${params.threads}
	tabix out.onlyvariants.vt.vcf.gz
	bcftools query -H -r ${params.chr} -f '${params.annotations}' out.onlyvariants.vt.vcf.gz | bgzip -c > unfilt_annotations.vt.tsv.gz
	"""
}

cutoff_list = Channel.from( params.cutoff.split(',') )

cutoff_list.set { cutoff_values}

process apply_model {
	/*
	Process to read-in the serialized ML model created by running MLfilter_trainmodel.nf
	and to apply this model on the unfiltered VCF
	*/
	tag "Apply model with $cutoff"

	memory { 5.GB * task.attempt }
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

	errorStrategy 'retry'
        maxRetries 5

	input:
	file unfilt_annotations
	each cutoff from cutoff_values

	output:
	file 'predictions.tsv' into predictions
	val cutoff into cutoff_ch

	"""
	#!/usr/bin/env python

	from VCF.VCFfilter.MLclassifier import MLclassifier

	ML_obj=MLclassifier(fitted_model = '${params.model}')

	ML_obj.predict(outprefix="predictions", annotation_f='${unfilt_annotations}', cutoff=${cutoff})
	"""
}

process compress_predictions {
	/*
	Process to compress and index the 'predictions.tsv' file generated by process 'apply_model'
	*/
	tag "Compress predictions with $cutoff"

	memory '1.GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

	input:
	file predictions
	val cutoff from cutoff_ch

	output:
	set file('predictions.tsv.gz'), file('predictions.tsv.gz.tbi'), val(cutoff) into predictions_table

	"""
	bgzip -c ${predictions} > 'predictions.tsv.gz'
	tabix -f -s1 -b2 -e2 'predictions.tsv.gz'
	"""
}

reannotate_parameters = gtps_vts.combine(predictions_table)

process reannotate_vcf {
	/*
	Process to reannotate the unfiltered VCF with the information generated after applying the classifier
	*/
	tag "reannotate_vcf with $cutoff"

	memory { 8.GB * task.attempt }
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

        errorStrategy 'retry'
        maxRetries 5

	publishDir "results_${params.chr}", mode: 'copy', overwrite: true

	exec:
	def selected
	def non_selected

	if ("${params.vt}"=="snps") {
           selected="snps"
           non_selected="indels"
        } else if ("${params.vt}"=="indels") {
           selected="indels"
           non_selected="snps"
        }

	input:
        set file("out.snps.gtps.vcf.gz"), file("out.indels.gtps.vcf.gz"), file(predictions_table), file(predictions_table_tabix), val(cutoff) from reannotate_parameters

        output:
        file(output_cutoff)

        script:
        output_cutoff="filt.${cutoff}".replace('.', '_')+".vcf.gz"
        output_cutoff_tabix="filt.${cutoff}".replace('.', '_')+".vcf.gz.tbi"

	"""
        bcftools annotate -a ${predictions_table} out.${selected}.gtps.vcf.gz -c CHROM,POS,FILTER,prob_TP -o reannotated.vcf.gz --threads ${params.threads} -Oz
        bcftools concat reannotated.vcf.gz out.${non_selected}.gtps.vcf.gz -o out.merged.vcf.gz --threads ${params.threads} -Oz
        bcftools sort -T tmpdir/ out.merged.vcf.gz -o ${output_cutoff} -Oz
        tabix ${output_cutoff}
        """
}
