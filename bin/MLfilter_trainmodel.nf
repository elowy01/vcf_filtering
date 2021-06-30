/* 
 * Workflow to Train a Logistic Regression ML model that can be used for filtering a variant call set
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false
params.threads = 1
params.queue = 'production-rh74'
params.rfe = false
params.region = false

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to train a Logistic Regression binary classifier'
    log.info '---------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow MLfilter_trainmodel.nf --vcf VCF --true VCF --vt snps --annotations ANNOTATION_STRING --threads 5 --outprefix out'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--vcf VCF    Path to the VCF file that will be used for training.'
    log.info '  --true VCF  Path to the VCF file containing the gold-standard sites.'
    log.info '  --vt  VARIANT_TYPE   Type of variant to use for training the model. Poss1ible values are 'snps'/'indels'.'
    log.info '  --region chr1,chr2  Comma separated list of chromosomes to analyse. If omitted then all chros will be analysed.'
    log.info '  --annotations ANNOTATION_STRING	String containing the annotations to filter, for example:'
    log.info '	%CHROM\t%POS\t%INFO/DP\t%INFO/RPB\t%INFO/MQB\t%INFO/BQB\t%INFO/MQSB\t%INFO/SGB\t%INFO/MQ0F\t%INFO/ICB\t%INFO/HOB\t%INFO/MQ\n.' 
    log.info '  --threads INT Number of threads used in the different BCFTools processes. Default=1.'
    log.info '  --outprefix OUTPREFIX Prefix for output files.'
    log.info ''
    exit 1
}

log.info 'Starting the analysis.....'

chrs_splitmultiallelic_chr=Channel.empty()

if (params.region) {
    log.info '\t--region provided'

    chrList = Channel.from( params.region.split(',') )

    chrList.set { chrs_splitmultiallelic_chr}
}

// Normalization
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

process allelic_primitives {
        /*
        Process to run vcflib vcfallelicprimitives to decompose of MNPs

        Returns
        -------
        Path to decomposed VCF
        */

        memory '9 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
	file out_splitted from out_splitted

        output:
        file "out.splitted.decomp.vcf.gz" into out_decomp

	when:
        !params.region

        """
        tabix -f ${out_splitted}
        vcfallelicprimitives -k -g ${out_splitted} |bgzip -c > out.splitted.decomp.vcf.gz
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
        file out_decomp

        output:
        file "out.${params.vt}.vcf.gz" into out_vts

	when:
        !params.region

        """
        bcftools view -v ${params.vt} ${out_decomp} -o out.${params.vt}.vcf.gz -O z --threads ${params.threads}
        """
}

process run_bcftools_sort {
        /*
        Process to run bcftools sort

        Returns
        -------
        Path to sorted VCF
        */

        memory '9 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
        file out_vts

        output:
        file "${params.outprefix}.sort.vcf.gz" into out_sort

	when:
        !params.region

        """
	mkdir tmpdir
	bcftools sort -T tmpdir/ ${out_vts} -o ${params.outprefix}.sort.vcf.gz -Oz
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
        file out_sort

        output:
        file "${params.outprefix}.normalized.vcf.gz" into out_uniq

	when:
        !params.region

        """
        vt uniq ${out_sort} | bgzip -c > ${params.outprefix}.normalized.vcf.gz
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

	when:
        !params.region

        """
        bcftools view -c1 ${out_uniq} -o out.onlyvariants.vcf.gz --threads ${params.threads} -Oz
        """
}

process intersecionCallSets {
        /*
        Process to find the intersection between out_sites_vts and the Gold standard call set
        */

        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
        file out_onlyvariants from out_onlyvariants

        output:
        file 'dir/' into out_intersect

	when:
        !params.region

        """
        tabix ${out_onlyvariants}
        bcftools isec -c ${params.vt}  -p 'dir/' ${out_onlyvariants} ${params.true}
        """
}

process compressIntersected {
        /*
        Process to compress the files generated by bcftools isec
        */
        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
	file out_intersect from out_intersect

        output:
        file 'FP.vcf.gz' into fp_vcf
        file 'FN.vcf.gz' into fn_vcf
        file 'TP.vcf.gz' into tp_vcf

	when:
        !params.region

        """
        bgzip -c ${out_intersect}/0000.vcf > FP.vcf.gz
        bgzip -c ${out_intersect}/0001.vcf > FN.vcf.gz
        bgzip -c ${out_intersect}/0002.vcf > TP.vcf.gz
        """
}

process get_variant_annotations {
        /*
        Process to get the variant annotations for training files
        and for VCF file to annotate (for a single chromosome in this case)
        */

        memory '2 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
        file tp_vcf
        file fp_vcf

        output:
        file 'TP_annotations.tsv.gz' into tp_annotations_train, tp_annotations_rfe
        file 'FP_annotations.tsv.gz' into fp_annotations_train, fp_annotations_rfe

	when:
        !params.region

        """
        bcftools query -H -f '${params.annotations}' ${tp_vcf} | bgzip -c > TP_annotations.tsv.gz
        bcftools query -H -f '${params.annotations}' ${fp_vcf} | bgzip -c > FP_annotations.tsv.gz
        """
}

process train_model {
        /*
        Process that takes TP_annotations.tsv and FP_annotations.tsv created above and will train the Logistic
        Regression binary classifier
        */

        memory '15 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        publishDir "trained_model", mode: 'copy', overwrite: true

        input:
        file tp_annotations from tp_annotations_train
        file fp_annotations from fp_annotations_train

        output:
        file 'fitted_logreg_vts.sav' into trained_model
        file 'fitted_logreg_vts.score' into trained_model_score

        when:
        !params.rfe && !params.region

        """
        #!/usr/bin/env python

        from VCF.VCFfilter.MLclassifier import MLclassifier

        ML_obj=MLclassifier()

        outfile=ML_obj.train(outprefix="fitted_logreg_vts",
                        tp_annotations='${tp_annotations}',
                        fp_annotations='${fp_annotations}')
        """
}
 
        process rfe {
        /*
        Process to do Recursive Feature Elimination in order to
        select the annotations that are more relevant for prediction'
        */

        memory '15 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        publishDir "selected_feats", mode: 'copy', overwrite: true

        input:
        file tp_annotations from tp_annotations_rfe
        file fp_annotations from fp_annotations_rfe

        output:
        file 'selected_feats.txt' into selected_feats

        when:
        params.rfe && !params.region

        """
        #!/usr/bin/env python
        from VCF.VCFfilter.MLclassifier import MLclassifier

        ML_obj=MLclassifier()

        outfile=ML_obj.rfe(
                        tp_annotations='${tp_annotations}',
                        fp_annotations='${fp_annotations}',
                        n_features='${params.no_features}',
                        outreport="selected_feats.txt")
        """
        }

process split_multiallelic_chr {
        /*
        This process will split the multiallelic variants for a particular chr by using BCFTools

        Returns
        -------
        Path to splitted VCF
        */

	tag {"Processing: "+chr}

        memory '10 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

        input:
        val chr from chrs_splitmultiallelic_chr

        output:
        file "out.splitted.${chr}.vcf.gz" into out_splitted_chr
	val chr into chr1

	when:
	params.region

        """
        bcftools norm -r ${chr} -m -any ${params.vcf} -o out.splitted.${chr}.vcf.gz -Oz --threads ${params.threads}
        """
}

process allelic_primitives_chr {
        /*
        Process to run vcflib vcfallelicprimitives to decompose MNPs for a particular chr

        Returns
        -------
        Path to decomposed VCF
        */

	tag {"Processing: "+chr}

        memory '9 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
	file out_splitted from out_splitted_chr
	val chr from chr1

        output:
        file "out.splitted.decomp.vcf.gz" into out_decomp_chr
	val chr into chr2

	when:
	params.region

        """
        tabix -f ${out_splitted}
        vcfallelicprimitives -k -g ${out_splitted} |bgzip -c > out.splitted.decomp.vcf.gz
        """
}

process select_variants_chr {
        /*
        Process to select the desired variants type (snps/indels) for a particular chr
        */
	
	tag {"Processing: "+chr}

        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

        input:
        file out_decomp_chr
	val chr from chr2

        output:
        file "out.${params.vt}.vcf.gz" into out_vts_chr
	val chr into chr3
	
	when:
	params.region

        """
        bcftools view -v ${params.vt} ${out_decomp_chr} -o out.${params.vt}.vcf.gz -O z --threads ${params.threads}
        """
}


process run_bcftools_sort_chr {
        /*
        Process to run bcftools sort for a particular chr

        Returns
        -------
        Path to sorted VCF
        */

	tag {"Processing: "+chr}

        memory '9 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
        file out_vts_chr
	val chr from chr3

        output:
        file "${params.outprefix}.sort.vcf.gz" into out_sort_chr
	val chr into chr4

	when:
	params.region

        """
	mkdir tmpdir
	bcftools sort -T tmpdir/ ${out_vts_chr} -o ${params.outprefix}.sort.vcf.gz -Oz
        """
}

process run_vt_uniq_chr {
        /*
        Process to run vt uniq for a particular region

        Returns
        -------
        Path to final normalized file for a particular chr
        */

	tag {"Processing: "+chr}

        memory '9 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
        file out_sort_chr
	val chr from chr4

        output:
        file "${params.outprefix}.normalized.vcf.gz" into out_uniq_chr
	val chr into chr5

	when:
	params.region

        """
        vt uniq ${out_sort_chr} | bgzip -c > ${params.outprefix}.normalized.vcf.gz
        """
}

process excludeNonVariants_chr {
        /*
        This process will select the only the variants for a particular chr
        */

	tag {"Processing: "+chr}

        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus "${params.threads}"

        input:
        file out_uniq_chr
	val chr from chr5

        output:
        file "out.onlyvariants.vcf.gz" into out_onlyvariants_chr
	val chr into chr6

	when:
	params.region

        """
        bcftools view -c1 ${out_uniq_chr} -o out.onlyvariants.vcf.gz --threads ${params.threads} -Oz
        """
}

process intersecionCallSets_chr {
        /*
        Process to find the intersection between out_sites_vts and the Gold standard call set
	for a particular region
        */

	tag {"Processing: "+chr}	

        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
        file out_onlyvariants from out_onlyvariants_chr
        val chr from chr6

        output:
        file 'dir/' into out_intersect_chr
	val chr into chr7

	when:
        params.region

        """
        tabix ${out_onlyvariants}
        bcftools isec -r ${chr} -c ${params.vt}  -p 'dir/' ${out_onlyvariants} ${params.true}
        """
}

process compressIntersected_chr {
        /*
        Process to compress the files generated by bcftools isec for a particular chr
        */
	
	tag {"Processing: "+chr}

        memory '500 MB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
	file out_intersect from out_intersect_chr
	val chr from chr7

        output:
        file 'FP.vcf.gz' into fp_vcf_chr
        file 'FN.vcf.gz' into fn_vcf_chr
        file 'TP.vcf.gz' into tp_vcf_chr
	val chr into chr8

	when:
        params.region

        """
        bgzip -c ${out_intersect}/0000.vcf > FP.vcf.gz
        bgzip -c ${out_intersect}/0001.vcf > FN.vcf.gz
        bgzip -c ${out_intersect}/0002.vcf > TP.vcf.gz
        """
}

process get_variant_annotations_chr {
        /*
        Process to get the variant annotations for training files
        and for VCF file to annotate (for a single chromosome in this case)
        */

	tag {"Processing: "+chr}

        memory '2 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        input:
        file tp_vcf_chr
        file fp_vcf_chr
	val chr from chr8

        output:
        file 'TP_annotations.tsv.gz' into tp_annotations_train_chr, tp_annotations_rfe_chr
        file 'FP_annotations.tsv.gz' into fp_annotations_train_chr, fp_annotations_rfe_chr
	val chr into train_model_chr, rfe_chr

	when:
        params.region

        """
        bcftools query -H -f '${params.annotations}' ${tp_vcf_chr} | bgzip -c > TP_annotations.tsv.gz
        bcftools query -H -f '${params.annotations}' ${fp_vcf_chr} | bgzip -c > FP_annotations.tsv.gz
        """
}

process train_model_chr {
        /*
        Process that takes TP_annotations.tsv and FP_annotations.tsv created above and will train the Logistic
        Regression binary classifier for a particular chr
        */

	tag {"Processing: "+chr}

        memory '15 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        publishDir "trained_model_${chr}", mode: 'copy', overwrite: true

        input:
        file tp_annotations from tp_annotations_train_chr
        file fp_annotations from fp_annotations_train_chr
        val chr from train_model_chr

        output:
        file 'fitted_logreg_vts.sav' into trained_model_chr
        file 'fitted_logreg_vts.score' into trained_model_score_chr

        when:
        !params.rfe && params.region

        """
        #!/usr/bin/env python

        from VCF.VCFfilter.MLclassifier import MLclassifier

        ML_obj=MLclassifier()

        outfile=ML_obj.train(outprefix="fitted_logreg_vts",
                        tp_annotations='${tp_annotations}',
                        fp_annotations='${fp_annotations}')
        """
}
 
process rfe_chr {
        /*
        Process to do Recursive Feature Elimination in order to
        select the annotations that are more relevant for prediction'
        */

	tag {"Processing: "+chr}

        memory '15 GB'
        executor 'lsf'
        queue "${params.queue}"
        cpus 1

        publishDir "selected_feats_${chr}", mode: 'copy', overwrite: true

        input:
        file tp_annotations from tp_annotations_rfe_chr
        file fp_annotations from fp_annotations_rfe_chr
        val chr from rfe_chr

        output:
        file 'selected_feats.txt' into selected_feats_chr

        when:
        params.rfe && params.region

        """
        #!/usr/bin/env python
        from VCF.VCFfilter.MLclassifier import MLclassifier

        ML_obj=MLclassifier()

        outfile=ML_obj.rfe(
                        tp_annotations='${tp_annotations}',
                        fp_annotations='${fp_annotations}',
                        n_features='${params.no_features}',
                        outreport="selected_feats.txt")
        """
}
