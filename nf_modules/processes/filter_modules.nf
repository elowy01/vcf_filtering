process TRAIN_MODULE {
    /*
    Process to train the Logistic Regression binary classifier

    Parameters
    ----------
    tsv.gz file with annotations of true positive variants
    tsv.gz file with annotations of false positive variants
    */
    input:
    path tp_annotations
    path fp_annotations

    output:
    path 'fitted_logreg_vts.sav', emit: trained_model
    path 'fitted_logreg_vts.score', emit: trained_model_score

    """
    #!/usr/bin/env python

    from VCF.VCFfilter.MLclassifier import MLclassifier

    ML_obj=MLclassifier()

    outfile=ML_obj.train(outprefix="fitted_logreg_vts",
                        tp_annotations='${tp_annotations}',
                        fp_annotations='${fp_annotations}')
    """
}

process APPLY_MODEL {
    /*
    Process to apply the fitted module obtained using TRAIN_MODEL

    Parameters
    ----------
    tsv.gz file : Annotations obtained using the BCFT_QUERY process
    model file : File containing the serialized model obtained by running 'train.nf'
    cutoff value : Sensitivity cutoff for filtering

    Output
    ------
    .tsv file with predictions for each variant position
    */
    input:
    path annotations
    path model
    val(cutoff)

    output:
    path 'predictions.tsv'

    """
	#!/usr/bin/env python

	from VCF.VCFfilter.MLclassifier import MLclassifier

	ML_obj=MLclassifier(fitted_model = '${model}')

	ML_obj.predict(outprefix="predictions", annotation_f='${annotations}', cutoff=${cutoff})
	"""
}

process COMPRESS_PREDICTIONS {
    /*
    Process to compress the predictions generated by the process named 'APPLY_MODEL'
    It will also generate a tabix index for the .tsv file

    Parameters
    ----------
    .tsv : file with predictions

    Output
    ------
    .tsv.gz : compressed file with predictions
    .tsv.gz.tbi : tabix index of predictions
    */
    input:
    path predictions

    output:
    path 'predictions.tsv.gz', emit: predictions
    path 'predictions.tsv.gz.tbi', emit: predictions_tbi

    """
	bgzip -c ${predictions} > 'predictions.tsv.gz'
	tabix -f -s1 -b2 -e2 'predictions.tsv.gz'
	"""
}

process RFE {
    /*
    Process to do Recursive Feature Elimination in order to
    select the annotations that are more relevant for prediction

    Parameters
    ----------
    tsv.gz : file with annotations of true positive variants
    tsv.gz : file with annotations of false positive variants
    int : number of features to select
    */
    input:
    path tp_annotations
    path fp_annotations
    val(no_features)

    output:
    path 'selected_feats.txt', emit: sel_feats

    """
    #!/usr/bin/env python
    from VCF.VCFfilter.MLclassifier import MLclassifier

    ML_obj=MLclassifier()

    outfile=ML_obj.rfe(
                    tp_annotations='${tp_annotations}',
                    fp_annotations='${fp_annotations}',
                    n_features='${no_features}',
                    outreport="selected_feats.txt")
    """
}

process MODIFY_HEADER {
    /*
    Process to modify the header of the unfiltered VCF

    Parameters
    ----------
    header.txt : old header that will be modified
    */
    input:
    path header_f

    output:
    path 'newheader.txt'

     """
    #!/usr/bin/env python

    from VCF.VcfUtils import VcfUtils

    vcf_object=VcfUtils(vcf='${params.vcf}')

    vcf_object.add_to_header(header_f='${header_f}', outfilename='newheader1.txt',
                                 line_ann='##FILTER=<ID=MLFILT,Description="Binary classifier filter">')
    vcf_object.add_to_header(header_f='newheader1.txt', outfilename='newheader.txt',
                            line_ann='##INFO=<ID=prob_TP,Number=1,Type=Float,Description="Probability of being a True positive">')
    """
}