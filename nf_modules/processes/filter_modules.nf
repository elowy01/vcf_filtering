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