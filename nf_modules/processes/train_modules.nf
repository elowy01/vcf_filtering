process TRAIN_MODULE {
    /*
    Process to train the Logistic Regression binary classifier
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