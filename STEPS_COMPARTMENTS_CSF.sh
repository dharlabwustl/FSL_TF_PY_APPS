echo " I AM WORKING"
## STEPS 1: Registration of template ct (COLIHM62) to a session CT.
#/software/lin_transform_before_deepreg_COLIHM62_as_moving.sh $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST /input1 /output
#/software/lin_transform_before_deepreg_COLIHM62_as_moving_mat_on_mask.sh $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST /input1 /output
#/software/deepregbasedregis_csf_cistern_midline_separation_with_COLIHM62_NO_RAPIDS.sh $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST /input1 /output
#/software/compartment_separation_with_vent_afterdeepreg_with_cistern_with_COLIHM62_with_midline.sh  $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST
#/software/csf_compartments_vols_N_display_for_NON_SAH.sh ${SESSION_ID} $XNAT_USER $XNAT_PASS $XNAT_HOST
##if [[ ${TYPE_OF_PROGRAM} == 'TRANFORM_BEFORE_DEEPREG_COLIHM62_MOVING' ]]; then
##  echo " I AM AT TRANFORM_BEFORE_DEEPREG_COLIHM62_MOVING" >> /software/ERROR.txt
##  echo " I AM AT TRANFORM_BEFORE_DEEPREG_COLIHM62_MOVING"
##  /software/lin_transform_before_deepreg_COLIHM62_as_moving.sh $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST /input1 /output
##fi
##if [[ ${TYPE_OF_PROGRAM} == 'APPLY_MAT_TRANFORM_BEFORE_DEEPREG_COLIHM62_MOVING' ]]; then
##  echo " I AM AT APPLY_MAT_TRANFORM_BEFORE_DEEPREG_COLIHM62_MOVING" >> /software/ERROR.txt
##  echo " I AM AT APPLY_MAT_TRANFORM_BEFORE_DEEPREG_COLIHM62_MOVING"
##  /software/lin_transform_before_deepreg_COLIHM62_as_moving_mat_on_mask.sh $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST /input1 /output
##fi
## STEP 2 : Deepreg application to find the deformation matrix and further align template to the session CT
##if [[ ${TYPE_OF_PROGRAM} ==   'NO_RAPIDS_APPLYDEEPREG_CSF_CISTERN_MIDLINE_SEP_COLIHM62' ]]; then
##  echo " I AM AT TYPE_OF_PROGRAM==NO_RAPIDS_APPLYDEEPREG_CSF_CISTERN_MIDLINE_SEP_COLIHM62"
##  /software/deepregbasedregis_csf_cistern_midline_separation_with_COLIHM62_NO_RAPIDS.sh $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST /input1 /output
##fi
## STEP 3 : Separate compartments based on this alignment
##if [[ ${TYPE_OF_PROGRAM} == "VENT_BOUND_IN_SNIPR_CSF_WITH_CISTERN_MIDLINE_WITH_COLI_HM62" ]] ;
##then
##echo " I AM IN TYPE_OF_PROGRAM == VENT_BOUND_IN_SNIPR_CSF_WITH_CISTERN_MIDLINE_WITH_COLI_HM62"
##/software/compartment_separation_with_vent_afterdeepreg_with_cistern_with_COLIHM62_with_midline.sh  $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST
###/software/compartment_separation_with_vent_afterdeepreg_include_hemorrhagemask.sh  $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST
##fi
## STEP 4 : Produce PDFs and CSVs for these compartment
##if [[ ${TYPE_OF_PROGRAM} == 'PDF_AFTER_CSF_COMPARTMENT_WITH_DEEPREG_FOR_NON_SAH' ]]; then
##  /software/csf_compartments_vols_N_display_for_NON_SAH.sh ${SESSION_ID} $XNAT_USER $XNAT_PASS $XNAT_HOST
##fi