echo " I AM WORKING"
SESSION_ID=${1}
function clean_directories(){
rm -r /workinginput/*
rm -r /workingoutput/*
rm -r /input1/*
rm -r /output/*
rm -r /working/*
rm -r /ZIPFILEDIR/*
rm -r /NIFTIFILEDIR/*
rm -r /DICOMFILEDIR/*
rm -r /maskonly/*
rm -r /software/data.h5

}
## STEPS 1: Registration of template ct (COLIHM62) to a session CT.
/software/lin_transform_before_deepreg_COLIHM62_as_moving.sh $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST /input1 /output
clean_directories
/software/lin_transform_before_deepreg_COLIHM62_as_moving_mat_on_mask.sh $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST /input1 /output
clean_directories
/software/deepregbasedregis_csf_cistern_midline_separation_with_COLIHM62_NO_RAPIDS.sh $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST /input1 /output
#clean_directories
#/software/compartment_separation_with_vent_afterdeepreg_with_cistern_with_COLIHM62_with_midline.sh  $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST
#clean_directories
#/software/csf_compartments_vols_N_display_for_NON_SAH.sh ${SESSION_ID} $XNAT_USER $XNAT_PASS $XNAT_HOST
