#!/bin/bash
SESSION_ID=${1}
XNAT_USER=${2}
XNAT_PASS=${3}
TYPE_OF_PROGRAM=${4}
echo TYPE_OF_PROGRAM::${TYPE_OF_PROGRAM}
#'https://redcap.wustl.edu/redcap/api/' #
echo ${REDCAP_API}
#export REDCAP_API=${6}
#echo REDCAP_API::${REDCAP_API}
# The input string
input=$XNAT_HOST ##"one::two::three::four"
# Check if '::' is present
if echo "$input" | grep -q "+"; then
  # Set the delimiter
  IFS='+'

  # Read the split words into an array
  read -ra ADDR <<< "$input"
  export XNAT_HOST=$XNAT_HOST ##{ADDR[0]}
  SUBTYPE_OF_PROGRAM=${ADDR[1]} 
else
export XNAT_HOST=$XNAT_HOST ##{5}
    echo "'+' is not present in the string"
fi

echo ${TYPE_OF_PROGRAM}::TYPE_OF_PROGRAM::${SUBTYPE_OF_PROGRAM}::${ADDR[0]}::${ADDR[2]}::${ADDR[3]}
######################### MIDLINE SHIFT 3D PROFILE PIPELINE##############################################################
if [[ ${TYPE_OF_PROGRAM} == 'CSF_COMPARTMENTS' ]]; then
  echo ${SESSION_ID}
  /software/STEPS_COMPARTMENTS_CSF.sh $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST /input1 /output
fi
if [[ ${TYPE_OF_PROGRAM} == 'CSF_COMPARTMENTS_BATCH' ]]; then
  echo ${SESSION_ID}
  PROJECT_ID=${SESSION_ID}
  /software/CSF_COMPARTMENT_IN_BATCH.sh ${PROJECT_ID} $XNAT_USER $XNAT_PASS $XNAT_HOST /input1 /output
fi
if [[ ${TYPE_OF_PROGRAM} == 'PROJECT_DEEPREG_FOR_CSF_COMPARTMENTS' ]]; then
  echo ${SESSION_ID}
/software/CSF_COMPARTMENT_DEEPREG_IN_BATCH.sh $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST /input1 /output
fi
if [[ ${TYPE_OF_PROGRAM} == 'SESSION_DEEPREG_FOR_CSF_COMPARTMENTS' ]]; then
  echo ${SESSION_ID}
/software/deepregbasedregis_csf_cistern_midline_separation_with_COLIHM62_NO_RAPIDS.sh $SESSION_ID $XNAT_USER $XNAT_PASS $XNAT_HOST /input1 /output
fi



