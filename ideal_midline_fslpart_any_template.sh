##!/usr/bin/env bash
#template_file_name=${2} #scct_strippedResampled1.nii.gz
#mask_on_template=${3} #midlinecssfResampled1.nii.gz
#input_filename=${1}
#inv_transformmatrix_file=${4} #{output_filename}_${exten}lin1Inv.mat
#/usr/local/fsl/bin/convert_xfm -omat $inv_transformmatrix_file -inverse $transformmatrix_file
#transformed_output_file=$template_mask_output_directory/${mask_on_template%.nii*}$file_basename".nii.gz"
#echo $transformmatrix_file
#echo $inv_transformmatrix_file
#/usr/local/fsl/bin/flirt -in ${midline_ontemplate_nii} -ref $x -out ${transformed_output_file} -init ${inv_transformmatrix_file} -applyxfm
#echo " REGISTERED OUTPUT FILE "
#echo $transformed_output_file
##
##output_directory=$(dirname ${input_filename})
##echo "I am working in Docker"
##echo $input_filename
##result_output_directory=${output_directory} #/output
##output_for_BET=${output_directory} #/output #/BET_OUTPUT
##nifti_reg_output_directory=${output_directory} #/output #/LINEAR_REGISTRATION_OUTPUT
##template_directory=/templatenifti
##input_directory=$output_for_BET
##output_directory=$nifti_reg_output_directory
### template_file_name=$templatefilename
##extension="_brain_f.nii.gz"
##input_filename_dirname="$(dirname -- $input_filename)"
##input_filename_basename="$(basename -- $input_filename)"
##x=$input_filename_dirname/${input_filename_basename%.nii*}$extension
##template_mask_directory=/templatemasks
##output_directory=$nifti_reg_output_directory
##template_mask_output_directory=${output_directory}  #/output #/TEMPLATEMASKS_OUTPUT_LR
##filename=$x
##echo $filename
##template_image=$template_directory/$template_file_name
##template_basename="$(basename -- $template_image)"
##exten=${template_basename%.nii*}
##echo $template_basename
##img=${filename%.nii*}
##img_basename=$(basename -- $img)
##output_filename=$output_directory/$img_basename
##echo $output_filename
### mask_on_template=$mask_on_template
##midline_ontemplate_nii=$template_mask_directory"/"$mask_on_template
##file_basename=$(basename $x)
##file_basename=${file_basename%.nii*}
##filename_without_ext=${x%.nii*}
##transformmatrix_file=${output_filename}_${exten}lin1.mat
##inv_transformmatrix_file=${output_filename}_${exten}lin1Inv.mat
##/usr/local/fsl/bin/convert_xfm -omat $inv_transformmatrix_file -inverse $transformmatrix_file
##transformed_output_file=$template_mask_output_directory/${mask_on_template%.nii*}$file_basename".nii.gz"
##echo $transformmatrix_file
##echo $inv_transformmatrix_file
##/usr/local/fsl/bin/flirt -in ${midline_ontemplate_nii} -ref $x -out ${transformed_output_file} -init ${inv_transformmatrix_file} -applyxfm
##echo " REGISTERED OUTPUT FILE "
##echo $transformed_output_file