#!/bin/bash
#


# 2exploreDTI [part 1]
# run pipeline to exploreDTI [part 1] on the specified subject
#

PROG=`basename $0`

LIB="/data/Phillips2/projects/dtistudy/LAMS/scripts/avLib.sh"
if [ ! -e ${LIB} ]; then
  echo " Library not found!"
  exit 1
fi
. ${LIB}

# Parse Input options
HELP="no";
SUBJECTS="";
while getopts "hs:" OPTION; do
  case "$OPTION" in
    s) SUBJECTS="$OPTARG";;
    h) HELP="yes" ;;
  esac
done

if [ "-${HELP}-" == "-yes-" ]; then
  echo ""
  echo " Script: ${PROG}"
  echo ""
  echo " Options:"
  echo " -h: print this help and exit"
  echo " -s: coma separated list of subjects to be analyzed. It is required. To run on all subject use -s all"
  echo "     Example: if you want to run subjects 10 and 12, this is the command to be used:"
  echo "  		> 2exploreDTI.sh -s 10,12, [no space and comma after the last subject too]"
  echo ""
  exit 0
fi

echo "+++++++++++++++++++++++++++++++++++"
echo "+ ${PROG}"
echo "+ Run 2exploreDTI for the specific subject"
echo "+ Arguments: $*"
echo "+"
echo "+ Start: "`date`
echo "+++++++++++++++++++++++++++++++++++"

# check if we have the subject
echo -n "Checking Subject ${SUBJECTS}... "
if [ "-${SUBJECTS}-" == "--" ]; then
  echo "Error"
  echo " Please Specify Subject to prepare."
  exit 2
elif [ "-${SUBJECTS}-" == "-all-" ]; then
  echo "OK"
  echo " Running in all subjects"
  LIST1=`ls -d /data/Phillips2/projects/dtistudy/LAMS/data/exploreDTI/ORIGINAL/*.nii | sed "s/.nii//" |xargs -n1 basename `
elif [ "-`echo ${SUBJECTS} | grep ,`-" != "--" ]; then
  LIST2=`echo ${SUBJECTS} | sed "s/,/ /g"`
else
  if [ -e ${BSOURCE}/${SUBJECTS} ]; then
    LIST2=${SUBJECTS}
    echo "OK"
  else
    echo "Error"
    echo " Please Specify a valid Subject."
    exit 3
  fi
fi
# run in all the subjects in the list
for SUBJECT in $LIST1 $LIST2; do
OUT=${EXPLOREDTI}/REORDERED/${SUBJECT}.bval
if [ -e ${OUT} ]; then
        echo " -- Data for subject ${SUBJECT} has already been done (if you want to re-run remove ${EXPLOREDTI}/REORDERED/${SUBJECT}.* first"
        echo "    Skip subject ${SUBJECT} "
        continue
fi
  # converts dicom to nii
  echo "Converting DTI dicom images in to nii format... " >> ${LOGS}
  echo "Command: dcm2nii ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/MR* " >> ${LOGS}
  echo "Output-----------" >> ${LOGS}
  echo "Converting DTI dicom images in to nii format... Done" >> ${LOGS}
  ${DCM2NII} ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/MR* 
  echo "-----------------" >> ${LOGS}
  echo "Moving nii/bvec/bval files from dicom folder to destination folder ... Done" >> ${LOGS}
  echo "Command: mv ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/*.nii.gz ${EXPLOREDTI}/ORIGINAL/LAS." 
  mv ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/*.nii.gz ${EXPLOREDTI}/ORIGINAL/LAS.
  echo "Command: mv ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/*.bvec ${EXPLOREDTI}/ORIGINAL/LAS." 
  mv ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/*.bvec ${EXPLOREDTI}/ORIGINAL/LAS.  
  echo "Command: mv ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/*.bval ${EXPLOREDTI}/ORIGINAL/LAS." 
  mv ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/*.bval ${EXPLOREDTI}/ORIGINAL/LAS.
  
# applying bet to data
${MRI_CONVERT} ${EXPLOREDTI}/ORIGINAL/LAS/${SUBJECT}.nii.gz ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}.nii.gz --out_orientation RAS
  echo " Brain Extraction Tool... w/o -F bet just te 1st volume==B0 image "
  echo " - Command: ${BET} ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}.nii.gz ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}.nii.gz -f 0.3 -g 0 -v -m" 
  echo " - Output-----------"
  ${BET} ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}.nii.gz ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}_brain.nii.gz -f 0.3 -g 0 -v -m
  echo " - -----------------"
  echo " - Done"
  echo " - Command:${FSLMATHS} ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}.nii.gz -mas ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}_brain_mask.nii.gz ${EXPLOREDTI}/ORIGINAL/${SUBJECT}.nii.gz"
${FSLMATHS}	${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}.nii.gz -mas ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}_brain_mask.nii.gz ${EXPLOREDTI}/ORIGINAL/${SUBJECT}.nii.gz
gzip -d ${EXPLOREDTI}/ORIGINAL/${SUBJECT}.nii.gz

# NOTE4EUGENE: Reorder ${EXPLOREDTI}/ORIGINAL/${SUBJECT}.nii with  /home/versacea/bin/ExploreDTI_2011/Source/Plugins/E_DTI_reorderB0s_to_beginning.p 

echo " - Command:head ${EXPLOREDTI}/ORIGINAL/${SUBJECT}.bvec |sed "s/ 0 / /g" | sed "s/^0 /0 0 0 0 0 0 0 /g" > ${EXPLOREDTI}/REORDERED/${SUBJECT}.bvec"
echo " - Output-----------"
head ${EXPLOREDTI}/ORIGINAL/${SUBJECT}.bvec |sed 's/ 0 / /g' | sed 's/^0 /0 0 0 0 0 0 0 /g' > ${EXPLOREDTI}/REORDERED/${SUBJECT}.bvec
head ${EXPLOREDTI}/ORIGINAL/${SUBJECT}.bval |sed 's/ 0 / /g' | sed 's/^0 /0 0 0 0 0 0 0 /g' > ${EXPLOREDTI}/REORDERED/${SUBJECT}.bval

${MRI_CONVERT} ${BFS}/${SUBJECT}/mri/brain.mgz ${EXPLOREDTI}/T1_1mm_RAS/${SUBJECT}_T1_masked_RAS.nii --out_orientation RAS
# NOTE4EUGENE: 
#from now on /home/versacea/bin/ExploreDTI/Source/Plugins/

# Resample ${EXPLOREDTI}/T1_1mm_RAS/${SUBJECT}_T1_masked_RAS.nii to [2 2 2] mm^3:
# script: home/versacea/bin/ExploreDTI/Source/Plugins/E_DTI_Resample_nii.p (or ???  E_DTI_Resample_nii_ex.p  or E_DTI_Resample_nii_ex_me.p or E_DTI_resample_nii_file.p  ??? NOT SURE WHICH ONE IS THE CORRECT ONE) 
# input: ${EXPLOREDTI}/T1_1mm_RAS/${SUBJECT}_T1_masked_RAS.nii 
# destination/output: ${EXPLOREDTI}/REORDERED/${SUBJECT}_T1.nii
# voxel size: [2 2 2]


# Convet bvec and bval in B-matrix 
# script: home/versacea/bin/ExploreDTI/Source/Plugins/E_DTI_convert_nii_dic_2_txt.p
# input and output destination folder: ${EXPLOREDTI}/REORDERED/.


# Calculate DTI*.mat file --> convert ${EXPLOREDTI}/REORDERED/${SUBJECT}.nii to ${EXPLOREDTI}/REORDERED/${SUBJECT}.mat file with the follwoing parameters  [permute (y x z); flip (-x y z) B value (1000); voxel size (2 2 2); N no-DWI (7); N DWI (61); matrix size (128 128 64)] 

echo	"done"
done
echo 	"+++++++++++++++++++++++++++++++++++"
echo 	"+ End: "`date`
echo 	"+++++++++++++++++++++++++++++++++++"

echo 	"+++++++++++++++++++++++++++++++++++"
echo 	"WHAT TO DO NEXT "
echo 	"+++++++++++++++++++++++++++++++++++"

