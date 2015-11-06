#!/bin/bash
# run DTI sequence
# author: Amelia Versace, Eugene Tseytlin 
# University of Pittsburgh

#
# define some variables
#

BASE=~/Data/dti/
BFS=$BASE/freesurfer
EXPLOREDTI=${BDATA}/exploreDTI

SUBJECT=cinci_20120808.22130
MGZ_FILE="${BFS}/${SUBJECT}/mri/orig/001.mgz"
BVEC=""
BVAL=""
NII=""

# make sure that required mgz file is there
if [ ! -f "${MGZ_FILE}" ]; then
	echo "Error: ${MGZ_FILE} doesn't exist"
	exit 1 
fi

echo "Command: recon-all -subjid ${BFS}/${SUBJECT} -autorecon1 -no-isrunning"
recon-all -subjid ${BFS}/${SUBJECT} -autorecon1 -no-isrunning 

# NOTE4EUGENE.After 'autorecon1' ExploreDTI can start.
# can put those into background
echo "Command: recon-all -subjid ${BFS}/${SUBJECT} -autorecon2 -no-isrunning"
recon-all -subjid ${BFS}/${SUBJECT} -autorecon2 -no-isrunning 

echo "Command: recon-all -subjid ${BFS}/${SUBJECT} -autorecon3 -no-isrunning -qcache"
recon-all -subjid ${BFS}/${SUBJECT} -autorecon3 -ncd ..o-isrunning -qcache


# prepare dataset to run ExploreDTI
OUT=${EXPLOREDTI}/REORDERED/${SUBJECT}.bval
if [ -e ${OUT} ]; then
        echo "Warning: Data for subject ${SUBJECT} has already been done "
	echo "   (if you want to re-run remove ${EXPLOREDTI}/REORDERED/${SUBJECT}.* first"
        exit 1
fi


#
# need a handle to DTI nii,bvec and bvol
#

# copy bvec/bval files to ORIGINAL directory
#echo "Moving nii/bvec/bval files from dicom folder to destination folder ... " 
#echo "Command: mv ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/*.nii.gz ${EXPLOREDTI}/ORIGINAL/LAS." 
#mv ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/*.nii.gz ${EXPLOREDTI}/ORIGINAL/LAS
#echo "Command: mv ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/*.bvec ${EXPLOREDTI}/ORIGINAL/LAS." 
#mv ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/*.bvec ${EXPLOREDTI}/ORIGINAL/LAS 
#echo "Command: mv ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/*.bval ${EXPLOREDTI}/ORIGINAL/LAS." 
#mv ${BSOURCE}/${SUBJECT}/*_dti_68_1024x1024.*/*.bval ${EXPLOREDTI}/ORIGINAL/LAS

#
# Assumed that BVEC/BVAL and NII are already in $EXPLOREDTI/ORIGINAL/LAS
#

# re-orient LAS to RAS
echo "Convert LAS to RAS .."
mriconvert ${EXPLOREDTI}/ORIGINAL/LAS/${SUBJECT}.nii.gz ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}.nii.gz --out_orientation RAS


# applying bet to data 
echo " Brain Extraction Tool... w/o -F bet just te 1st volume==B0 image "
echo " - Command: ${BET} ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}.nii.gz ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}.nii.gz -f 0.3 -g 0 -v -m" 
echo " - Output-----------"
bet ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}.nii.gz ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}_brain.nii.gz -f 0.3 -g 0 -v -m
echo " - -----------------"
echo " - Done"

# apply mask to of first image to the rest of images
echo " - Command: fslmaths ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}.nii.gz -mas ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}_brain_mask.nii.gz ${EXPLOREDTI}/ORIGINAL/${SUBJECT}.nii.gz"
fslmaths ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}.nii.gz -mas ${EXPLOREDTI}/ORIGINAL/2bet/${SUBJECT}_brain_mask.nii.gz ${EXPLOREDTI}/ORIGINAL/${SUBJECT}.nii.gz
gzip -d ${EXPLOREDTI}/ORIGINAL/${SUBJECT}.nii.gz

# reording bvec and bval for some reason
echo " - Command:head ${EXPLOREDTI}/ORIGINAL/${SUBJECT}.bvec |sed "s/ 0 / /g" | sed "s/^0 /0 0 0 0 0 0 0 /g" > ${EXPLOREDTI}/REORDERED/${SUBJECT}.bvec"
echo " - Output-----------"
head ${EXPLOREDTI}/ORIGINAL/${SUBJECT}.bvec |sed 's/ 0 / /g' | sed 's/^0 /0 0 0 0 0 0 0 /g' > ${EXPLOREDTI}/REORDERED/${SUBJECT}.bvec
head ${EXPLOREDTI}/ORIGINAL/${SUBJECT}.bval |sed 's/ 0 / /g' | sed 's/^0 /0 0 0 0 0 0 0 /g' > ${EXPLOREDTI}/REORDERED/${SUBJECT}.bval

# re-orient T1 image (struct) 
mriconvert ${BFS}/${SUBJECT}/mri/brain.mgz ${EXPLOREDTI}/T1_1mm_RAS/${SUBJECT}_T1_masked_RAS.nii --out_orientation RAS

# invoke MATLAB with explore DTI scripts
matlab -r "dti_sequence($EXPLOREDTI,$SUBJECT); exit" 


# NOW TRACULA bit
# this is specially crafted script customized for each subject, need to genereate it
DMRIRC = /data/Phillips2/projects/dtistudy/BIOS/data/tracula/scripts/dmrirc.${SUBEJCT}
trac-all –c $DMRIRC –prep
trac-all –c $DMRIRC –bedp
trac-all –c $DMRIRC –path

# we are DONE!

