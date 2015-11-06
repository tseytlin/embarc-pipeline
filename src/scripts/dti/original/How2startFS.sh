#!/bin/bash
#

# data folders
BDATA=${BASE}"/data"
BSOURCE=${BASE}"/rawdata"
BFS=${BDATA}"/freesurfer"

# Original Data Location. ### NOTE4EUGENE. I WOULD LIKE TO USE THIS FOLDER INSTEAD OF ${BSOURCE} TO AVOID MULITPLE COPIES OF THE DICOMs... NEED TO TALK TO LISA AND RICKI FOR THIS####
ODL="/data/Phillips1/rawdata/BIOS"

LIB="/data/Phillips2/projects/dtistudy/BIOS/scripts/avLib.sh"
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
    u) NSUBJ="$OPTARG";;
  esac
done

if [ "-${HELP}-" == "-yes-" ]; then
  echo ""
  echo " Script: ${PROG}"
  echo "Run pipeline to generate surface reconstruction on the specified subject"
  echo ""
  echo " Options:"
  echo " -h: print this help and exit"
  echo " -s: coma separated LIST1 of subjects to be analyzed. It is required. To run on all subject use -s all"
  echo "     Example: if you want to run subjects 10 and 12, this is the command to be used:"
  echo "  		> How2startFS.sh -s 10,12, or How2startFS.sh -s 10,"
  echo ""
  exit 0
fi

echo "+++++++++++++++++++++++++++++++++++"
echo "+ ${PROG}"
echo "+ Run surface reconstruction for the specific subject"
echo "+ Arguments: $*"
echo "+"
Start=`date`
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
  LIST1=`ls -d /data/Phillips2/projects/BIOS/${SUBJECTS}/*_mprage*_/MR* | xargs -n1 basename `
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
OUT=${BFS}/${SUBJECT}/mri
if [ -e ${OUT} ]; then
        echo " -- Data for subject ${SUBJECT} has already been done "
        echo "    Skip subject ${SUBJECT} "
        continue-autorecon1
fi

mkdir -p ${BFS}/${SUBJECT}/mri/orig
DICOM1=`ls ${BSOURCE}/${SUBJECT}/ANATOMY/MR* |head -n 1`
# NOTE4EUGENE. The 'MR*' should work for any local study; however, make a note that it might not for multicenter studies that do not have a Siemens scanner!!
mri_convert -it dicom -ot mgz ${DICOM1} ${BFS}/${SUBJECT}/mri/orig/001.mgz 
echo "Command: recon-all -subjid ${BFS}/${SUBJECT} -autorecon1 -no-isrunning"
recon-all -subjid ${BFS}/${SUBJECT} -autorecon1 -no-isrunning 
done
# NOTE4EUGENE.After 'autorecon1' ExploreDTI can start.
echo "Command: recon-all -subjid ${BFS}/${SUBJECT} -autorecon2 -no-isrunning"
recon-all -subjid ${BFS}/${SUBJECT} -autorecon2 -no-isrunning 
done
echo "Command: recon-all -subjid ${BFS}/${SUBJECT} -autorecon3 -no-isrunning -qcache"
recon-all -subjid ${BFS}/${SUBJECT} -autorecon3 -no-isrunning -qcache
done

echo " ...Done"

  echo "+++++++++++++++++++++++++++++++++++"
  echo "+ Time: "`date`
  echo "+++++++++++++++++++++++++++++++++++"

