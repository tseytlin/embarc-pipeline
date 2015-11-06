#
# library for DTI_68 analysis scripts
#
# Author: Max Novelli <man8@pitt.edu>
#
# Witten for Amlia Versace
#

BASE="/data/Phillips2/projects/dtistudy/BIOS"

# scripts folder
BINDIR=${BASE}"/scripts"

# scripts
# superseeded by avDCM2NII_Bet ### I don't use none of these TBSS scripts anymore####
avPrepSubject=${BINDIR}"/avPrepSubject.sh"
avDCM2NII_Bet=${BINDIR}"/avDCM2NII_Bet.sh"
avDCM2NII_Bet_All=${BINDIR}"/avDCM2NII_Bet_All.sh"
avPrep4dti=${BINDIR}"/avPrep4dti.sh"
avFirstLevel=${BINDIR}"/avFirstLevel.sh"
avL1Prep=${BINDIR}"/avL1Prep.sh"
avBedpostS="${BINDIR}/avBedpostS.sh"
avTbssReg="${BINDIR}/avTbssReg.sh"
avFAiP="${BINDIR}/avFAiP.sh"
avVbmBet="${BINDIR}/avVbmBet.sh"
avVbmSeg="${BINDIR}/avVbmSeg.sh"
avVbmPA="${BINDIR}/avVbmPA.sh"
avVbmRand="${BINDIR}/avVbmRand.sh"

# data folders
BDATA=${BASE}"/data"
BSOURCE=${BASE}"/rawdata"
BPROC=${BDATA}"/processing"
BTBSS=${BDATA}"/tbss"
BFS=${BDATA}"/freesurfer"

# Original Data Location
ODL="/data/Phillips1/rawdata/BIOS"


# Logs folder
BLOGS=${BASE}"/logs"

# Work folder
BWORK=${BASE}"/work"

# temp folder
BTEMP=${BASE}"/temp"

# Reference Brain for filrt step 2 in avL1Prep registration
# Fsl 3 reference
#L1RefBrain="/usr/local/pkg/fsl/etc/standard/avg152T1_brain"
# Fsl 4 reference
L1RefBrain="/usr/local/pkg/fsl-4.1.8-x64/data/standard/MNI152_T1_1mm_brain"

# target for tbss step 2
L2Tbss2Target=""

STANDARDS="${FSLDIR}/etc/standard"
# reference image for vbm segmentation
VbmSegRef_1mm_g=${FSLDIR}/etc/standard/tissuepriors/avg152T1_1mm_gray
VbmSegRef_2mm_g=${FSLDIR}/etc/standard/tissuepriors/avg152T1_gray

# Application REpository
MYAPPDIR="/home/versacea/bin"
APPDIR="/usr/local/pkg"
# FSL folders
FSLDIR="${APPDIR}/fsl-4.1.8-x64"
FSLBIN="${FSLDIR}/bin"
# FSL applications
DTIFIT="${FSLBIN}/dtifit"
BEDPOSTX="${FSLBIN}/bedpostx"
FLIRT="${FSLBIN}/flirt"
CONVERT_XFM="${FSLBIN}/convert_xfm"
TBSS_1_PREPROC="${FSLBIN}/tbss_1_preproc"
TBSS_2_REG="${FSLBIN}/tbss_2_reg"
TBSS_REG="${FSLBIN}/tbss_reg"
TBSS_3_POSTREG="${FSLBIN}/tbss_3_postreg"
TBSS_4_PRESTAT="${FSLBIN}/tbss_4_prestat"
EDDY_CORRECT="${FSLBIN}/eddy_correct"
BET="${FSLBIN}/bet"
FSLROI="${FSLBIN}/fslroi"
FSLSTATS="${FSLBIN}/fslstats"
FSLMATHS="${FSLBIN}/fslmaths"
FSLHD="${FSLBIN}/fslhd"
FSLCREATEHD="${FSLBIN}/fslcreatehd"
CLUSTER="${FSLBIN}/cluster"
ATLASQUERY="${FSLBIN}/atlasquery"
FSLMEANTS="${FSLBIN}/fslmeants"
IMTEST="${FSLBIN}/imtest"
FSLVAL="${FSLBIN}/fslval"
IMCP="${FSLBIN}/imcp"
FSLSLICE="${FSLBIN}/fslslice"
XFIBRES="${FSLBIN}/xfibres"
ZEROPAD="${FSLBIN}/zeropad"
IMGLOB="${FSLBIN}/imglob"
FSLMERGE="${FSLBIN}/fslmerge"
MAKE_DYADIC_VECTORS="${FSLBIN}/make_dyadic_vectors"
RANDOMISE="${FSLBIN}/randomise"
PROBTRACKX="${FSLBIN}/probtrackx"
IMCP="${FSLBIN}/imcp"
JACOBIAN="${FSLBIN}/DR/jacobian"
REMOVE_EXT="${FSLBIN}/remove_ext"
IMREF="${FSLBIN}/imref"
FAST="${FSLBIN}/fast"
FSLSWAPDIM="${FSLBIN}/fslswapdim"
STD_SPACE_ROI="${FSLBIN}/standard_space_roi"
SLICESDIR="${FSLBIN}/slicesdir"
AREG="${FSLBIN}/DR/areg"
NREG="${FSLBIN}/DR/nreg"
TRANSFORMATION="${FSLBIN}/DR/transformation"
AFFINECCPAR="${FSLBIN}/DR/affineCC.par"
DOF2GIPL="${FSLBIN}/DR/dof2gipl"
DOF2WARP="${FSLBIN}/dof2warp"
FINDTHEBIGGEST="${FSLBIN}/find_the_biggest"

# freesurfer folders
ALTAPPDIR="/usr/local"
FREESURFER_HOME="${ALTAPPDIR}/freesurfer"
FREEBIN="${FREESURFER_HOME}/bin"
# freesurfer applications
MRI_CONVERT="${FREEBIN}/mri_convert"
MRIS_CONVERT="${FREEBIN}/mris_convert"
MRI_BINARIZE="${FREEBIN}/mri_binarize"
TKREGISTER2="${FREEBIN}/tkregister2"
MRIS_SAMPLE_PARC="${FREEBIN}/mris_sample_parc"
MRI_ANNOTATION2LABEL="${FREEBIN}/mri_annotation2label"
MRI_MERGELABELS="${FREEBIN}/mri_mergelabels"
MRI_LABEL2VOL="${FREEBIN}/mri_label2vol"
MRI_LABEL_VOLUME="${FREEBIN}/mri_label_volume"
MRI_MERGELABELS="${FREEBIN}/mri_mergelabels"


# exploreDTI folder
EXPLOREDTI=${BDATA}"/exploreDTI"

# dtk folders
DTKBIN="/usr/local/dtk"
HARDI_MAT="${DTKBIN}/hardi_mat"
ODF_RECON="${DTKBIN}/odf_recon"
DTI_RECON="${DTKBIN}/dti_recon"
ODF_TRACKER="${DTKBIN}/odf_tracker"
DTI_TRACKER="${DTKBIN}/dti_tracker"
SPINE_FILTER="${DTKBIN}/spline_filter"
TRACK_TRANSFORM="${DTKBIN}/track_transform"

# other applications
DCM2NII=${MYAPPDIR}"/mricron/dcm2nii"

# atlases in fsl 4
ATLASES[1]="Harvard-Oxford Cortical Structural Atlas"
ATLASES[2]="Harvard-Oxford Subcortical Structural Atlas"
ATLASES[3]="JHU ICBM-DTI-81 White-Matter Labels"
ATLASES[4]="JHU White-Matter Tractography Atlas"
ATLASES[5]="Juelich Histological Atlas"
ATLASES[6]="MNI Structural Atlas"
ATLASES[7]="Oxford Thalamic Connectivity Probability Atlas"
ATLASES[8]="Talairach Daemon Labels"
ATLASES_MAX=8
ATLASES_NUM="1 2 3 4 5 6 7 8"


# useful functions
#
# return 1 if the path is absolute
#
function IsAbsolutePath {
  dir=$1
  if [ "-`echo $dir | grep -E "^/"`-" != "--" ]; then
    # is absolute
    echo 1
  else
    echo 0
  fi
}

#
# Make absolut path
#  - original in fsl 4.
function AbsolutePath {
  file=$1
  def=$2
  path=$3
  #echo "AbsolutePath File =>${file}<="
  #echo "AbsolutePath Def  =>${def}<="
  #echo "AbsolutePath Path =>${path}<="
  # prepare directory name
  if [ "-${file}-" == "--" ]; then
    #echo "AbsolutePath File empty. Using Default"
    file="${def}"
  elif [ `IsAbsolutePath ${file}` -eq 0 ]; then
    #echo "AbsolutePath File is not absolute."
    if [ "-${path}-" != "--" ] && [ -d ${path} ]; then
      #echo "AbsolutePath Using Path"
      file="${path}/${file}"
    else
      #echo "AbsolutePath Using current folder"
      file="`pwd`/${file}"
    fi
  fi
  #echo "AbsolutePath File =>${file}<="
  # check if it is a valid directory
  if [ -e ${file} ]; then
    #echo "AbsolutePath Valid folder"
    ret=${file}
  else
    #echo "AbsolutePath Invalid Folder"
    if [ "-${def}-" != "--" ]; then
      ret=${def}
    else
      ret=${dir}
    fi
  fi
  # eliminates double slash
  ret=`echo ${ret} | sed "s#//#/#g"`
  #echo "AbsolutePath Result =>${ret}<="
  echo ${ret}
}


#
# echo time stamp
#
function PrintTimestamp {
  text=$1
  if [ "-$text-" == "--" ]; then
    text="Time"
  fi
  echo "+++++++++++++++++++++++++++++++++++"
  echo "+ ${text}:"`date`
  echo "+++++++++++++++++++++++++++++++++++"
}

#
# print and execute command
#
function Command {
  cmd=$1
  echo " - Command: ${cmd}"
  echo " - Output---------"
  $cmd
  res=$?
  echo " - ---------------"
  echo " - Result: $res"
  echo " - ---------------"
}

#
# print and execute in the background
#
function CommandBg {
  cmd=$1
  logs=$2
  echo " - Command: ${cmd}"
  echo " - Logs: ${logs}"
  if [ "-${logs}-" == "--" ]; then
    $cmd &
  else
    $cmd 1>${logs} 2>&1 &
  fi
  res=$?
  echo " - Result : ${res}"
}

