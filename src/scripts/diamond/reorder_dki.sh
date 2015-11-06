#!/bin/bash
# author: Eugene Tseytlin (University of Pittsburgh)
#
export FSLOUTPUTTYPE=NIFTI

if [ -z $1 ]; then
	echo "Usage: reorder_dki.sh <original folder> <reordered folder> <subject id>"
	exit 1
fi

ORIGINAL=$(cd $1; pwd)
REORDERD=$(cd $2; pwd)
SUBJECT=$3

cd $ORIGINAL
NII_FILE="${SUBJECT}.nii"
BVEC_FILE="${SUBJECT}.bvec"
BVAL_FILE="${SUBJECT}.bval"

# split nifti file
fslsplit $NII_FILE vol -t

# save files in arrays
string=$(cat $BVAL_FILE)
BVAL_ARRAY=(${string// / })
NII_ARRAY=(vol*.nii)

# create a list of items as a spreadsheet
rm list >& /dev/null
for i in "${!BVAL_ARRAY[@]}"
do
	BVAL=${BVAL_ARRAY[$i]}
	NII=${NII_ARRAY[$i]}
	j=$(($i+1))
	BVEC=$(awk -v col=$j '{print $col}' $BVEC_FILE | tr  '\n' '|')
	echo "$BVAL|$NII|$BVEC" >> list
done

# sort list 
sort -n list > sorted.list

# merge files back in different order
FILES=$(awk -F '|' '{print $2}' sorted.list | tr '\n' ' ')
fslmerge -t $REORDERD/$NII_FILE $FILES

# reorganize bval and bvec files
awk -F '|' '{print $1}' sorted.list | tr '\n' ' ' > $REORDERD/$BVAL_FILE
awk -F '|' '{print $3}' sorted.list | tr '\n' ' ' > $REORDERD/$BVEC_FILE
echo "" >> $REORDERD/$BVEC_FILE
awk -F '|' '{print $4}' sorted.list | tr '\n' ' ' >> $REORDERD/$BVEC_FILE
echo "" >> $REORDERD/$BVEC_FILE
awk -F '|' '{print $5}' sorted.list | tr '\n' ' ' >> $REORDERD/$BVEC_FILE
echo "" >> $REORDERD/$BVEC_FILE

# cleanup
rm vol*.nii list sorted.list
