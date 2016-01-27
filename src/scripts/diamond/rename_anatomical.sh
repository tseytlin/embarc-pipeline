#!/bin/bash
# rename anatomical slides
ANAT=anat

function rename () {
	SOURCE=$1
	TARGET=$2
	if [ $SOURCE != $TARGET ]; then
		SB=$(basename $SOURCE .nii)
		TB=$(basename $TARGET .nii)
		echo "  $SOURCE -> $TARGET"
		mv $SOURCE tmp/$TARGET
		mv $SB.acpc.csv tmp/$TB.acpc.csv 2> /dev/null
	fi
}

# go over all parameters and find all anat directory
for DIR in $(find $@ -name ${ANAT} -type d)
do
	# get the files of interest
	cd $DIR	
	ORIG=$(ls *_anat.nii*)
	ORIENT=$(ls *_anat_orient.nii*)
	CROP=$(ls *_anat_crop.nii*)

	if [ -z $CROP ]; then
		continue;
	fi

	# remember names
	NM[0]=$ORIG
	NM[1]=$ORIENT
	NM[2]=$CROP

	# get orientation info
	OR[0]=$(fslhd $ORIG | grep  "qform_.orient")
	OR[1]=$(fslhd $ORIENT | grep  "qform_.orient")
	OR[2]=$(fslhd $CROP | grep  "qform_.orient")

	# find the CROP image as it is the smallest file
	SZ[0]=$(stat --printf="%s" $ORIG)
	SZ[1]=$(stat --printf="%s" $ORIENT)
	SZ[2]=$(stat --printf="%s" $CROP)

	# find crop image offset
	NEW_CROP=0
	MIN_SZ=${SZ[0]}
	for i in {1..2}
	do
		if [ ${SZ[i]} -lt $MIN_SZ ]; then
			NEW_CROP=$i
			MIN_SZ=${SZ[i]}
		fi
	done

	# now orient image will have the same orientation as crop
	for i in {0..2} 
	do
		# if we have a different offset then crop, but same orientation
		if [ $i -ne $NEW_CROP ]; then
			if [ "${OR[$NEW_CROP]}" == "${OR[$i]}" ]; then
				NEW_ORIENT=$i
			else
				NEW_ORIG=$i
			fi
		fi
	done

	# now do a rename
	SUBJECT=$(basename $(dirname $DIR))
	echo "processing $SUBJECT .."	
	mkdir tmp	
	
	rename ${NM[$NEW_ORIG]} ${NM[0]} 
	rename ${NM[$NEW_ORIENT]} ${NM[1]}
	rename ${NM[$NEW_CROP]} ${NM[2]}
	
	mv tmp/* . 2> /dev/null
	rmdir tmp
done



