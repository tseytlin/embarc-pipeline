#!/bin/bash

for line in `cat  ROI-map.csv`
do
	OLD=`echo $line | cut -d, -f1 | cut -d. -f1`
	NEW=`echo $line | cut -d, -f2 `
	NEW=${NEW::-10}
	
	ls $OLD*.[ni]* &> /dev/null
	if [ $? -eq 0 ]; then
		echo "$OLD -> $NEW"
		fslmaths $OLD*.[ni]* -mul 1 ${NEW}_mni.nii
		rm $OLD*.[nhi]*
	fi

done
