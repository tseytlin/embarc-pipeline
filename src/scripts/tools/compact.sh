#!/bin/sh

if [ -z $@ ]; then
	echo "Usage:  compact [-zip] <dir1> <dir2> ... <dirN>"
	echo "\tReduce size of analyzed fMRI dataset by deleting residual images, "
	echo "\tNiPype analysis directory and optionally compressing fMRI images. "
	exit 1
fi

# go over input files
ZIP=0	
for DIR in $@
do
	if [ "$DIR" = "-zip" ]; then
		ZIP = 1
	else	
		BEFORE=`du -sh ${DIR} | awk '{print $1; }'`
		echo -n "Compacting ${DIR} from ${BEFORE} to ... "
		# removing analysis folder for embarc-2.0
		rm -rf $(find ${DIR} -name 'analysis' -type d)
		# remove residual images
		find ${DIR}  -name 'ResI_*' -exec rm -f '{}' \; 
		# archive all images
		if [ $ZIP ]; then		
			find ${DIR}  -name '*.[nimh]?[igtr]' -exec gzip -f '{}' \; 
		fi
		AFTER=`du -sh ${DIR} | awk '{print $1; }'`
		echo "$AFTER"
	fi
done
