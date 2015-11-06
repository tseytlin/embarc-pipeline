#!/bin/sh

# go over input files
for DIR in $@
do
	BEFORE=`du -sh ${DIR} | awk '{print $1; }'`
	echo -n "Uncompacting ${DIR} from ${BEFORE} to ... "
	find ${DIR} -name '*.gz' -exec gunzip -f '{}' \; 
	AFTER=`du -sh ${DIR} | awk '{print $1; }'`
	echo "$AFTER"
done
