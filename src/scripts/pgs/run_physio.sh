#!/bin/sh 
# Eugene Tseytlin (University of Pittsburgh)
# run physio script

SCRIPT_DIR=`dirname $0`
SCRIPT_DIR=`(cd $SCRIPT_DIR;pwd)`
PROG=`basename $0` 
export MATLABPATH=$SCRIPT_DIR:$MATLABPATH
CWD=`pwd`

# do some parameter checking
if [ -z $1 ];
then
    echo "Usage: ${PROG} <data directory 1> ...[data directory N]"
    echo "Run an PHYSIO analysis on a set of subject data directories for cry study." 
exit 1
fi

# go over each dataset one at a time
for DATA_DIR in $@
do
	
	# get some info
	NAME=`basename $DATA_DIR`
	USERNAME=`echo $NAME | awk -F '.' '{print $2}'`
	DT=`date`
	DATA_DIR_PATH=`(cd $DATA_DIR;pwd)`
	LOG_DIR=${DATA_DIR_PATH}/logs
	
	# make sure that this is a valid data directory
	if [ ! -d ${DATA_DIR}/physio ];
	then
		echo "Error: ${DATA_DIR} does not seem to be a valid CRY subject directory" 2>&1 | tee -a $LOG_FILE
		exit 1;
	fi
	
	# create a log file
	if [ ! -d $LOG_DIR ];
	then
	    mkdir -p $LOG_DIR
	fi
	# transfer.subject.N.log
	# find empty log file
	i=1
	LOG_FILE="${LOG_DIR}/physio.${USERNAME}.${i}.log" 
	while [ -e $LOG_FILE ] 
	do
	   let "i=i+1"
	   LOG_FILE="${LOG_DIR}/physio.${USERNAME}.${i}.log" 
	done

	
	echo "+------------------------------------------------------------" > $LOG_FILE 2>&1
	echo "|  Script:      ${PROG} " >> $LOG_FILE 2>&1
	echo "|  Description: run an PHYSIO analysis on a set of subject data directories for cry study" >> $LOG_FILE 2>&1
	echo "|  Start:       ${DT} " >> $LOG_FILE 2>&1
	echo "|  Argument:    ${DATA_DIR_PATH} " >> $LOG_FILE 2>&1
	echo "|  Subject:     ${USERNAME} " >> $LOG_FILE 2>&1
	echo "+------------------------------------------------------------" >> $LOG_FILE 2>&1
	echo " " >> $LOG_FILE 2>&1
	echo "Running Physio Script on CRY dataset for ${USERNAME}.." 2>&1 | tee -a $LOG_FILE
	
	# go into directory in question
	cd ${DATA_DIR}/physio
	
	# extract required values
	# copy/paste from clear_stats.sh 
	sed '/ECG/,$d' wpc*_${USERNAME}_*.puls  > new_pulseox.txt
	sed '/ECG/,$d' wpc*_${USERNAME}_*.resp  > new_resp.txt
	#
	grep LogStartMPCUTime: wpc*_${USERNAME}_*.puls > start_pulse.txt
	grep LogStartMPCUTime: wpc*_${USERNAME}_*.puls > start_resp.txt
	#
	more start_pulse.txt | cut -d ':' -f 2 > pulse_onset_time.txt 
	more start_resp.txt | cut -d ':' -f 2 > resp_onset_time.txt

	rm start_pulse.txt
	rm start_resp.txt
    # copy/paste from clear_stats.sh
		
	# call an PHYSIO script
	(matlab -nosplash -nodesktop -r "run_physio, quit") 2>&1 | tee -a $LOG_FILE
		
	# make an error summary
	echo "Done!" | tee -a $LOG_FILE
	
	cd $CWD
done
reset
