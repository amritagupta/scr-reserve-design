#!/bin/bash

JOBS_DIR=$1
JOBS_FILE=$2
PBS_FILE=$3

re='^[0-9]+$'

mkdir -p JOBS/${JOBS_DIR}

while read CMD_STR; do
  echo ${CMD_STR}
  #CMD_STR="./${PROGRAM} --instance=${p}.mps --exp=${EXP_NAME} ${argsStr} > /dev/null"
  CMD_STR=$(echo -n ${CMD_STR} | tr -d "\r")
  
  cat ${PBS_FILE} > JOBS/${JOBS_DIR}/tempjob
  echo ${CMD_STR} >> JOBS/${JOBS_DIR}/tempjob
  
  SUB=""
  while ! [[ ${SUB} =~ ${re} ]]
  do
	#SUB="ERROR: cannot migrate job 'Moab.108133' to PBS - qsub: submit error (No default queue specified MSG=requested queue not found)"
	#SUB=$'\n12345'
	SUB=$(qsub JOBS/${JOBS_DIR}/tempjob)
	if grep -q "error" <<< ${SUB};
	then
		continue
	fi
	SUB=$(grep -o "[0-9]*" <<< ${SUB})
  done
  
  echo "${CMD_STR},${SUB}" >> JOBS/${JOBS_DIR}/job_id.csv
done < ${JOBS_FILE}