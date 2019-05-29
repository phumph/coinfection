#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
		sbatch raxml_run.sh ${line}
		sleep 2
	done < ${1} # take list of files to process
