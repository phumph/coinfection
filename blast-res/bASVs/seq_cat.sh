#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
		cat ${line}* > ${line}_all_seqs.fasta
	done < ${1} # take list of files to process
