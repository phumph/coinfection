#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
		mafft "${line}" > "${line}.aln"
		fasta2phylip.pl "${line}.aln" "${line}.phy"
	done < ${1} # take list of files to process


