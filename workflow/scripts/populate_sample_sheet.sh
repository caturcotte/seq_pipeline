#!/bin/bash

INPUT=$(pwd)
ONT_FOLDERS=$(find "$INPUT" -type d -maxdepth 1)
TOTAL="${i}"
KIT="SQK-NBD114_24"
for ((i = 1; i < TOTAL; ++i))
do
	for j in $ONT_FOLDERS
	do
		SAMPLES=$(find "$i/fastq_pass" -type d -maxdepth 1)
		IFS='_' read -ra ADDR <<< "$j"
		# DATE=${ADDR[0]}
		EXPERIMENT_ID=${ADDR[1]}
		POSITION_ID=${ADDR[2]}
		FLOW_CELL_ID=${ADDR[3]}
		PROTOCOL_RUN_ID=${ADDR[4]}
		for k in $SAMPLES
		do
			printf "WT,%03d,%s,%s,%s,%s,%s,%s\n" "$i" "$EXPERIMENT_ID" "$POSITION_ID" "$FLOW_CELL_ID" "$PROTOCOL_RUN_ID" "$k" "$KIT" >> "$FILE"
		done
	done
done
