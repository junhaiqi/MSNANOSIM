#!/bin/bash
num_of_args=$#

if [ $num_of_args -eq 2 ]; then

    python module/simulateNanoSigs/generatNoiseSignal.py -i $1 -o $2 -rootName timeSeries -seed 0
        
else
  echo "Usage: bash $0 *.fasta output_folder."
  echo "The first parameter indicates the fasta file containing the DNA sequence, and the second parameter indicates the folder path of the output simulated nanopore signal."
fi
 