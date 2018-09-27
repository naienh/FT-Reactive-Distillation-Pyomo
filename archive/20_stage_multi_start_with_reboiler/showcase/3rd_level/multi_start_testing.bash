#!/bin/bash
rand_file_name=$1
low=$2
high=$3
output=$4

# loop from low to high
for ((i=$low; i<=$high; i++)); do
       callcmd="./multi_start_onestep.py"
       # generate a testing command
       callcmd="$callcmd $rand_file_name $i $output"
       echo $callcmd
       $callcmd
done
