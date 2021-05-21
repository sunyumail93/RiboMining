#!/usr/bin/env bash
#GetRPF5end.sh
#This scripts gets the 5end 1nt of RPF reads, especially for shifted reads
#Version: Yu Sun, 

if [ ! -n "$1" ]
then
  echo "    This scripts gets the 5end 1nt of RPF reads, especially for shifted reads (after running Shift5endOffset_BED13RPF.py)"
  echo "    Usage: `basename $0` Data.RPF.bed13"
  echo "    Output: Data.RPF.bed13.5end"
else
  awk 'BEGIN{OFS="\t"}{if ($6=="+") print $1,$2,$2+1,$4,$5,$6,$2,$2+1,$9,"1","1","0",$13;else print $1,$3-1,$3,$4,$5,$6,$3-1,$3,$9,"1","1","0",$13}' $1 > ${1}.5end
fi
