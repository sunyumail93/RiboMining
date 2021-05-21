#!/usr/bin/env bash
#This is a pipeline to get bed12 transcript length individually.
#Input: just bed12 file
#Output: two columns: TranscriptName|Length(nt)
#Version: Yu Sun, 2018-03-09
#Version: Yu Sun, 2019-08-06
#Version: Yu Sun, 2019-08-08, Fix bug that caused by col11, some records have "," at the end, but some don't

if [ ! -n "$1" ]
then
  echo "    This is a pipeline to get bed12 transcript length individually."
  echo "    Input: just bed12 file"
  echo "    Usage: `basename $0` [bed12]"
  echo "    Output: bed12.len.txt"
else
  Data=$1
  awk '{OFS="\t";split($11,a,",");sum=0;if ($11~/,$/) {for (i=1;i<length(a);i++) sum=sum+a[i];print $4,sum}else {for (i=1;i<=length(a);i++) sum=sum+a[i];print $4,sum}}' $Data > ${Data}.len.txt
fi
