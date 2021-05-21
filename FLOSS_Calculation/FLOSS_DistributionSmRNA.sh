#!/usr/bin/env bash
#FLOSS_DistributionSmRNA.sh
#This script calculates the individual SmRNA length distribution of the input data on a set of transcripts
#Not like Floss_Reference.sh, this script will output FLOSS counts and freq for individual transcript
#Version: Yu Sun, 2018-11-06

if [ ! -n "$1" ]
then
  echo "    This script calculates the individual smRNA length distribution of the input data on a set of transcripts"
  echo "    Not like Floss_Reference.sh, this script will output FLOSS counts and freq for individual transcript"
  echo "    Usage: `basename $0` [Data.bed13] [Reference.bed12] [OutputPrefix] [Offset|default 0]"
  echo "    Output: OutputPrefix.summary.count and OutputPrefix.summary.freq"
else
  Data=$1
  Annotation=$2
  OutputPrefix=$3
  Seed=$RANDOM$RANDOM
  if [ ! -n "$4" ];then POffset=0;else POffset=$4;fi

  rm -rf ${OutputPrefix}.summary.count
  rm -rf ${OutputPrefix}.summary.freq
  echo "1. Starting Floss individual length distribution estimation"
  echo "  Data: "$Data
  echo "  Annotation: "$Annotation
  echo "  OutputPrefix: "$OutputPrefix
  echo "  P offset: "$POffset

  echo "2. Looping through each transcript..."
  awk '{print $4}' ${Annotation} > ${Annotation}.${Seed}.list
  for i in {26..32};do echo $i;done > ${Annotation}.${Seed}.ref
  echo -e "Name\t26\t27\t28\t29\t30\t31\t32" > ${OutputPrefix}.summary.count
  echo -e "Name\t26\t27\t28\t29\t30\t31\t32" > ${OutputPrefix}.summary.freq
  for name in `cat ${Annotation}.${Seed}.list`
  do
    echo $name
    awk -v t=$name '{OFS="\t";if ($4==t) print $0}' ${Annotation} > ${OutputPrefix}.${name}.bed12
    bedtools intersect -wa -split -b ${OutputPrefix}.${name}.bed12 -a ${Data} > ${OutputPrefix}.${name}.Intersect
    awk '{print length($13)}' ${OutputPrefix}.${name}.Intersect|sort|uniq -c|awk '{OFS="\t"}{print $2,$1}'|awk '$1>25 && $1<33' > ${OutputPrefix}.${name}.Intersect.temp
    awk '{print $1}' ${OutputPrefix}.${name}.Intersect.temp > ${OutputPrefix}.${name}.Intersect.temp.col1
    cat ${Annotation}.${Seed}.ref ${OutputPrefix}.${name}.Intersect.temp.col1 |sort|uniq -c|awk '$1==1'|awk '{OFS="\t";print $2,0}' > ${OutputPrefix}.${name}.Intersect.fill
    cat ${OutputPrefix}.${name}.Intersect.temp ${OutputPrefix}.${name}.Intersect.fill |sort -k1,1 -n > ${OutputPrefix}.${name}.Intersect.count
    SumNum=`awk '{SUM+=$2}END{printf "%.0f\n",SUM}' ${OutputPrefix}.${name}.Intersect.count`
    if [ $SumNum == "0"  ];then
      cp ${OutputPrefix}.${name}.Intersect.count ${OutputPrefix}.${name}.Intersect.freq
    else
      awk -v SumNum=$SumNum '{OFS="\t";print $1,$2/SumNum*100}' ${OutputPrefix}.${name}.Intersect.count > ${OutputPrefix}.${name}.Intersect.freq
    fi

    #Write to the master table:
    awk '{print $2}' ${OutputPrefix}.${name}.Intersect.count|tr '\n' ' '|awk -v name=$name '{OFS="\t";print name,$1,$2,$3,$4,$5,$6,$7}' >> ${OutputPrefix}.summary.count
    awk '{print $2}' ${OutputPrefix}.${name}.Intersect.freq|tr '\n' ' '|awk -v name=$name '{OFS="\t";print name,$1,$2,$3,$4,$5,$6,$7}' >> ${OutputPrefix}.summary.freq
    #Cleanning up
    rm -rf ${OutputPrefix}.${name}.Intersect.temp ${OutputPrefix}.${name}.Intersect.temp.col1 ${OutputPrefix}.${name}.Intersect.fill
    rm -rf ${OutputPrefix}.${name}.bed12 ${OutputPrefix}.${name}.Intersect ${OutputPrefix}.${name}.Intersect.count ${OutputPrefix}.${name}.Intersect.freq ${OutputPrefix}.${name}.Intersect
  done
  rm -rf ${Data}.TempProIndi${Seed}.Shifted ${Data}.TempProIndi${Seed}.Shifted.5end ${Annotation}.${Seed}.list
  rm -rf ${Annotation}.${Seed}.ref
  echo "5. Done"
fi
