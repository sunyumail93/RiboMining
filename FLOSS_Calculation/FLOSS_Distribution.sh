#!/usr/bin/env bash
#FLOSS_Distribution.sh
#This script calculates the individual RPF length distribution of the input data on a set of transcripts
#Not like FLOSS_Reference.sh, this script will output FLOSS counts and freq for individual transcript
#Version: Yu Sun, 2018-11-06
#Version: Yu Sun, 2019-06-20, fix a main bug when counting distribution: add awk '$1>25 && $1<33' when generating temp file. This is not very necessary because RPF data already filter out <26, >32 reads
#                         This is just for when you run small RNA data using this script, previous bug will show because small RNA contains <26 reads.

if [ ! -n "$1" ]
then
  echo "    This script calculates the individual RPF length distribution of the input data on a set of transcripts"
  echo "    Not like FLOSS_Reference.sh, this script will output FLOSS counts and freq for individual transcript"
  echo "    This script will only look at [26,32], and calculate read counts/freq based on that length region"
  echo "    For small RNA data, please run the other script: FLOSS_DistributionSmRNA.sh"
  echo "    Usage: `basename $0` [Data.bed13] [Reference.bed12] [OutputPrefix] [Offset|default 12]"
  echo "    Output: OutputPrefix.summary.count and OutputPrefix.summary.freq"
else
  Data=$1
  Annotation=$2
  OutputPrefix=$3
  Seed=$RANDOM$RANDOM
  if [ ! -n "$4" ];then POffset=12;else POffset=$4;fi

  rm -rf ${OutputPrefix}.summary.count
  rm -rf ${OutputPrefix}.summary.freq
  echo "1. Starting FLOSS individual length distribution estimation"
  echo "  Data: "$Data
  echo "  Annotation: "$Annotation
  echo "  OutputPrefix: "$OutputPrefix
  echo "  P offset: "$POffset

  echo "2. Shifting RPF 5end:"
  Shift5endOffset_BED13RPF.py $Data $POffset ${Data}.TempProIndi${Seed}.Shifted

  echo "3. Get 5end"
  awk 'BEGIN{OFS="\t"}{if ($6=="+") print $1,$2,$2+1,$4,$5,$6,$2,$2+1,$9,"1","1","0",$13;else print $1,$3-1,$3,$4,$5,$6,$3-1,$3,$9,"1","1","0",$13}' ${Data}.TempProIndi${Seed}.Shifted > ${Data}.TempProIndi${Seed}.Shifted.5end

  echo "4. Looping through each transcript..."
  awk '{print $4}' ${Annotation} > ${Annotation}.${Seed}.list
  for i in {26..32};do echo $i;done > ${Annotation}.${Seed}.ref
  echo -e "Name\t26\t27\t28\t29\t30\t31\t32" > ${OutputPrefix}.summary.count
  echo -e "Name\t26\t27\t28\t29\t30\t31\t32" > ${OutputPrefix}.summary.freq
  for name in `cat ${Annotation}.${Seed}.list`
  do
    echo $name
    awk -v t=$name '{OFS="\t";if ($4==t) print $0}' ${Annotation} > ${OutputPrefix}.${name}.bed12
    bedtools intersect -u -wa -split -b ${OutputPrefix}.${name}.bed12 -a ${Data}.TempProIndi${Seed}.Shifted.5end > ${OutputPrefix}.${name}.Intersect
    awk '{for (i=1;i<=$4;i++) print $0}' ${OutputPrefix}.${name}.Intersect > ${OutputPrefix}.${name}.Intersect.bed1
    awk '{print length($13)}' ${OutputPrefix}.${name}.Intersect.bed1|sort|uniq -c|awk '{OFS="\t"}{print $2,$1}'|awk '$1>25 && $1<33' > ${OutputPrefix}.${name}.Intersect.bed1.temp
    awk '{print $1}' ${OutputPrefix}.${name}.Intersect.bed1.temp > ${OutputPrefix}.${name}.Intersect.bed1.temp.col1
    cat ${Annotation}.${Seed}.ref ${OutputPrefix}.${name}.Intersect.bed1.temp.col1 |sort|uniq -c|awk '$1==1'|awk '{OFS="\t";print $2,0}' > ${OutputPrefix}.${name}.Intersect.bed1.fill
    cat ${OutputPrefix}.${name}.Intersect.bed1.temp ${OutputPrefix}.${name}.Intersect.bed1.fill |sort -k1,1 -n > ${OutputPrefix}.${name}.Intersect.bed1.count
    SumNum=`awk '{SUM+=$2}END{printf "%.0f\n",SUM}' ${OutputPrefix}.${name}.Intersect.bed1.count`
    if [ $SumNum == "0"  ];then
      cp ${OutputPrefix}.${name}.Intersect.bed1.count ${OutputPrefix}.${name}.Intersect.bed1.freq
    else
      awk -v SumNum=$SumNum '{OFS="\t";print $1,$2/SumNum*100}' ${OutputPrefix}.${name}.Intersect.bed1.count > ${OutputPrefix}.${name}.Intersect.bed1.freq
    fi

    #Write to the master table:
    awk '{print $2}' ${OutputPrefix}.${name}.Intersect.bed1.count|tr '\n' ' '|awk -v name=$name '{OFS="\t";print name,$1,$2,$3,$4,$5,$6,$7}' >> ${OutputPrefix}.summary.count
    awk '{print $2}' ${OutputPrefix}.${name}.Intersect.bed1.freq|tr '\n' ' '|awk -v name=$name '{OFS="\t";print name,$1,$2,$3,$4,$5,$6,$7}' >> ${OutputPrefix}.summary.freq
    #Cleanning up
    rm -rf ${OutputPrefix}.${name}.Intersect.bed1.temp ${OutputPrefix}.${name}.Intersect.bed1.temp.col1 ${OutputPrefix}.${name}.Intersect.bed1.fill
    rm -rf ${OutputPrefix}.${name}.bed12 ${OutputPrefix}.${name}.Intersect ${OutputPrefix}.${name}.Intersect.bed1.count ${OutputPrefix}.${name}.Intersect.bed1.freq ${OutputPrefix}.${name}.Intersect.bed1
  done
  rm -rf ${Data}.TempProIndi${Seed}.Shifted ${Data}.TempProIndi${Seed}.Shifted.5end ${Annotation}.${Seed}.list
  rm -rf ${Annotation}.${Seed}.ref
  echo "5. Done"
fi
