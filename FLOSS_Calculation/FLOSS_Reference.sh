#!/usr/bin/env bash
#FLOSS_Reference.sh
#This scipt takes RPF bed13 data and a reference annotation as input, calculate a reference length distribution
#The reference can be CDS, UTR3, UTR5 of a set of transcripts
#The method is based on Ingolia's FLOSS paper
#Version: Yu Sun, 2018-11-06

if [ ! -n "$1" ]
then
  echo "    This scipt takes RPF bed13 data and a reference annotation as input, calculate a reference length distribution"
  echo "    The reference can be CDS, UTR3, UTR5 of a set of transcripts"
  echo "    The method is based on Ingolia's FLOSS paper"
  echo "    Length range: [26,32], by reads"
  echo "    Usage: `basename $0` [Data.bed13] [Reference.bed12] [OutputPrefix] [Offset|default 12]"
  echo "    Output: OutputPrefix.count and OutputPrefix.freq files"
else
  Data=$1
  Annotation=$2
  OutputPrefix=$3
  Seed=$RANDOM$RANDOM
  if [ ! -n "$4" ];then POffset=12;else POffset=$4;fi

  echo "1. Starting Floss reference length distribution estimation"
  echo "  Data: "$Data
  echo "  Annotation: "$Annotation
  echo "  OutputPrefix: "$OutputPrefix
  echo "  P offset: "$POffset

  echo "2. Shifting RPF 5end:"
  Shift5endOffset_BED13RPF.py $Data $POffset ${Data}.TempPro${Seed}.Shifted
  
  echo "3. Get 5end"
  awk 'BEGIN{OFS="\t"}{if ($6=="+") print $1,$2,$2+1,$4,$5,$6,$2,$2+1,$9,"1","1","0",$13;else print $1,$3-1,$3,$4,$5,$6,$3-1,$3,$9,"1","1","0",$13}' ${Data}.TempPro${Seed}.Shifted > ${Data}.TempPro${Seed}.Shifted.5end
  #awk '{for (i=1;i<=$4;i++) print $0}' ${Data}.TempPro${Seed}.Shifted.5end > ${Data}.TempPro${Seed}.Shifted.5end.bed1
  
  echo "4. Intersecting Annotation"
  bedtools intersect -u -wa -split -b $Annotation -a ${Data}.TempPro${Seed}.Shifted.5end > ${Data}.TempPro${Seed}.Shifted.5end.Intersect
  
  echo "5. Converting to reads"
  awk '{for (i=1;i<=$4;i++) print $0}' ${Data}.TempPro${Seed}.Shifted.5end.Intersect > ${Data}.TempPro${Seed}.Shifted.5end.Intersect.bed1

  echo "6. Counting..."
  awk '{print length($13)}' ${Data}.TempPro${Seed}.Shifted.5end.Intersect.bed1|sort|uniq -c|awk '{OFS="\t";print $2,$1}' > ${OutputPrefix}.FlossRef.temp
  for i in {26..32};do echo $i;done > ${OutputPrefix}.ref
  awk '{print $1}' ${OutputPrefix}.FlossRef.temp > ${OutputPrefix}.FlossRef.temp.col1
  cat ${OutputPrefix}.ref ${OutputPrefix}.FlossRef.temp.col1|sort|uniq -c|awk '$1==1'|awk '{OFS="\t";print $2,0}' > ${OutputPrefix}.fill
  cat ${OutputPrefix}.FlossRef.temp ${OutputPrefix}.fill|sort -k1,1 -n > ${OutputPrefix}.count
  SumNum=`awk '{SUM+=$2}END{printf "%.0f\n",SUM}' ${OutputPrefix}.count`
  awk -v SumNum=$SumNum '{OFS="\t";print $1,$2/SumNum*100}' ${OutputPrefix}.count > ${OutputPrefix}.freq

  rm -rf ${Data}.TempPro${Seed}.Shifted ${Data}.TempPro${Seed}.Shifted.5end ${Data}.TempPro${Seed}.Shifted.5end.Intersect ${Data}.TempPro${Seed}.Shifted.5end.Intersect.bed1
  rm -rf ${OutputPrefix}.FlossRef.temp ${OutputPrefix}.FlossRef.temp.col1 ${OutputPrefix}.fill ${OutputPrefix}.ref 

fi
