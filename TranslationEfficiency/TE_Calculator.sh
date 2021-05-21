#!/usr/bin/env bash
#TE_Calculator.sh
#This scripts takes RPF.bed13 and RNAseq.unique.bed12 files as input, calculates Translation Efficiency (TE) using CDS regions of a set of transcripts
#Formula of TE: (density of RPF within ORF) / (RNA expression of ORFs transcript), optional: adding a pseudo count
#The R script uses Bioconductor ORFik package, following this manual: https://rdrr.io/github/JokingHero/ORFik/man/translationalEff.html
#Version: Yu Sun, 2019-07-17
#Version: Yu Sun, 2019-08-05. I realized that this Bioconductor function cannot detect the stranded RNAseq type, and only count the reads mapped to the sense of the trasncript.
#                             This caused errors when single-end dUTP stranded RNAseq reads bam file containing majority of the reads mapped to the reverse strand
#                             I decided to use RNAseq unique bed12 file as input, so that I can control the strand info

array=( "$@" )

if [ ! -n "$1" ]
then
  echo "    This scripts takes RPF.bed13 and RNAseq.unique.bed12 files as input, calculates Translation Efficiency (TE) using CDS regions of a set of transcripts"
  echo "    Formula: (density of RPF within ORF in rpkm +1) / (RNA expression of ORFs transcript in rpkm +1), the +1 is the pseudo count, and it's adjustable"
  echo "    The bed12 gene annotation should have ORF information in col7-8. If col7=col8, those transcripts will be ignored in analysis."
  echo "    Usage: `basename $0` -rpf [RPF.bed13] -g [Transcript.bed12] -rna [RNAseq.unique.bed12] -o [OutputPrefix]"
  echo "        Optional: -genome: The default is mm10, but if you are running other species' data, you have to specify the genome, such as galGal6, anoCar2, etc."
  echo "        Optional: -pseudo: The default pseudo count is 1 added on RPKM values. You can adjust it."
  echo "        `basename $0` -rpf [RPF.unique.bed13] -g [Transcript.bed12] -rna [RNAseq.unique.bed12] -o [OutputPrefix] -genome [mm10/galGal6/etc.] -pseudo [Number]"
  echo "    Output: OutputPrefix.TE.txt, OutputPrefix.RPKM.TE.txt"
else
  module load bedtools
  module load samtools
  module load r/3.5.0/b1
  RPFseq="unassigned"
  RNAseq="unassigned"
  Transcript="unassigned"
  OutputName="unassigned"
  genome="unassigned"
  pseudocount="unassigned"
  TEMP=$RANDOM
  SECONDS=0

echo "1. Getting inputs"
for arg in "$@"
do
  if [[ $arg == "-rpf" ]]
    then
      RPFseq=${array[$counter+1]}
      echo '   RPFseq: '$RPFseq
  elif [[ $arg == "-rna" ]]
    then
      RNAseq=${array[$counter+1]}
      echo '   RNAseq: '$RNAseq
  elif [[ $arg == "-g" ]]
    then
      Transcript=${array[$counter+1]}
      echo '   Transcript: '$Transcript
  elif [[ $arg == "-o" ]]
    then
      OutputName=${array[$counter+1]}
      echo '   OutputPrefix: '$OutputName
  elif [[ $arg == "-genome" ]]
    then
      genome=${array[$counter+1]}
      echo '   genome: '$genome
  elif [[ $arg == "-pseudo" ]]
    then
      pseudocount=${array[$counter+1]}
      echo '   pseudocount: '$pseudocount
  fi
  let counter=$counter+1
done

if [[ $genome == "unassigned" ]];then
    genome="mm10"
    echo '   use default genome: mm10'
fi
if [[ $pseudocount == "unassigned" ]];then
    pseudocount=1
    echo '   use default pseudocount: 1'
fi

if [[ $RPFseq == "unassigned" || $Transcript == "unassigned" || $OutputName == "unassigned" || $RNAseq == "unassigned" ]];then
    echo "     >>> [Error]: Input missing! Exitting..."
    exit 1
fi

echo "2. Extracting ORF CDS"
./bin/BED12Extractor.sh -a cds -i $Transcript -o ${OutputName}.TE.cds 1> /dev/null
if [ ! -s ${OutputName}.TE.cds ];then
    echo "   Problem occurred... Exitting..."
    exit 1
fi
CDSLines=`wc -l ${OutputName}.TE.cds|awk '{print $1}'`
TxLines=`wc -l $Transcript|awk '{print $1}'`
if [[ $CDSLines != $TxLines ]];then
    echo "   [Warning!] Some ORFs were not extracted (the bed12 file is probably wrong)... Those transcripts will be missed..."
fi
awk '{print $4}' ${OutputName}.TE.cds > ${OutputName}.TE.cds.list

echo "   Getting Tx length, outputting ${OutputName}.TE.cds.len.txt"
#This program has been modified, and it is fast now:
./bin/CountBed12_Length.sh ${OutputName}.TE.cds

echo "3. Converting RPF.bed13 to RPF.bam"
  echo "   bed22bed1"
  awk '{for (i=1;i<=$4;i++) print $0}' $RPFseq > ${RPFseq}.${TEMP}.bed1
  echo "   Get bed12"
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ${RPFseq}.${TEMP}.bed1  > ${RPFseq}.${TEMP}.bed1.bed12
  RPFDepth=`wc -l ${RPFseq}.${TEMP}.bed1.bed12|awk '{print $1}'`
  echo "   Converting to bam"
  bedToBam -bed12 -g ./data/${genome}.ChromInfo.txt -i ${RPFseq}.${TEMP}.bed1.bed12 > ${RPFseq}.${TEMP}.bed1.bed12.bam
  echo "   Sorting bam"
  samtools sort ${RPFseq}.${TEMP}.bed1.bed12.bam > ${RPFseq}.${TEMP}.bed1.bed12.sorted.bam
  echo "   Indexing"
  samtools index ${RPFseq}.${TEMP}.bed1.bed12.sorted.bam

echo "4. Converting RNAseq.bed12 to RNAseq.bam"
#I didn't modify the col4 since this won't affect further analysis. No need to do bed22bed1
  echo "   Converting to bam"
  RNADepth=`wc -l ${RNAseq}|awk '{print $1}'`
  bedToBam -bed12 -g ./data/${genome}.ChromInfo.txt -i ${RNAseq} > ${RNAseq}.${TEMP}.bam
  echo "   Sorting bam"
  samtools sort ${RNAseq}.${TEMP}.bam > ${RNAseq}.${TEMP}.sorted.bam
  echo "   Indexing"
  samtools index ${RNAseq}.${TEMP}.sorted.bam

echo "4. Calculating TE"
echo "   Running: Rscript CalculateTE_BH.R ${RPFseq}.${TEMP}.bed1.bed12.sorted.bam ${RNAseq}.${TEMP}.sorted.bam ${OutputName}.TE.cds ${OutputName}.TE.cds.len.txt ${RPFDepth} ${RNADepth} $pseudocount ${OutputName}"
echo "   Using pseudocount: "$pseudocount
Rscript ./bin/CalculateTE_BH.R ${RPFseq}.${TEMP}.bed1.bed12.sorted.bam ${RNAseq}.${TEMP}.sorted.bam ${OutputName}.TE.cds ${OutputName}.TE.cds.len.txt ${RPFDepth} ${RNADepth} $pseudocount ${OutputName} 2> /dev/null
touch ${OutputName}.TE.temp.txt
LineOut=`wc -l ${OutputName}.TE.temp.txt|awk '{print $1}'`
#echo $LineOut
if [[ $CDSLines == $LineOut ]];then
    echo "   Results generated"
    paste ${OutputName}.TE.cds.list ${OutputName}.TE.temp.txt | awk '{OFS="\t";print $1,$2,$3,$4,$5,$6,$7}' > ${OutputName}.TE.temp2.txt
    echo -e "Transcript\tLength\tRPF.Count\tRNAseq.Count\tRPF.RPKM\tRNAseq.RPKM\tTE" > ${OutputName}.header
    cat ${OutputName}.header ${OutputName}.TE.temp2.txt > ${OutputName}.RPKM.TE.txt
    awk 'NR>1' ${OutputName}.RPKM.TE.txt | awk '{OFS="\t";print $1,$7}' > ${OutputName}.TE.txt
else
    echo "   [Error!] Error occurred! Exitting..."
fi

echo "5. Finalizating results."
rm -rf ${OutputName}.header ${OutputName}.TE.cds.list ${OutputName}.TE.temp.txt ${OutputName}.TE.temp2.txt ${OutputName}.TE.cds.fa ${OutputName}.TE.cds.fa.len ${OutputName}.TE.cds.len.txt
rm -rf ${RPFseq}.${TEMP}.bed1 ${RPFseq}.${TEMP}.bed1.bed12 ${RPFseq}.${TEMP}.bed1.bed12.bam ${RNAseq}.${TEMP}.bam
rm -rf ${RPFseq}.${TEMP}.bed1.bed12.sorted.bam* ${RNAseq}.${TEMP}.sorted.bam* ${OutputName}.TE.cds.fa ${OutputName}.TE.cds

duration=$SECONDS
echo "   $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "   Done"

fi
