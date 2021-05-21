#CalculateRRS_BH.R
#This is an R script that calculates Ribosome Release Score (RRS) based on RPF.bam and Regions.bed12 inputs, for bluehive only
#Objects in this script:
#    GAlignments: The direct input from a bam file to GAlignments object
#    GRanges: An object that stores coordinate information of many transcripts
#    GRangesList: An list object that contains multiple Tx, but you must iterate each Tx as an independent record
#Inputs: 
#    Data.RPF.bam: Converted from a Data.bed13.RPF file, to Data.bed13.RPF.bed1.bed12.sorted.bam
#    Regions.bed12.utr3/cds: Annotation in bed12 Format, 3utr and cds
#    Data.RNAseq.bam: This is optional
#Example: Rscript CalculateRRS_BH.R $RPFseq $RNAseq ${OutputName}.RRS.cds.matched ${OutputName}.RRS.utr3 ${OutputName}
#Version: Yu Sun, 2018-06-26

args <- commandArgs(TRUE)
InputRPFBam <- args[1]
InputRNABam <- args[2]
AnnotationCDS <- args[3]
AnnotationUTR3 <- args[4]
OutputName <-  args[5]

library(GenomicAlignments)
library(Guitar)
library(ORFik)
print("Extracting BAM from RPF...")
RPFaln <- readGAlignments(InputRPFBam)
if (InputRNABam != "unassigned"){
  print("Extracting BAM from RNAseq...")
  RNAaln <- readGAlignments(InputRNABam)
}
AnnoCDS <-  BED12toGRangesList(AnnotationCDS,header = F)     #This is a GRangesList object
AnnoUTR3 <- BED12toGRangesList(AnnotationUTR3,header = F)
Length <- length(AnnoCDS)
names(AnnoCDS) <- seq(1:Length)
names(AnnoUTR3) <- seq(1:Length)
print("Calculating RRS...")
if (InputRNABam == "unassigned"){
  final <- ribosomeReleaseScore(AnnoCDS, RPFaln, AnnoUTR3, RNA = NULL)
}else{
  final <- ribosomeReleaseScore(AnnoCDS, RPFaln, AnnoUTR3, RNAaln)
}
write.table(final,file=paste(OutputName,".RRS.temp.txt",sep = ""),quote=F,col.names = F, row.names = F)
