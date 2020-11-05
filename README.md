# RiboMining

Mining of Ribo-seq data: Translation efficiency (TE), FLOSS score and Ribosome release score (RRS) calculation

## Translation efficiency (TE)

Translation efficiency (TE) was first defined by Li et. al in their [Cell paper](https://www.sciencedirect.com/science/article/pii/S0092867414002323)

The rate of protein synthesis per mRNA (TE), as measured by protein synthesis rates (from ribosome profiling) divided by mRNA levels (from mRNA-seq)

Here we provide one script to calculate TE:

```
[Usage]
```

## Fragment length organization similarity score (FLOSS)

Fragment length organization similarity score (FLOSS) was first defined by Ingolia et.al in their [Cell Reports paper] (https://www.sciencedirect.com/science/article/pii/S2211124714006299?via%3Dihub)

FLOSS looks at the similarity between a given ribosome footprint (RFP) length distribution and a reference (RFP length distribution on mRNA CDS regions). It measures the magnitude of disagreement between these two distributions, with lower scores reflecting higher similarity.

```
[Usage]
```

## Ribosome release score (RRS)

Ribosome release score (RRS) was first defined by Guttman et.al in their [Cell paper](https://www.sciencedirect.com/science/article/pii/S0092867413007113)

RRS identifies functional protein-coding transcripts with greater sensitivity by detecting the termination of translation at the end of an ORF.

```
[Usage]
```

