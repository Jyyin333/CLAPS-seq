# Codes for CLAPS-seq

This repository stroes original codes used in ***Base resolution mapping of 8-oxoguanine unveils reduced occurrence at G-quadruplexes***.



# Software Requiremnents

## Python
+ Pyhton version >= 3.7.x
+ numpy >= 1.19.2
+ pandas >= 1.1.2
+ py2bit >= 0.3.0
+ pyBigWig >= 0.3.17
+ pyfaidx >= 0.5.9.1
+ pysam >= 0.16.0.1
+ scipy >= 1.6.0
+ re >= 2.2.1

## R
+ R version >= 3.6.1
+ ggplot2 >= 3.3.0
+ stringr >= 1.4.0

## Unix Tools
+ picard >= 2.25.0
+ FastQC >= 0.11.9
+ bgzip >= 1.9
+ tabix >= 1.9
+ samtools >= 1.9
+ deeptools >= 3.5.0
+ bedtools >= 2.29.2
+ bwa >= 0.7.17
+ snakemake >= 5.5.4


---

# Workflow Example

## 1. Quality Check and Alignment
You can downloaded the raw sequencing data from [GEO](https://www.ncbi.nlm.nih.gov/geo/) under GSE181312.  
Type `snakemake -s Align.snakefile`.  
Some configuration file formats can be viewed in [Example](https://github.com/Jyyin333/CLAPS-seq/tree/main/Data%20Example)

Immediately, deeptools function multiBamSummary was used for validating reproducibility as well as robustness:
```
# run multiBamSummary
multiBamSummary bins -b ${bamfiles[@]} -o bam.bins.npz -l ${labels[@]} -bs 10000 --outRawCounts AllBam.bins.rawCounts.mtx
```
At the same time, library sizes for each sample were computed using samtools or other appropriate tools, the lib.sizes file should following this format and in order of the bamfiles used in previous step:
```
114834689
60748386
90203070
109105413
125435253
...
```
Then run plot_corr.r, `Rscript plot_corr.r`.

## 2. Fetch predicted OG sites
After aligning reads to reference genome, we fetched the predicted OG sites.
Type `snakemake -s snOG-treat.up.snakefile` for Treat samples or `snakemake -s snOG-treat.up.snakefile` for control samples.

For this two steps, you can replace BWA with any other align-software. Also you can break up the snakefile and run each step separately.  
You will end up with a BED file, a BAM file as well as a BigWig file. This files record the OG sites or OG signal.

Here we present part of BED file.
```
chr1	10190	10200	TTAGGGTTAG	940	-
chr1	10522	10532	TGAAGGCGGA	940	-
chr1	16346	16356	AGAAGGGGTC	940	-
chr1	30794	30804	GTTAGATGAG	60	+
chr1	30978	30988	CTTCGCATCC	940	-
chr1	51328	51338	AAAAGCCCTC	60	+
chr1	52033	52043	GCATGGTAAA	60	+
```

Furthermore, base context infomation is also counted during this process.  
There are two statistics files. One looks like:
```
	A	C	G	T
0	12798546	9066306	10478653	12061170
1	12559260	8615754	11154373	12075288
2	12119596	8407238	12300852	11576989
3	12656947	3814228	14015600	13917900
4	5489463	1869801	32315569	4729842
5	14521092	10210376	13236254	6436953
6	13789846	10777407	11225875	8611547
7	17217005	8842695	10215494	8129481
8	16081684	11684762	8548633	8089596
9	15131825	11283308	8293559	9695983
```
and the other:
```
GGG	2957222
AGG	2702839
AGA	3244275
CGC	420545
AGC	2402683
TGG	3256454
... ...

```

## 3. OG and G/GC relationship
Type `computeGbias.py -b treat.bam ctrl.bam -g genome.fa -o treat.res` to get an intermediate file:
```
# T_plus  T_minus C_plus  C_minus  T_all  C_all  G_ratio  C_ratio
0.0     0.0     0.0     0.0     0.0     0.0     0.502   0.196
0.0     0.0     0.0     0.0     0.0     0.0     0.518   0.18
6.0     4.0     11.0    2.0     10.0    13.0    0.492   0.202
6.0     5.0     29.0    18.0    11.0    47.0    0.39    0.291
4.0     5.0     21.0    30.0    9.0     51.0    0.335   0.343
5.0     6.0     15.0    23.0    11.0    38.0    0.266   0.408
4.0     2.0     25.0    34.0    6.0     59.0    0.314   0.359
```
Then, using plotGbias.py to plot relationship between OG & G/GC
```
plotGbias.py -i treat.res -o treat.pdf -gc true
```
Note that you can choose which type of figure to plot by setting the -gc parameter.

## 4. Chromatin state
To annotate OG with chromatin states. You should first obtain the annotation file, here we downloaded the [Hela-S3 chromatin state annotation file](https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeAwgSegmentation&db=hg19) from UCSC:
```
chr1	10200	10309	WE	1000	.	10200	10309	255,252,4
chr1	10400	10409	R	1000	.	10400	10409	127,127,127
chr1	10509	10701	R	1000	.	10509	10701	127,127,127
chr1	13960	14850	CTCF	1000	.	13960	14850	10,190,254
chr1	14850	15800	R	1000	.	14850	15800	127,127,127
```
Then run `ChromstateAnno.py -a treat.bam -b ctrl.bam --chromstate Annotationfile.bed.gz -g genome.2bit --outmatrix res.tsv -o boxplot.pdf`, the res.tsv:
```
Chrom	Start	End	Name	gA_mean_counts	gB_mean_counts	GC_ratio
chr1	10200	10309	WE	0.0	0.0	0.5412844036697247
chr1	10400	10409	R	0.0	0.0	0.3333333333333333
chr1	10509	10701	R	1.0	1.0	0.7135416666666666
chr1	13960	14850	CTCF	0.0	1.0	0.5797752808988764
chr1	14850	15800	R	0.0	5.0	0.6105263157894737
chr1	16000	16430	T	1.0	13.0	0.5581395348837209
```

## 5. Metaprofiles
For genes, you can use refpointMatrix.py or scaleMatrix.py to generate two different type of matrices.  
In scaleMatrix.py, it will first shrunk or strech total region into samle length.
```
# make score matrix
refpointMatrix.py -b treat.bam -r genes.bed -refpoint tss -a 3000 -d 3000 -bs 50 -o treat.rp.mtx.gz
scaleMatrix.py -b treat.bam -r genes.bed -a 2000 -d 2000 -m 5000 -bs 50 -o treat.sc.mtx.gz

# plotprofile
plotMatrix.py -t treat.rp.mtx.gz -o treat.rp.pdf
plotMatrix.py -t treat.sc.mtx.gz -o treat.sc.pdf
```
Using following command to normalize OG signal when you have control sample
```
plotMatrix.py -t treat.rp.mtx.gz -c ctrl.rp.mtx.gz -o treat.rp.pdf
```

For other type of Regions, such as DHS, using [deeptools](https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html) to generate profiles.

## 6. Background G/GC content
There are also two type of scripts to visualize background G/GC content. With [reference point](https://github.com/Jyyin333/CLAPS-seq/blob/main/Scripts/bgG_rp.py) mode, you can focus on a specific point(e.g. TSS) and see how G/GC content distributing around it.  As for [scale](https://github.com/Jyyin333/CLAPS-seq/blob/main/Scripts/bgG_sc.py) mode, it's similar to scaleMatrix.py in Step5. See the script for details.
```
bgG_rp.py -r genes.bed -g genome.2bit -refpoint tss -a 3000 -d 3000 -bs 50 -o genes.backgroundGC.rp.pdf

# or

bgG_sc.py -r genes.bed -g genome.2bit -a 2000 -d 2000 -m 5000 -bs 50 -o gene.backgroundG.sc.pdf
```
For PQS region, due to the large amount of data, running this script directly can be very time-consuming. We suggest splitting the files and merging the results.

Note: Either in Step5 or in Step6, in scale mode, these scripts are only applicable to Region files with strand information such as [Gene BED](https://github.com/Jyyin333/CLAPS-seq/blob/main/Data%20Example/genes.bed) file. In other words, the images they produce must be in strand-split format.

## 7. Wilcoxon test
In order to investigate the effects of G4, we classified Promoters(defined as TSS ± 2kb) into two categories depending on whether they overlap with G4-peak or not. Type command line as following:
```
# classify
intersectBed -wa genes.ol.bed -u -a <(awk '{if($6=="+"){print $1"\t"$2-2000"\t"$2+2000} else if($6=="-"){print $1"\t"$3-2000"\t"$3+2000}}' genes.bed) -b G4.narrowpeak
intersectBed -wa genes.uol.bed -v -a <(awk '{if($6=="+"){print $1"\t"$2-2000"\t"$2+2000} else if($6=="-"){print $1"\t"$3-2000"\t"$3+2000}}' genes.bed) -b G4.narrowpeak

# calculate scores
refpointMatrix.py -b treat.bam -r genes.ol.bed -refpoint center -a 3000 -d 3000 -bs 50 -o treat.ol.mtx.gz
refpointMatrix.py -b treat.bam -r genes.uol.bed -refpoint center -a 3000 -d 3000 -bs 50 -o treat.uol.mtx.gz

# replace treat.bam with ctrl.bam and run again

# run Rscript
# specify an interval around TSS for comparison and test through setting parameter --interval, default 1000
wtest.r --treat treat.ol.mtx.gz treat.uol.mtx.gz --control ctrl.ol.mtx.gz ctrl.uol.mtx.gz --outfig boxplot.pdf --interval 1000
```

## 8. Others

### PQS profile
```
# seperate PQS by G4
intersectBed -wa olG4.PQS.bed -u -a ALL_type.PQS.bed -b G4.narrowpeak
intersectBed -wa uolG4.PQS.bed -v -a ALL_type.PQS.bed -b G4.narrowpeak

# bed2bam
bedtools bedtobam -i ALL_type.PQS.bed -g genome > ALL_type.PQS.unsort.bam

# sort bam
samtools sort -o ALL_type.PQS.sort.bam ALL_type.PQS.unsort.bam && samtools index ALL_type.PQS.sort.bam

# make bigwig
bamCoverage -b ALL_type.PQS.sort.bam -o ALL_type.PQS.bw -of bigwig --normalizeUsing None

# compute matrix 
computeMatrix reference-point -S ALL_type.PQs.bw -R olG4.PQS.bed uolG4.PQS.bed -o ALL_type.PQs.inPQS.mtx.gz --referencePoint center -b 3000 -a 3000 -bs 50

# plotprofile
plotProfile -m ALL_type.PQs.inPQS.mtx.gz -o ALL_type.PQs.inPQS.pdf  --regionsLabel "G4+" "G4-"
```

### IGV snapshot
For aesthetic purposes, we remade the OG signal tracks using makeBW.py
```
makeBW.py -b treat.bam -o treat.remake.bw -r chr1:100-200
``` 

Note: This can be a time-consuming operation due to large number of repeated sampling, therefore we do not recommend the construction of genome-wide signal tracks
