# This is snakefile procedure to make BAM files

# load config file
# config file contains some path or tools environment variables
# .json file contains a description of the properties of the sample, such as Read group identifier, Library ...

include: "conf/config"
configfile: "conf/project.json"

# create sample name list
SAMPLES=[line.rstrip() for line in open("conf/sm.list")]


# run
rule all:
	input:
		expand("{sample}.done",sample=SAMPLES),
		expand("dataQC/Sample_{sample}/{sample}_combined_R{i}_fastqc.html",sample=SAMPLES,i=[1,2])


rule fastQC:
	input:
		fastq="fastq/Sample_{sample}/{sample}_combined_R{i}.fastq.gz"
	output:	 
		"dataQC/Sample_{sample}/{sample}_combined_R{i}_fastqc.html"
	params:
		outDir="dataQC/Sample_{sample}"
	threads: 5
	shell:
		"fastqc --noextract --format=fastq --threads={threads} --outdir {params.outDir} {input.fastq}"


rule bwa_map:
	input:
		"fastq/Sample_{sample}/{sample}_combined_R1.fastq.gz",
		"fastq/Sample_{sample}/{sample}_combined_R2.fastq.gz"
	output:
		temp("mapped/{sample}.bam")
	params:
		rg="@RG\\tID:{sample}\\tSM:{sample}",
		lb=lambda wildcards:config["samples"][wildcards.sample]["lb"],
		pl=lambda wildcards:config["samples"][wildcards.sample]["pl"],
		pu=lambda wildcards:config["samples"][wildcards.sample]["pu"],
		cn=lambda wildcards:config["samples"][wildcards.sample]["cn"]
	threads:5
	shell:
		"bwa mem -t {threads}"
		" -R '{params.rg}\\tLB:{params.lb}\\tPL:{params.pl}\\tPU:{params.pu}\\tCN:{params.cn}'"
		" {INDEX_PATH} {input}"
		" | samtools view -Sb - > {output}"


rule sort_bam:
	input:
		rules.bwa_map.output
	output:
		bam=temp("mapped/{sample}_sorted.bam"),
		bai=temp("mapped/{sample}_sorted.bai")
	threads: 5
	shell:
		"java -Xmx40g -jar {picard} SortSam"
		" -I {input}"
		" -O {output.bam}"
		" --SORT_ORDER coordinate"
		" --CREATE_INDEX true"
		" --VALIDATION_STRINGENCY SILENT"


rule dedup:
	input:
		bam="mapped/{sample}_sorted.bam",
		bai="mapped/{sample}_sorted.bai"
	output:
		bam="mapped/{sample}/{sample}_sorted_dedup.bam",
		bai="mapped/{sample}/{sample}_sorted_dedup.bai",
		metrics="mapped/{sample}/{sample}_sorted_dedup.metrics"
	threads: 5
	shell:
		"java -Xmx40g -jar {picard} MarkDuplicates \
			--INPUT {input.bam} \
			--OUTPUT {output.bam}  \
			--METRICS_FILE {output.metrics} \
			--REMOVE_DUPLICATES true \
			--ASSUME_SORTED true \
			--CREATE_INDEX true \
			--VALIDATION_STRINGENCY SILENT \
		"

rule collectMultipleMetrics:
	input:
		bam="mapped/{sample}_sorted.bam",
		bai="mapped/{sample}_sorted.bai"
	output:
		touch("dataQC/Sample_{sample}/{sample}.CollectMultipleMetrics.done")
	params:
		prefix="dataQC/Sample_{sample}/{sample}"
	threads: 5
	shell:
		"java -Xmx20g -jar {picard} CollectMultipleMetrics \
			-I {input.bam} \
			-O {params.prefix} \
			-R {REF_FA} \
		"


rule DONE:
	input:
		rules.dedup.output,
		rules.collectMultipleMetrics.output
	output:
		touch("{sample}.done")
