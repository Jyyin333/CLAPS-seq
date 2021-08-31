# up-analysis For TREAT samples

incude: "conf/config"
SAMPLES=[line.rstrip() for line in open("conf/sm.Treat.list")]

# Specifies which chromosomes to keep
contigs=list(map(str, range(1,23)))+['X']


rule all:
	input:
		expand("{sample}.AllDone",sample=SAMPLES)


rule fetch_damage_by_chrom:
	input:
		bam="mapped/{sample}/{sample}_sorted_dedup.bam",
		bai="mapped/{sample}/{sample}_sorted_dedup.bai"
	output:
		bed=temp("{sample}_damage.chr{i}.raw.bed")
	params:
		chrom="chr{i}"
	shell:
		"python {scripts}/OG_damage_fetch_by_chrom.py \
		--ref_fa {REF_FA} \
		--bam {input.bam} \
		--out_file {output.bed} \
		--chrom {params.chrom} \
		"


rule merge_bed:
	input:
		["{sample}_damage.chr"+ i + ".raw.bed" for i in contigs]
	output:
		temp("{sample}_damage.rmdup.bed")
	shell:
		"cat {input} > {output}"


rule count_single_base:
	input:
		"{sample}_damage.rmdup.bed"
	output:
		"Basecontext/10bp/{sample}.basecount.txt"
	shell:
		"python {scripts}/count_base.py {input} {output}"


rule fetch_G:
	input:
		"{sample}_damage.rmdup.bed"
	output:
		temp("{sample}_damage.rmdup.G.bed")
	shell:
		"python {scripts}/select_G.py {input} {output}"


rule bedSort:
	input:
		"{sample}_damage.rmdup.G.bed"
	output:
		temp("{sample}_damage.rmdup_sorted.G.bed")
	shell:
		"{bedSort} {input} {output}"


rule count_3bases:
	input:
		rules.bedSort.output
	output:
		"Basecontext/3bp/{sample}.3bases.txt"
	shell:
		"python {scripts}/count_3base.py {input} {output}"


rule bgzip:
	input:
		rules.bedSort.output
	output:
		"{sample}_damage.rmdup_sorted.G.bed.gz"
	shell:
		"{bgzip} -c {input} > {output}"


rule tabix:
	input:
		rules.bgzip.output
	output:
		"{sample}_damage.rmdup_sorted.G.bed.gz.tbi"
	shell:
		"{tabix} {input}"


rule bedtobam:
	input:
		bed=rules.bgzip.output,
		tbi=rules.tabix.output
	output:
		bam="{sample}_damage.rmdup_sorted.G.bam",
		bai="{sample}_damage.rmdup_sorted.G.bam.bai"
	shell:
		"bedtools bedtobam"
		" -g {Genome}"
		" -mapq 60"
		" -i {input.bed}"
		" > {output.bam} && sleep 5s && "
		" samtools index {output.bam}"


rule bamtobw:
	input:
		rules.bedtobam.output.bam
	output:
		"{sample}_damage.rmdup_sorted.G.rpkm.bw"
	threads: 6
	shell:
		"bamCoverage --Offset 5 -p {threads} --normalizeUsing RPKM -b {input} -of bigwig -o {output}"


rule done:
	input:
		rules.bamtobw.output,
		rules.count_single_base.output,
		rules.count_3bases.output
	output:
		touch("{sample}.AllDone")
