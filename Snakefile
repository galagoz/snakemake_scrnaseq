# alignment
# before running, format your fastq file names as the following:
# id_SRR..._S1_L001_R1_001.fastq.gz ("id" should be the same for all fastqs
# and must be used in count commmand below)

rule all:
    input:
        #"test1",
        #"/work/project/becstr_008/data/bams",
        #expand("{file}",file=[f.split("-")[0] for f in os.listdir("/work/project/becstr_008/data/bams") if f.endswith('-1.bam')])
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/counts/{file}/counts.txt",file=(f for f in os.listdir("/work/project/becstr_008/data/bams") if f.endswith('-1.bam'))),
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_2/{file}",file=(f for f in os.listdir("/work/project/becstr_008/data/bams") if f.endswith('-1.bam'))) #if doesn't work, try picard_1!!!!
        expand("/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_2/{bam}.bam",bam=(f.split(".")[0] for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_2/") if f.endswith('-1.bam')))

#rule alignment:
#    input:
#        fastqs="/work/project/becstr_008/data/fastq",
#        transcriptome="/work/project/becstr_008/refdata/refdata-cellranger-GRCh38-1.2.0"
#    output:
#        "test1"
#    threads:
#        16
#    shell:
#        """
#        module load ngs/cellranger/2.1.1
#        cellranger count --fastqs={input.fastqs} --sample=SRR6782109,SRR6782110,SRR6782111,SRR6782112 --id=test --transcriptome={input.transcriptome} --localmem=240
#        """

#rule bamCleave:
#    input:
#        bam="/work/project/becstr_008/poirion_snakemake/test/outs/possorted_genome_bam.bam"
#    output:
#        out="/work/project/becstr_008/data/bams/"
#    threads:
#        16
#    shell:
#        """
#        ./bamCleave -b {input.bam} -o {output.out} -t CB -c 1008
#        """

#rule featureCounts:
#    input:
#        gtfDir="/work/data/genomes/human/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf",
#        inDir="/work/project/becstr_008/data/bams/{file}"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/counts/{file}/counts.txt"
#    threads:
#        16
#    shell:
#        """
#        /work/project/becstr_008/subread-1.4.6-p5-Linux-x86_64/bin/featureCounts -p -B -C -M --primary -T 16 -t exon -g gene_id -a {input.gtfDir} -o {output} {input.inDir}
#        """

#rule addOrReplaceReadGroups:
#    input:
#        "/work/project/becstr_008/data/bams/{file}"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_1/{file}"
#    threads:
#        8
#    shell:
#        """
#        java -jar /work/project/becstr_008/picard.jar AddOrReplaceReadGroups I={input} O={output} SORT_ORDER=coordinate RGID={wildcards.file} RGLB=Homo_sapiens RGPL=illumina RGPU={wildcards.file} RGSM={wildcards.file}
#        """

#rule markDuplicates:
#    input:
#        "/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_1/{file}"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_2/{file}"
#    threads:
#        8
#    shell:
#        """
#        java -jar /work/project/becstr_008/picard.jar MarkDuplicates I={input} O={output} CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={output}.marked_dup_metrics.txt
#        """

rule buildBamIndex:
    input:
        "/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_2/{file}.bam"
    output:
        "/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_3/{file}.bai"
    threads:
        8
    shell:
        """
        java -jar /work/project/becstr_008/picard.jar BuildBamIndex I={input} O={output} TMP_DIR="/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_3/tmpDir"
        """
