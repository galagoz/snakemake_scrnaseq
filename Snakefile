# alignment
# before running, format your fastq file names as the following:
# id_SRR..._S1_L001_R1_001.fastq.gz ("id" should be the same for all fastqs
# and must be used in count commmand below)

FILES = (f for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_2/") if f.endswith('-1.bam'))
BAMS = (f.split(".")[0] for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_2/") if f.endswith('-1.bam'))

rule all:
    input:
        #"test1",
        #"/work/project/becstr_008/data/bams",
        #expand("{file}",file=[f.split("-")[0] for f in os.listdir("/work/project/becstr_008/data/bams") if f.endswith('-1.bam')])
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/counts/{file}/counts.txt",file=FILES),
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_2/{file}",file=FILES), #if doesn't work, try picard_1!!!!
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_3/{bam}.bai",bam=BAMS),
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_4/{file}",file=FILES),
        #"/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dict",
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_1/{file}",file=(f for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_6/") if f.endswith('-1.bam'))),
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_1/{file}.ALLforIndelRealigner.intervals",file=(f for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_1/") if f.endswith('-1.bam'))),
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_3/{file}",file=(f for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_1/") if f.endswith('-1.bam'))),
        expand("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_3/{file}.recal_data.csv",file=(f for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_3/") if f.endswith('-1.bam')))

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
#        bam="/work/project/becstr_008/snakemake_scrnaseq/test/outs/possorted_genome_bam.bam"
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

#####################
#       PICARD
#####################

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

#rule buildBamIndex:
#    input:
#        "/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_2/{bam}.bam"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_3/{bam}.bai"
#    threads:
#        8
#    shell:
#        """
#        java -jar /work/project/becstr_008/picard.jar BuildBamIndex I={input} O={output} TMP_DIR=/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_3/tmpDir
#        """

#rule sortSam:
#    input:
#        "/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_2/{file}"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_4/{file}"
#    threads:
#        8
#    shell:
#        """
#        java -jar /work/project/becstr_008/picard.jar SortSam I={input} O={output} SORT_ORDER=coordinate TMP_DIR=/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_4/tmpDir CREATE_INDEX=TRUE
#        """

#rule createSequenceDictionary:
#    input:
#        "/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
#    output:
#        "/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dict"
#    threads:
#        16
#    shell:
#        """
#        java -jar /work/project/becstr_008/picard.jar CreateSequenceDictionary R={input} O={output}
#        """

#rule reorderSam:
#    input:
#        file="/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_4/{file}",
#        ref="/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
#        dict="/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.dict"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_6/{file}"
#    threads:
#        8
#    shell:
#        """
#        java -jar /work/project/becstr_008/picard.jar ReorderSam I={input.file} O={output} REFERENCE_SEQUENCE={input.ref} SEQUENCE_DICTIONARY={input.dict} CREATE_INDEX=TRUE
#        """

#####################
#       GATK
#####################

#rule splitNCigarReads:
#    input:
#        file="/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_6/{file}",
#        ref="/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_1/{file}"
#    threads:
#        16
#    shell:
#        """
#        java -jar /work/project/becstr_008/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T SplitNCigarReads -I {input.file} -o {output} -R {input.ref} -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
#        """

#rule realignerTargetCreator:
#    input:
#        file="/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_1/{file}",
#        ref="/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
#        vcf="/work/project/becstr_008/data/indel_calls/00-All.sorted.vcf"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_1/{file}.ALLforIndelRealigner.intervals"
#    threads:
#        16
#    shell:
#        """
#        java -jar /work/project/becstr_008/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T RealignerTargetCreator -I {input.file} --known {input.vcf} -o {output} -R {input.ref} -nt 16
#        """

#rule indelRealigner:
#    input:
#        file="/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_1/{file}",
#        ref="/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
#        vcf="/work/project/becstr_008/data/indel_calls/00-All.sorted.vcf"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_3/{file}"
#    threads:
#        16
#    shell:
#        """
#        java -jar /work/project/becstr_008/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T IndelRealigner -I {input.file} -o {output} -targetIntervals {input.file}.ALLforIndelRealigner.intervals -R {input.ref}
#        """

rule baseRecalibrator:
    input:
        file="/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_3/{file}",
        ref="/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        vcf="/work/project/becstr_008/data/indel_calls/00-All.sorted.vcf"
    output:
        "/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_3/{file}.recal_data.csv"
    threads:
        16
    shell:
        """
        java -jar /work/project/becstr_008/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T BaseRecalibrator -I {input.file} -o {output} -R {input.ref} -nct 20 -knownSites {input.vcf}
        """
