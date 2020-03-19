# alignment
# before running, format your fastq file names as the following:
# id_SRR..._S1_L001_R1_001.fastq.gz ("id" should be the same for all fastqs
# and must be used in count commmand below)

FILES = (f for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_2/") if f.endswith('-1.bam'))
BAMS = (f.split(".")[0] for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/picard/picard_2/") if f.endswith('-1.bam'))

#rule all:
    #input:
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
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_3/{file}.recal_data.csv",file=(f for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_3/") if f.endswith('-1.bam'))),
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_5/{file}",file=(f for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_3/") if f.endswith('-1.bam'))),
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_6/{file}.snp.vcf",file=(f for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_5/") if f.endswith('-1.bam'))),
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_7/{file}/snv_filtered.vcf",file=(f for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_6/") if f.endswith('.vcf'))),
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_8/{file}/results.table",file=(f for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_7/") if f.endswith('.vcf'))),
        #"/work/project/becstr_008/results/poirion_snakemake/counts/keys.txt",
        #"/work/project/becstr_008/results/poirion_snakemake/counts/values.txt"
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/counts/{file}/count_table.txt",file=(f for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/counts/") if f.endswith('.bam')))
        #expand("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_8/{vcfs}.snp.vcf/results.table",vcfs=(f for f in os.listdir("/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_8/") if f.endswith('.vcf')))

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

#rule baseRecalibrator:
#    input:
#        file="/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_3/{file}",
#        ref="/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
#        vcf="/work/project/becstr_008/data/indel_calls/00-All.sorted.vcf"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_3/{file}.recal_data.csv"
#    threads:
#        16
#    shell:
#        """
#        java -jar /work/project/becstr_008/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T BaseRecalibrator -I {input.file} -o {output} -R {input.ref} -nct 20 -knownSites {input.vcf}
#        """

#rule printReads:
#    input:
#        file="/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_3/{file}",
#        ref="/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_5/{file}"
#    threads:
#        16
#    shell:
#        """
#        java -jar /work/project/becstr_008/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T PrintReads -I {input.file} --out {output} -R {input.ref} -BQSR {input.file}.recal_data.csv -nct 20
#        """

#rule haplotypeCaller:
#    input:
#        file="/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_5/{file}",
#        ref="/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
#        vcf="/work/project/becstr_008/data/indel_calls/00-All.sorted.vcf"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_6/{file}.snp.vcf"
#    threads:
#        16
#    shell:
#        """
#        java -jar /work/project/becstr_008/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T HaplotypeCaller -I {input.file} -o {output} -R {input.ref} --dbsnp {input.vcf} -dontUseSoftClippedBases
#        """

#rule variantFiltration:
#    input:
#        file="/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_6/{file}",
#        ref="/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_7/{file}/snv_filtered.vcf"
#    threads:
#        16
#    shell:
#        """
#        java -jar /work/project/becstr_008/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantFiltration -V {input.file} -o {output} -R {input.ref} -cluster 3 --filterExpression "FS > 30.0" --filterName "FS" --filterExpression "QD < 2.0" --filterName "QD"
#        """

#rule variantToTable:
#    input:
#        file="/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_7/{file}/snv_filtered.vcf",
#        ref="/work/project/becstr_008/data/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
#    output:
#        "/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_8/{file}/results.table"
#    threads:
#        16
#    shell:
#        """
#        java -jar /work/project/becstr_008/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T VariantsToTable -V {input.file} -R {input.ref} -F CHROM -F POS -F ID -F QUAL -F AC -o {output}
#        """

#rule keysAndValues:
#    input:
#        counts="/work/project/becstr_008/results/poirion_snakemake/results/counts/bams_sel_TTTGTCACATGCCACG-1.bam/counts.txt"
#    output:
#        keys="/work/project/becstr_008/results/poirion_snakemake/results/counts/keys.txt",
#        values="/work/project/becstr_008/results/poirion_snakemake/results/counts/values.txt"
#    shell: # Save gene IDs as "keys.txt" & other values as "values.txt"
#        """
#        awk 'NR>=3 {{print $1}}' {input.counts} > {output.keys}
#        awk 'NR>=3 {{$1=""; print $0}}' {input.counts} > {output.values}
#        """

rule countTable:
    input:
        file="/work/project/becstr_008/results/poirion_snakemake/results/counts/bams_sel_TGAGCCGTCGCAAACT-1.bam"
    output:
        count_table="/work/project/becstr_008/results/poirion_snakemake/results/counts/bams_sel_TGAGCCGTCGCAAACT-1.bam/count_table.txt"
    shell: # Generate count table from counts.txt files & generate variant table from results.table files + merge all snv_filtered.vcf files
           # init the txt file with geneID, geneLength and geneCounts = count_table.txt
           # append the geneCount column of each counts.txt to the end of the count_table2.txt
           # list and save all cell "vcfToTable" paths in gatk_8
           # list and save all SNV_dictionary.txt paths in gatk_8
        """
        shuf -n1 -e * = /work/project/becstr_008/results/poirion_snakemake/results/counts/$randomFolder
        awk 'NR>=2 {{print $1, $6}}' /work/project/becstr_008/results/poirion_snakemake/results/counts/$randomFolder/counts.txt > /work/project/becstr_008/results/poirion_snakemake/results/counts/$randomFolder/count_table.txt
        awk 'NR==FNR{{a[NR]=$7;next}}{{a[FNR]=a[FNR+1];a[1]="cellID"}}{{print $0,a[FNR]}}' {input.file}/test_counts.txt {input.file}/count_table.txt > {input.file}/tmp_counts.txt; mv {input.file}/tmp_counts.txt {output.count_table}

        ls -d /work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_8/*/ > /work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_8/vcfList.txt
        ls -d /work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_8/*/SNV_dictionary.txt > /work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_8/dictList.txt
        """

#vcfs="/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_8/{vcfs}.snp.vcf"
#awk '{{print $0, "chr" $1 "_" $2}}' {input.vcfs}/results.table > {input.vcfs}/tmp_res.table; mv {input.vcfs}/tmp_res.table {output.res_table}
#res_table="/work/project/becstr_008/results/poirion_snakemake/results/gatk/gatk_8/{vcfs}.snp.vcf/results.table"
