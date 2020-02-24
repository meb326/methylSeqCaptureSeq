#!/bin/bash 

#bam2fastq
samtools sort -n in.bam | samtools fastq -1 sample_1.fastq.gz -2 sample_2.fastq.gz -0 sample_0.fastq.gz name-collate.bam

#fastqc
fastqc -o fastqc_out/ sample_1.fastq.gz sample_2.fastq.gz 

#readtrimming  
trim_galore --paired --trim1 sample_1.fastq.gz sample_2.fastq.gz 

#genomeindexing 
biscuit index <path_to_genome_folder> hg19.fa

#indexbwa-meth
bwameth.py index $REFERENCE

#alignment
biscuit align -t 30 <path_to_genome_folder> hg19.fa A1/*1.fq.gz A1/*2.fq.gz | samtools view -Sb > A1.bam

#alignbwa-meth
bwameth.py --threads 16 \
     --reference $REFERENCE \
     $sample_1.fastq.gz $sample_2.fastq.gz > some.sam | samtools view -Sb > A1_bwameth.bam

#flagstat 
samtools flagstat -@ 8 SAMPLE.bam -o samtoolsflagstat_out/
samtools flagstat -@ 8 SAMPLE_bwameth.bam -o samtoolsflagstat_bwamethout/

#mergebams_technicalreplicates
samtools merge all.bam `find /basedir/ -name "*myfiles*.bam"`

#sortbam
for i in $(ls *.bam);do
samtools sort "$i" > "$i".sort
done ;

#rename
for i in *bam.sort
do
mv -- "$i" "${i/%.bam.sort/.sort.bam}"
done

#indexbam
for i in *.sort.bam
do samtools index "$i"
done 

#filterbam
for i in *sort.bam
do 
bamtools filter -isMapped true -isPaired true -isProperPair true -forceCompression -in $i -out $i.filter
done

#sort
for i in $(ls *.filter.bam);do
samtools sort "$i" > "$i".sort
done ;

#rename
for i in *sort.bam.filter
do
mv -- "$i" "${i/%.sort.bam.filter/_filter_sort.bam}"
done

#picard_mark_duplicates
java –Xmx4g –Xms4g -jar /path/to/Picard/MarkDuplicates.jar MarkDuplicates \
      I=input.bam \
      O=_marked_duplicates.bam \
      M=marked_dup_metrics.txt
      ASSUME_SORTED=true
      Index=true 

#index
for i in *filter_sort_marked_duplicates.bam
do samtools index "$i"
done

#flagstat 
samtools flagstat -@ 8 SAMPLE_filter_sort_marked_duplicates.bam -o samtoolsflagstat_out/

#Methylationbias_methyldackel
for i in *filter_sort_marked_duplicates.bam
do MethylDackel mbias <path_to_genome_folder> hg19.fa --minDepth 10 "$i" -o "$i"
done 

#ClipOverlap
for i in *filter_sort_marked_duplicates.bam
do bam clipOverlap --stats --in "$i" --out "$i".clipped.bam
done

#MethylDackel_extract 
for i in *filter_sort_marked_duplicates.bam
do MethylDackel extract <path_to_genome_folder> hg19.fa --minDepth 10 "$i" -o "$i"_methyldackelextract 
done

#MethylDackel_logit 
for i in *filter_sort_marked_duplicates.bam
do MethylDackel extract <path_to_genome_folder> hg19.fa --minDepth 10 --logit "$i" -o "$i"_methyldackel 
done

#MethylDackel_perRead 
for i in *filter_sort_marked_duplicates.bam
do MethylDackel extract <path_to_genome_folder> hg19.fa --minDepth 10 --perRead "$i" -o "$i"_methyldackel 
done

#epiread 
for i in *filter_sort_marked_duplicates.bam
do biscuit epiread <path_to_genome_folder> hg19.fa "$i" > "$i".epiread
done

#pileup
for i in *filter_sort_marked_duplicates.bam
do biscuit pileup -p <path_to_genome_folder> hg19.fa "$i" > "$i".vcf
done 

#betabed
for i in *.vcf
do biscuit vcf2bed -k 10 "$i" > "$i".beta.bed
done

#snpbed
for i in *.vcf
do biscuit vcf2bed -k 10 -t snp "$i" > "$i".snp.bed
done

#wgbs_tools_bam2pat
cd wgbs_tools
for i in *filter_sort_marked_duplicates.bam
do python3 wgbs_tools.py bam2pat --genome hg19 --out_dir bam2pat_out/ "$i"
done