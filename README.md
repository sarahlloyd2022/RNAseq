## RNA-seq Processing in terminal

**1. Align fastq files with HISAT2**  
   This loop will run hisat2 for all files ending in .fastq  
   genome should be built ahead of time with bowtie2-build and placed in working directory  
      ie. run bowtie2-build on hg19, where bt2_index_base is 'genome'  
      Usage: bowtie2-build [options]* <reference_in> <bt2_index_base>  
```
for file in *.fastq
do
        hisat2 -p 16 --rna-strandness R -x genome -U $file -S AL_$file
done
```
**2. Convert from Sam to Bam format and filter by MAPQ score**  
   Sam files are much larger, human readable files  
   Bam files are binary and easier to work with and store  
   -q option sets MAPQ alignment score cut off  
   
```
for file in AL_*
do
  samtools view -S -b -q 30 -o Q30_$file $file
done
```
**3. Sort bam files**

```
for file in Q30_*
do
  samtools sort -@ 16 -o Sort_$file $file
done
```
**4. Generate bigwig files to make genome browser tracks with**  
   This is an optional step only necessary if browser tracks will be made  
   Bigwig files can be stored on google drive and visualized with ucsc genome browser
   
```
for file in Sort_*
do
bedtools genomecov -ibam $file -g hg19.chrom.sizes -bg -strand - -split > BGM_$file
        bedtools genomecov -ibam $file -g hg19.chrom.sizes -bg -strand + -split > BGP_$file
        sort -k1,1 -k2,2n BGM_$file
        sort -k1,1 -k2,2n BGP_$file
        ./bedGraphToBigWig BGM_$file hg19.chrom.sizes GTM_$file.bw
        ./bedGraphToBigWig BGP_$file hg19.chrom.sizes GTP_$file.bw
done
```

**5. Separate strands and re-sort**

```
for file in Sort_*
do
  samtools view -bh -F 16 $file > Plus_$file
  samtools view -bh -f 16 $file > Minus_$file
done
for file in Plus_*
do
  samtools sort -@ 16 -o S$file $file
done
for file in Minus_*
do
  samtools sort -@ 16 -o S$file $file
done
```

**6. Count total number of reads aligning at each gene**
    Refseq files for genomes are publically available for download
```
for file in SMinus_*
do
  htseq-count -s reverse -f bam $file --additional-attr=gene_name refseq.gtf > Counts_$file
done

for file in SPlus_*
do
  htseq-count -s reverse -f bam $file --additional-attr=gene_name refseq.gtf > Counts_$file
done
```

**7. Paste counts tables together as desired**  
  some extra columns will need to be removed  
  final table should be column 1 as gene names and rest of columns as total counts/gene  

```
paste ControlRep1 ControlRep2 TreatmentRep1 TreatmentRep2 > Final Table
```

#### Differential Expression analysis can now be completed in R with software such as DESeq2


