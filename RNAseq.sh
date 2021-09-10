for file in *.fastq
do
        hisat2 -p 16 --rna-strandness R -x genome -U $file -S AL_$file
done

for file in AL_*
do
  samtools view -S -b -q 30 -o Q30_$file $file
done

for file in Q30_*
do
  samtools sort -@ 16 -o Sort_$file $file
done

for file in Sort_*
do
bedtools genomecov -ibam $file -g hg19.chrom.sizes -bg -strand - -split > BGM_$file
        bedtools genomecov -ibam $file -g hg19.chrom.sizes -bg -strand + -split > BGP_$file
        sort -k1,1 -k2,2n BGM_$file
        sort -k1,1 -k2,2n BGP_$file
        ./bedGraphToBigWig BGM_$file hg19.chrom.sizes GTM_$file.bw
        ./bedGraphToBigWig BGP_$file hg19.chrom.sizes GTP_$file.bw
done

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


for file in SMinus_*
do
  htseq-count -s reverse -f bam $file --additional-attr=gene_name refseq.gtf > Counts_$file
done

for file in SPlus_*
do
  htseq-count -s reverse -f bam $file --additional-attr=gene_name refseq.gtf > Counts_$file
done
