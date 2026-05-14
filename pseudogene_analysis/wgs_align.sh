bwa mem genome.fa reads_1.fastq reads_2.fastq > aln.sam
samtools sort -o aln.bam aln.sam
samtools index aln.bam
