for R1 in *1.fq.gz
do
sample=${R1//_CKDL*}
R2=${R1//1.fq/2.fq}
mkdir trimmed_fastqs
trimmomatic SE -threads 3 $R1 ${sample}_L001_R1_001.fastq.gz CROP:90
trimmomatic SE -threads 3 $R2 ${sample}_L001_R2_001.fastq.gz CROP:90
mv *001.fastq.gz trimmed_fastqs/
done