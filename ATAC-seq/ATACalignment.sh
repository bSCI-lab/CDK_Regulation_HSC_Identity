#####ATAC alignment
config="/List.txt"
# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
pair1=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
pair2=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

#module in use
module load bowtie2/2.3.5.1
module load BEDTools/2.31.1-GCC-11.2.0
module load SAMtools/1.19-foss-2021b
module load deepTools/3.5.1-foss-2021b

#strip suffix to give filename
f=$(basename "$pair1")
FILENAME="${f%%_1.fq.gz}"

# Create SAM file name from sample name
SAMFILE="${FILENAME}.sam"
BAMFILE="${FILENAME}_tsrmChM_bkl.bam"

bowtie2 -p 8 --very-sensitive --local -X 1000 --no-mixed \
  -x /Reference/hg38_bowtie/grch38_1kgmaj \
  -1 "$pair1" -2 "$pair2" | \
samtools view -@ 8 -bS - | \
samtools sort -@ 8 -n - | \
samtools fixmate -m - - | \
samtools sort -@ 8 - | \
samtools markdup - "${FILENAME}.bam"
samtools index "${FILENAME}.bam"

samtools flagstat "${FILENAME}.bam" > "${FILENAME}_flagstat.txt"
#ATAC shift
alignmentSieve --ATACshift --bam ${FILENAME}.bam -o ${FILENAME}_tsNoAl.bam
samtools sort -O bam -o ${FILENAME}_ts.bam ${FILENAME}_tsNoAl.bam
samtools index ${FILENAME}_ts.bam
rm ${FILENAME}.bam
rm ${FILENAME}.bam.bai
rm ${FILENAME}_tsNoAl.bam
#remove mytochondrial DNA, intersect with blacklist to take it out
samtools view -h ${FILENAME}_ts.bam | awk 'substr($0,1,1) == "@" || (length($3) <= 5 && $3 != "chrM" && $3 != "*") {print $0}' | samtools view -bS - > ${FILENAME}_tsrmChM.bam
samtools index ${FILENAME}_tsrmChM.bam
bedtools intersect -abam ${FILENAME}_tsrmChM.bam -b /Reference/hg38-blacklist.v2.bed -v > $BAMFILE
samtools index $BAMFILE
rm ${FILENAME}_tsrmChM.bam.bai
rm ${FILENAME}_tsrmChM.bam
rm ${FILENAME}_ts.bam
rm ${FILENAME}_ts.bam.bai
#generate bigwig to upload on UCSC browser
bamCoverage --bam $BAMFILE --binSize 5  --normalizeUsing RPGC  --effectiveGenomeSize 2913022398 --maxFragmentLength 1000 -o ${FILENAME}.bw --ignoreForNormalization chrX --ignoreDuplicates
