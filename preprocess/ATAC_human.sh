#BSUB -q normal
#BSUB -J ATAC
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 54


# input reference and tools [[choose reference]]
reference=/share/home/hanxiaoping/tools/BWA_Reference_Human_hg19/hg19.fa
dropseq_root=/share/home/hanxiaoping/tools/Drop-seq_tools-2.5.1/
picard_jar=${dropseq_root}/3rdParty/picard/picard.jar
bwa_exec=/share/home/hanxiaoping/tools/bwa-0.7.15/bwa
correctBC_script=/share/home/hanxiaoping/tools_new/MWATAC/mw_ATAC_correct.py
barcodepath=/share/home/hanxiaoping/tools/Microwellseq_barcode/

sample_name=$(basename `pwd`)
outdir=bowtie_out
tmpdir=tmp

# create new files
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

if [ ! -d $tmpdir ]; then
    mkdir $tmpdir
fi


# fastq.gz --> bam
java -jar ${picard_jar} FastqToSam F1=H_R1.fq.gz F2=H_R2.fq.gz O=H.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name TMP_DIR=$tmpdir 

# tag R1&R2 with barcodes 
${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular.bam_summary_R1.txt \
    BASE_RANGE=1-18:25-34 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=true DISCARD_READ=false TAG_NAME=CB NUM_BASES_BELOW_QUALITY=5 \
    INPUT=H.bam OUTPUT=$tmpdir/unaligned_tagged_Cell_R1.bam COMPRESSION_LEVEL=0
${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular.bam_summary_R2.txt \
    BASE_RANGE=1-18:25-34 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=false DISCARD_READ=false TAG_NAME=CB NUM_BASES_BELOW_QUALITY=5 \
  	INPUT=$tmpdir/unaligned_tagged_Cell_R1.bam OUTPUT=$tmpdir/unaligned_tagged_Cell.bam COMPRESSION_LEVEL=0 

# filter bam
${dropseq_root}/FilterBam TAG_REJECT=XQ \
  INPUT=$tmpdir/unaligned_tagged_Cell.bam \
  OUTPUT=$tmpdir/unaligned_tagged_Cell_filtered.bam 

# barcode correct
python $correctBC_script $barcodepath $tmpdir/unaligned_tagged_Cell_filtered.bam $tmpdir/unaligned_tagged_Cell_corrected.bam 

# # add barcode to qname
samtools view $tmpdir/unaligned_tagged_Cell_corrected.bam -H > $tmpdir/unaligned_tagged_qname.sam
samtools view $tmpdir/unaligned_tagged_Cell_corrected.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); printf "%s:%s\n", td["CB"], $0 } } }' >> $tmpdir/unaligned_tagged_qname.sam 

# Sam to Fastq
java -Xmx100g -jar ${picard_jar} SamToFastq INPUT=$tmpdir/unaligned_tagged_Cell.snap.sam READ1_TRIM=53 FASTQ=$tmpdir/unaligned_R1.fastq  SECOND_END_FASTQ=$tmpdir/unaligned_R2.fastq

# Alignment
snaptools align-paired-end \
    --input-reference=${reference} \
    --input-fastq1=$tmpdir/unaligned_R1.fastq  \
    --input-fastq2=$tmpdir/unaligned_R2.fastq  \
    --output-bam=$outdir/Aligned.out.bam \
    --aligner=$(basename $bwa_exec) \
    --path-to-aligner=$(dirname $bwa_exec) \
    --min-cov=0 \
    --num-threads=20 \
    --if-sort=False \ 
    --tmp-folder=$tmpdir \
    --overwrite=TRUE

# shift
samtools sort $outdir/Aligned.out.bam -o $outdir/sort.bam
samtools index $outdir/sort.bam
alignmentSieve --numberOfProcessors 8 --ATACshift -b $outdir/sort.bam -o $shift/shift.bam
samtools sort -n $shift/shift.bam -o $shift/shift.sort.bam

# Pre-processing.
samtools view $shift/shift.sort.bam -H  |cut -f 2,3 |grep -v VN |grep -v bwa |grep -v samtools |awk '{printf"%s\t%s\n", substr($1,4,length($1)-2), substr($2,4,length($2)-2)}' >$outdir/h.chrom.sizes

# bam --> snap
snaptools snap-pre \
    --input-file=$outdir/Aligned.out.bam \
    --output-snap=$outdir/H.snap \
    --genome-name=hg19 
    --genome-size=$outdir/h.chrom.sizes \
    --overwrite=True \
    --min-mapq=30 \
    --min-flen=0 \
    --max-flen=1000 \
    --keep-chrm=True \
    --keep-single=True \
    --keep-secondary=True \
    --keep-discordant=True \
    --max-num=50000 \ 
    --min-cov=10 \
    --verbose=True

# snap to fragment file
python /share/home/hanxiaoping/tools_new/MWATAC/fraglist.py $outdir/H.snap 

# sort fragment file
sort -k1,1V -k2,2n -k3,3n barcode_fragment_human.bed > H.sort.bed
bgzip H.sort.bed



