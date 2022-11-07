#BSUB -q normal
#BSUB -J 25-ATAC
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=36]"
#BSUB -n 36


# input reference and tools [[choose reference]]
reference=/share/home/hanxiaoping/tools/BWA-Reference-Human-Mouse-Merged/hg19_mm10_transgenes.fasta
dropseq_root=/share/home/hanxiaoping/tools/Drop-seq_tools-2.5.1/
picard_jar=${dropseq_root}/3rdParty/picard/picard.jar
bwa_exec=/share/home/hanxiaoping/tools/bwa-0.7.15/bwa
correctBC_script=/share/home/hanxiaoping/tools_new/MWATAC/mw_ATAC_correct.py
barcodepath=/share/home/hanxiaoping/tools/Microwellseq_barcode/

sample_name=$(basename `pwd`)
outdir=bowtie_out
tmpdir=tmp
shift=shift

# create new files
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

if [ ! -d $tmpdir ]; then
    mkdir $tmpdir
fi

if [ ! -d $shift ]; then
    mkdir $shift
fi


# fastq.gz --> bam
java -jar ${picard_jar} FastqToSam F1=H_R1.fq.gz F2=H_R2.fq.gz O=H.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name TMP_DIR=$tmpdir 

# tag R1&R2 with barcodes 
${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular.bam_summary_R1.txt \
   BASE_RANGE=1-18:19-28 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=true DISCARD_READ=false TAG_NAME=CB NUM_BASES_BELOW_QUALITY=5 \
   INPUT=H.bam OUTPUT=$tmpdir/unaligned_tagged_Cell_R1.bam COMPRESSION_LEVEL=0
${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular.bam_summary_R2.txt \
   BASE_RANGE=1-18:19-28 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=false DISCARD_READ=false TAG_NAME=CB NUM_BASES_BELOW_QUALITY=5 \
  	INPUT=$tmpdir/unaligned_tagged_Cell_R1.bam OUTPUT=$tmpdir/unaligned_tagged_Cell.bam COMPRESSION_LEVEL=0 

# filter bam
${dropseq_root}/FilterBam TAG_REJECT=XQ \
 INPUT=$tmpdir/unaligned_tagged_Cell.bam \
  OUTPUT=$tmpdir/unaligned_tagged_Cell_filtered.bam 

# barcode correct
python $correctBC_script $barcodepath $tmpdir/unaligned_tagged_Cell_filtered.bam $tmpdir/unaligned_tagged_Cell_corrected.bam 

# add barcode to qname
samtools view $tmpdir/unaligned_tagged_Cell_corrected.bam -H > $tmpdir/unaligned_tagged_qname.sam
samtools view $tmpdir/unaligned_tagged_Cell_corrected.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); printf "%s:%s\n", td["CB"], $0 } } }' >> $tmpdir/unaligned_tagged_qname.sam 

# Sam to Fastq
java -Xmx100g -jar ${picard_jar} SamToFastq INPUT=$tmpdir/unaligned_tagged_qname.sam  READ1_TRIM=28 FASTQ=$tmpdir/unaligned_R1.fastq  SECOND_END_FASTQ=$tmpdir/unaligned_R2.fastq

#Step 3. Alignment
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

# # filter
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/FilterBam I=$outdir/Aligned.out.bam O=$outdir/mouse.bam REF_SOFT_MATCHED_RETAINED=MOUSE
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/FilterBam I=$outdir/Aligned.out.bam O=$outdir/human.bam REF_SOFT_MATCHED_RETAINED=HUMAN

# shift
#human
samtools sort $outdir/human.bam -o $outdir/human_sort.bam
samtools index $outdir/human_sort.bam
alignmentSieve --numberOfProcessors 8 --ATACshift -b $outdir/human_sort.bam -o $shift/shift_human.bam
samtools sort -n $shift/shift_human.bam -o $shift/shift_human.sort.bam

#mouse
samtools sort $outdir/mouse.bam -o $outdir/mouse_sort.bam
samtools index $outdir/mouse_sort.bam
alignmentSieve --numberOfProcessors 8 --ATACshift -b  $outdir/mouse_sort.bam -o $shift/shift_mouse.bam
samtools sort -n $shift/shift_mouse.bam -o $shift/shift_mouse.sort.bam

# # Step 4. Pre-processing.
samtools view $outdir/human_sort.bam -H |cut -f 2,3 |grep -v VN |grep -v bwa |grep -v samtools |awk '{printf"%s\t%s\n", substr($1,4,length($1)-2), substr($2,4,length($2)-2)}' >$outdir/h.chrom.sizes
samtools view $outdir/mouse_sort.bam -H |cut -f 2,3 |grep -v VN |grep -v bwa |grep -v samtools |awk '{printf"%s\t%s\n", substr($1,4,length($1)-2), substr($2,4,length($2)-2)}' >$outdir/m.chrom.sizes

# mouse
snaptools snap-pre \
    --input-file=$shift/shift_mouse.sort.bam \
	--output-snap=$shift/mouse.snap \
	--genome-name=mm10 \
	--genome-size=$outdir/m.chrom.sizes \
	--overwrite=True \
	--min-mapq=30 \
	--min-flen=0 \
	--max-flen=1000 \
	--keep-chrm=True \
    --keep-single=False \
	--keep-secondary=False \
    --keep-discordant=False \
	--max-num=20000 \
	--min-cov=10 \
	--verbose=True

# human
snaptools snap-pre \
    --input-file=$shift/shift_human.sort.bam \
	--output-snap=$shift/human.snap \
	--genome-name=hg19 \
	--genome-size=$outdir/h.chrom.sizes \
	--min-mapq=30 \
	--min-flen=0 \
	--max-flen=1000 \
	--keep-chrm=True \
    --keep-single=False \
	--keep-secondary=False \
    --keep-discordant=False \
	--overwrite=True \
	--max-num=50000 \
	--min-cov=10 \
	--verbose=True

# snap to fragment file
python /share/home/hanxiaoping/tools_new/MWATAC/fraglist.py $shift/human.snap $shift/human.barcode_fragment.bed
python /share/home/hanxiaoping/tools_new/MWATAC/fraglist.py $shift/mouse.snap $shift/mouse.barcode_fragment.bed

# sort fragment file
sort -k1,1V -k2,2n -k3,3n $shift/human.barcode_fragment.bed > $shift/human.sort.bed
sort -k1,1V -k2,2n -k3,3n $shift/mouse.barcode_fragment.bed > $shift/mouse.sort.bed

# modify chr
zcat $shift/human.sort.bed | awk '{print substr($1,7,length($1)-6) "\t" $2 "\t" $3 "\t" $4}' > $shift/human.sort2.bed
zcat $shift/mouse.sort.bed | awk '{print substr($1,7,length($1)-6) "\t" $2 "\t" $3 "\t" $4}' > $shift/mouse.sort2.bed

# # bgzip 
bgzip $shift/human.sort2.bed
bgzip $shift/mouse.sort2.bed

# frag_length
samtools view $shift/shift_human.sort.bam | awk '{if ($9 > 0){print substr($1,1,28) "\t" $9}}' > $shift/frag_length.txt
cat $shift/frag_length.txt | grep -f $shift/bc.txt > frag_length.txt

# rm tmpdir if final output exist
if [ -f shift/H.snap ]; then
    rm -r $tmpdir
fi

