#BSUB -q normal
#BSUB -J mouse
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 36

cd `pwd`
mkdir tmp
tmpdir=tmp
sample_name=$(basename `pwd`)
## fastq --> bam
java -jar /share/home/hanxiaoping/tools/Drop-seq_tools-1.12/3rdParty/picard/picard.jar FastqToSam F1=H_R1.fq.gz F2=H_R2.fq.gz  O=H.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name TMP_DIR=$tmpdir && rm H_R1.fq.gz && rm H_R2.fq.gz

## Tag Cell Barcode:
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/TagBamWithReadSequenceExtended INPUT=H.bam OUTPUT=unaligned_tagged_Cell.bam  SUMMARY=unaligned_tagged_Cellular.bam_summary.txt BASE_RANGE=1-18 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1

## Tag Molecular Barcode:
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/TagBamWithReadSequenceExtended INPUT=unaligned_tagged_Cell.bam OUTPUT=unaligned_tagged_CellMolecular.bam SUMMARY=unaligned_tagged_Molecular.bam_summary.txt BASE_RANGE=19-23 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1

## FilterBAM:
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/FilterBAM TAG_REJECT=XQ INPUT=unaligned_tagged_CellMolecular.bam OUTPUT=unaligned_tagged_filtered.bam

## PolyATrimmer
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/PolyATrimmer INPUT=unaligned_tagged_filtered.bam OUTPUT=unaligned_mc_tagged_polyA_filtered.bam OUTPUT_SUMMARY=polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=10

## corrected BC for one mismatch
barcodepath=/share/home/hanxiaoping/tools/Microwellseq_barcode/

python $barcodepath/mw_RNA_correct.py $barcodepath unaligned_mc_tagged_polyA_filtered.bam filtered.bam
echo "correct sam files done"

## bam --> fastq
java -Xmx100g -jar /share/home/hanxiaoping/tools/Drop-seq_tools-1.12/3rdParty/picard/picard.jar SamToFastq INPUT=filtered.bam FASTQ=unaligned_mc_tagged_polyA_filtered.fastq

## Alignment STAR
/share/home/hanxiaoping/tools/STAR-2.5.2a/source/STAR --genomeDir /share/home/hanxiaoping/tools/STAR_Reference_Human/ --readFilesIn unaligned_mc_tagged_polyA_filtered.fastq --outFileNamePrefix star

## SortSam
java -Xmx100g -jar /share/home/hanxiaoping/tools/Drop-seq_tools-1.12/3rdParty/picard/picard.jar SortSam I=starAligned.out.sam O=aligned.sorted.bam SO=queryname

## MergeBamAlignment
java -Xmx100g -jar /share/home/hanxiaoping/tools/Drop-seq_tools-1.12/3rdParty/picard/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=/share/home/hanxiaoping/tools/STAR_Reference_Human/Homo_sapiens.GRCh38.fa UNMAPPED_BAM=filtered.bam ALIGNED_BAM=aligned.sorted.bam OUTPUT=merged.bam INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false

## TagReadWithGeneExon
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/TagReadWithGeneExon I=merged.bam O=star_gene_exon_tagged.bam ANNOTATIONS_FILE=/share/home/hanxiaoping/tools/STAR_Reference_Human/Homo_sapiens.GRCh38.fa.gtf TAG=GE

## Digital Gene Expression
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/DigitalExpression I=star_gene_exon_tagged.bam O=_dge.txt.gz SUMMARY=_dge.summary.txt NUM_CORE_BARCODES=30000

files_to_delete="unaligned_tagged_Cell.bam unaligned_tagged_CellMolecular.bam unaligned_tagged_filtered.bam unaligned_mc_tagged_polyA_filtered.fastq merged.bam aligned.sorted.bam unaligned_mc_tagged_polyA_filtered.bam filtered.bam starAligned.out.sam"

# rm tmpdir if final output exist
if [ -f ./_dge.summary.txt ]; then
    rm $files_to_delete
fi



