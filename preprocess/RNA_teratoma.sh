#BSUB -q normal
#BSUB -J mix
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=10]"
#BSUB -n 10

cd `pwd`
mkdir tmp
tmpdir=tmp
sample_name=$(basename `pwd`)
# fastq --> bam
java -jar /share/home/hanxiaoping/tools/Drop-seq_tools-1.12/3rdParty/picard/picard.jar FastqToSam F1=H_R1.fq.gz F2=H_R2.fq.gz  O=H.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name TMP_DIR=$tmpdir

## Tag Cell Barcode:
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/TagBamWithReadSequenceExtended INPUT=H.bam OUTPUT=unaligned_tagged_Cell.bam  SUMMARY=unaligned_tagged_Cellular.bam_summary.txt BASE_RANGE=1-18 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1

## Tag Molecular Barcode:
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/TagBamWithReadSequenceExtended INPUT=unaligned_tagged_Cell.bam OUTPUT=unaligned_tagged_CellMolecular.bam SUMMARY=unaligned_tagged_Molecular.bam_summary.txt BASE_RANGE=19-23 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1

## FilterBAM:
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/FilterBAM TAG_REJECT=XQ INPUT=unaligned_tagged_CellMolecular.bam OUTPUT=unaligned_tagged_filtered.bam

## PolyATrimmer
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/PolyATrimmer INPUT=unaligned_tagged_filtered.bam OUTPUT=unaligned_mc_tagged_polyA_filtered.bam OUTPUT_SUMMARY=polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=10

## corrected bam for one mismatch
barcodepath=/share/home/hanxiaoping/tools/Microwellseq_barcode/

python $barcodepath/mw_RNA_correct.py $barcodepath unaligned_mc_tagged_polyA_filtered.bam filtered.bam
echo "correct sam files done"

## bam --> fastq
java -Xmx100g -jar /share/home/hanxiaoping/tools/Drop-seq_tools-1.12/3rdParty/picard/picard.jar SamToFastq INPUT=filtered.bam FASTQ=unaligned_mc_tagged_polyA_filtered.fastq

## Alignment STAR
/share/home/hanxiaoping/tools/STAR-2.5.2a/source/STAR --genomeDir /share/home/hanxiaoping/tools/Human-Mouse-Merged/  --readFilesIn unaligned_mc_tagged_polyA_filtered.fastq --outFileNamePrefix star

## SortSam
java -Xmx100g -jar /share/home/hanxiaoping/tools/Drop-seq_tools-1.12/3rdParty/picard/picard.jar SortSam I=starAligned.out.sam O=aligned.sorted.bam SO=queryname

## MergeBamAlignment
java -Xmx100g -jar /share/home/hanxiaoping/tools/Drop-seq_tools-1.12/3rdParty/picard/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=/share/home/hanxiaoping/tools/Human-Mouse-Merged/hg19_mm10_transgenes.fasta UNMAPPED_BAM=filtered.bam ALIGNED_BAM=aligned.sorted.bam OUTPUT=merged.bam INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false

## TagReadWithGeneExon
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/TagReadWithGeneExon I=merged.bam O=star_gene_exon_tagged.bam ANNOTATIONS_FILE=/share/home/hanxiaoping/tools/Human-Mouse-Merged/hg19_mm10_transgenes.gtf TAG=GE

## FilterBAM
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/FilterBAM  I=./star_gene_exon_tagged.bam O=./mouse.bam REF_SOFT_MATCHED_RETAINED=MOUSE
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/FilterBAM  I=./star_gene_exon_tagged.bam O=./human.bam REF_SOFT_MATCHED_RETAINED=HUMAN

## DigitalExpression
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/DigitalExpression I=./mouse.bam O=./mouse.dge.txt.gz  SUMMARY=mouse.dge.summary.txt NUM_CORE_BARCODES=20000
/share/home/hanxiaoping/tools/Drop-seq_tools-1.12/DigitalExpression I=./human.bam O=./human.dge.txt.gz  SUMMARY=human.dge.summary.txt NUM_CORE_BARCODES=20000

files_to_delete="unaligned_tagged_Cell.bam unaligned_tagged_CellMolecular.bam unaligned_tagged_filtered.bam unaligned_mc_tagged_polyA_filtered.fastq merged.bam aligned.sorted.bam unaligned_mc_tagged_polyA_filtered.bam filtered.bam starAligned.out.sam"

# rm tmpdir if final output exist
if [ -f ./human.dge.summary.txt ]; then
    rm $files_to_delete
fi



