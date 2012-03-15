#!/bin/bash

# SOPRA manual and soures are available here:
# http://www.physics.rutgers.edu/~anirvans/SOPRA/SOPRA_v1.4.5.zip

if [ ! -e reads1.fastq ] || [ ! -e reads2.fastq ] || [ ! -e contigs.fasta ] || [ $# -ne 1 ]
then
  echo "Usage: put files reads1.fastq, reads2.fastq, and contigs.fasta into the current directory and run $0 <insertsize>"
  echo "Bye!"
  exit
fi

INSERTSIZE=$1

READS1=reads1.fastq
READS2=reads2.fastq
CONTIGS=contigs
OUTDIR=outdir
C=5 # If the number of times a read and its reverse complement appear in the library is equal to or more than this value, the pairing information from that read will be disregarded
# E=0 # if set equal to 1, the empirical value for the insert size will not be used
# W=1 # Minimum number of links between two contigs
# H=20.2 # High coverage contigs (above mean coverage + h×std coverage) are not considered in the scaffold assembly mainly to exclude reads from repetitive regions


# !!!HACK!!! to be replaced with a check whether the directory already exists
rm -rf ${OUTDIR}

# merge paired reads into one file
perl shuffleSequences_fastq.pl ${READS1} ${READS2} reads.fastq

# convert reads from fastq to fasta
# WARNING: THIS SHOULD BE REPLACED WITH A MORE ACCURATE SCRIPT AS @ SIGN MAY BE PRESENT AS A QUALITY VALUE ALSO
cat reads.fastq | perl -e '$i=0;while(<>){if(/^\@/&&$i==0){s/^\@/\>/;print;}elsif($i==1){print;$i=-3}$i++;}' > reads.fasta

# STEP 1: prepare contigs and paired reads
perl s_prep_contigAseq_v1.4.5.pl -contig ${CONTIGS}.fasta -mate reads.fasta -a ${OUTDIR}

# align reads to contigs
cd ${OUTDIR}
bwa index -a is ${CONTIGS}_sopra.fasta
bwa aln ${CONTIGS}_sopra.fasta reads_sopra.fasta > readstocontigs.sai
bwa samse ${CONTIGS}_sopra.fasta readstocontigs.sai reads_sopra.fasta > readstocontigs.sam
cd ..

# STEP 2: remove the pairs where at least one of the reads either do not map or maps in multiple places
perl s_parse_sam_v1.4.5.pl -sam ${OUTDIR}/readstocontigs.sam -a ${OUTDIR}

# STEP 3
perl s_read_parsed_sam_v1.4.5.pl -parsed ${OUTDIR}/readstocontigs.sam_parsed -d ${INSERTSIZE} -a ${OUTDIR} -c ${C} # -e ${E}

# STEP 4: scaffold assembly
perl s_scaf_v1.4.5.pl -o ${OUTDIR}/orientdistinfo_c${C} -a ${OUTDIR} # -w ${W} -h ${H}

cp ${OUTDIR}/scaffold* .

