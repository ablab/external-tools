#!/bin/bash

# SSPACE manual and soures are available here:
# http://www.baseclear.com/landingpages/sspacev11/

if [ ! -e reads1.fastq ] || [ ! -e reads2.fastq ] || [ ! -e contigs.fasta ] || [ $# -ne 1 ]
then
  echo "Usage: put files reads1.fastq, reads2.fastq, and contigs.fasta into the current directory and run $0 <insertsize>"
  echo "Bye!"
  exit
fi

LIBRARYINSERTSIZE=$1

READS1=reads1.fastq
READS2=reads2.fastq
CONTIGS=contigs.fasta
OUTFILEPREFIX=scaffolds
MINIMUMALLOWEDERROR=0.1
EXTENSION=0
PRINTDOT=1

# vars
LIBRARYFILENAME=library.txt

echo "lib1 ${READS1} ${READS2} ${LIBRARYINSERTSIZE} ${MINIMUMALLOWEDERROR} 0" > ${LIBRARYFILENAME}

perl ../SSPACE_v1-1.pl -l ${LIBRARYFILENAME} -s ${CONTIGS} -k 5 -a 0.7 -x ${EXTENSION} -p ${PRINTDOT} -b ${OUTFILEPREFIX}




