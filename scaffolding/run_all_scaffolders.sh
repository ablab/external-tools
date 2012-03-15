#!/bin/bash

if [ ! -e reads1.fastq ] || [ ! -e reads2.fastq ] || [ ! -e contigs.fasta ] || [ $# -ne 1 ]
then
  echo "Usage: put files reads1.fastq, reads2.fastq, and contigs.fasta into the current directory and run $0 <insertsize>"
  echo "Bye!"
  exit
fi

WD=`pwd`

for myfile in reads1.fastq reads2.fastq contigs.fasta
do
  for mydir in sopra sspace/scaffolding opera/bin
  do
    ln -sf ${WD}/${myfile} ./${mydir}/${myfile}
  done
done

echo "Running SOPRA..."
cd ./sopra/
sh ./scaffolding.sh $1 &> log_sopra
cd ..
cp -f ./sopra/scaffolds_h2.2_L150_w4.fasta scaffolds_sopra.fasta
cp -f ./sopra/log_sopra log_sopra
echo "...OK"

echo "Running SSPACE..."
cd ./sspace/scaffolding/
sh ./scaffolding_sspace.sh $1 &> log_sspace
cd ../..
cp -f ./sspace/scaffolding/scaffolds.final.scaffolds.fasta scaffolds_sspace.fasta
cp -f ./sspace/scaffolding/log_sspace log_sspace
echo "...OK"


echo "Running OPERA..."
cd ./opera/bin/
sh ./scaffolding_opera.sh &> log_opera
cd ../..
cp -f ./opera/bin/scaffolds.fasta scaffolds_opera.fasta
cp -f ./opera/bin/log_opera log_opera
echo "...OK"
