#!/bin/bash
#SBATCH --account=bgmp               ###
#SBATCH --partition=bgmp       ### Partition
#SBATCH --output=deduper%j.out
#SBATCH --nodes=1              ### Number of Nodes
#SBATCH --mail-type=END              ### Mail events (NONE, BEGIN, END, FA$
#SBATCH --mail-user=sobermil@uoregon.edu  ### Where to send mail
#SBATCH --cpus-per-task=1
#SBATCH --error=deduper%j.err

conda activate bgmp_py310
/usr/bin/time -v /projects/bgmp/sobermil/bioinfo/Bi624/deduper/obermiller_deduper.py \
-f /projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam -u /projects/bgmp/sobermil/bioinfo/Bi624/deduper/STL96.txt -o \
/projects/bgmp/sobermil/bioinfo/Bi624/deduper/C1_SE_uniqAlign_output.sam