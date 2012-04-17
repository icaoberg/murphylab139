#!/bin/bash
#
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -o job.stdout
#$ -e job.stderr
cd $PBS_O_WORKDIR
 
perl yeast_image_download.pl
