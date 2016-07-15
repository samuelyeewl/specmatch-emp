#!/bin/bash

#PBS -N specmatchemp
#PBS -q default
#PBS -l nodes=24
#PBS -l walltime=12:00:00
#PBS -V
#PBS -m bae -M syee@caltech.edu

#change the working directory (default is home directory)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

# Write out some information on the job
echo Running on host `hostname`
echo Time is `date`

### Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS cpus

# Tell me which nodes it is run on
echo ” ”
echo This jobs runs on the following processors:
echo `cat $PBS_NODEFILE`
echo ” “


parallel --slf $PBS_NODEFILE -a script.txt
