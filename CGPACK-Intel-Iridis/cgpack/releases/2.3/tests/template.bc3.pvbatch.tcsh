#!/bin/tcsh
#
#$Id: template.bc3.pvbatch.tcsh 11 2014-12-01 10:07:20Z mexas $
#
#PBS -l walltime=00:10:00,nodes=1:ppn=16
#PBS -j oe
####PBS -m abe

cd $HOME/nobkp/cgpack/branches/coarray/tests
setenv OUTFILE z.out
#echo > $OUTFILE

module load apps/paraview-4.0.1
echo "module list" >> $OUTFILE
`module list` >> $OUTFILE 

######echo "LD_LIBRARY_PATH: " $LD_LIBRARY_PATH >> $OUTFILE
######echo "which mpirun: " `which mpirun` >> $OUTFILE
######export I_MPI_DAPL_PROVIDER=ofa-v2-ib0 >> $OUTFILE


mpirun -n 16 pvbatch --use-offscreen-rendering z.xdmf.py
#mpirun -n 16 pvbatch z.xdmf.py
