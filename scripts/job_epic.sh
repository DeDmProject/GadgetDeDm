#/bin/bash
#PBS -W group_list=astronomy167
#PBS -l nodes=64:ppn=8,mem=1280gb
#PBS -l walltime=11:59:00
#PBS -q routequeue

. /etc/profile.d/modules.sh
cd /home/aknebe/gadget_exec

module load intel-compilers
module load fftw/2.1.5-intel-11.1-mpi
module load gsl/1.14
module load mpi/openmpi-1.4.3-11-intel-11.1

mpirun ./P-Gadget2-dedm1024 parameterfiles/1024/cde099.param 
