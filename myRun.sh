# example3.pbs
# invoke with
# qsub -l nodes=2:ppn=1,pmem=1800mb,mem=1800mb,ncpus=2,cput=5:00:00 example3.pbs
# note that we are using one cpu per node. That :ppn=1 is critical.

# comment these out if you wish
# echo "qsub host is"
# echo $PBS_O_HOST
# echo "original queue is"
# echo $PBS_O_QUEUE
# echo "qsub working directory absolute is"
# echo $PBS_O_WORKDIR
# echo "pbs environment is"
# echo $PBS_ENVIRONMENT
# echo "pbs batch id"
# echo $PBS_JOBID
# echo "pbs job name from me is"
# echo $PBS_JOBNAME
# echo "Name of file containing nodes is"
# echo $PBS_NODEFILE
# echo "contents of nodefile is"
# cat $PBS_NODEFILE
# echo "Name of queue to which job went is"
# echo $PBS_QUEUE

# make sure we are in the right directory in case writing files
cd $PBS_O_WORKDIR

# Run the mpi job
# Note: the MPICH mpirun command is in /usr/local/bin
# use the full path
mpirun -np 2 -machine $PBS_NODEFILE ./mpicppHello