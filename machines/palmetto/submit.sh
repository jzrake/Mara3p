#PBS -N example
#PBS -M username@clemson.edu
#PBS -o example.o
#PBS -e example.e
#PBS -l select=1:ncpus=40:phase=19a,walltime=01:00:00

cd $PBS_O_WORKDIR
./mara opt1=val1 opt2=val2 > example.out
