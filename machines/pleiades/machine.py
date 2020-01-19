machine = {
        'submit_command': 'qsub -q long',
        'submit_script':
        '''#PBS -S /bin/bash
#PBS -l select=1:ncpus=24:mpiprocs=1:model=has
#PBS -N {job_name}
#PBS -l walltime=120:00:00
#PBS -q long
#PBS -W group_list=s1994
cd $PBS_O_WORKDIR
{command} > {output}
'''
}
