machine = {
        'submit_command': 'qsub -V',
        'submit_script':
        '''#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -N {job_name}
cd $PBS_O_WORKDIR
{command} > {output}
'''
}
