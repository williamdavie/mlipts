

import re
def archer2_submission_template(nodes: int, ranks: int, time: str, account='e89-camm'):
    
    assert re.match(r'^\d{2}:\d{2}:\d{2}$', time), 'Time must have format XX:XX:XX'
    
    hours, minutes, seconds = time.split(":")
    
    if int(minutes) <= 20 and int(hours) == 0:
        
        qos='short'

    else:
        
        qos = 'standard'
    
    
    return f'''#!/bin/bash

#SBATCH --job-name=lammps_test
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={ranks}
#SBATCH --cpus-per-task=1
#SBATCH --time={time}

#SBATCH --account={account}
#SBATCH --partition=standard
#SBATCH --qos={qos}
    '''