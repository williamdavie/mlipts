'''
@author: William Davie

File containing vasp specific functionality. Used to build many vasp calculations.

some of these functions may be generalised if other codes are added e.g. write_vasp_submission_scipts
'''

from ase import Atoms
from ase.io import read, write
import numpy as np
import shutil


def build_vasp_calculation(vasp_base_dir: str, config: Atoms, calc_name: str, outdir: str) -> str: 
    '''
    Builds a vasp calculation directory for a given atomic configuration.
    
    Parameters
    ----------
    vasp_base_dir: str
        path to base directory, should contain POTCAR, KPOINTS and INCAR.
    config: :class:`ase.Atoms` 
        atomic configuration.
    outname: str
        name of calculation directory
    outdir: str
        output path of calculation directory.
        
    Returns
    -------
    new_calc_dir: str 
        vasp directory generated in outdir.
    '''
    
    poscar = write_POSCAR_str(config)
    new_calc_dir = outdir + '/' + f'{calc_name}'
    shutil.copytree(vasp_base_dir, new_calc_dir, dirs_exist_ok=True)
            
    with open(new_calc_dir +'/POSCAR','w') as f:
                f.write(poscar)
    
    return new_calc_dir



def write_vasp_submission_scripts(calc_dirs: list[str],vasp_cmd_line: str, 
                                  num_script_partitions: int,
                                  scipts_outdir: str='./QMsubmission_scripts',
                                  save_and_remove: bool=True,
                                  database_file: str=None, 
                                  pythonenv: str=None):
    '''
    Writes a number of submission scripts to run vasp calculations
    
    Parameters
    ----------
    calc_dirs : list[str]
        A list of directories containing a vasp calculation (POSCAR,INCAR,KPOINTS,POTCAR)
    num_script_partitions: int
        number of submission scripts to generate, calculations per script = len(calc_dirs) / num_script_partitions
    save_and_remove: bool
        If true, after each calculation is complete, the relevant configuration is saved and all remaining data is removed (used to save disk space). Default is True.
    pythonenv: bool
        If save_and_remove == True, a pythonenv must be specified.

    Returns
    -------
    None : None
        generates a number of submission scripts.
    
    '''
    
    print(f'Generating scripts for {len(calc_dirs)} vasp calculations.')
    
        
        # Leaving this as is for now
        assert self.DFTcount % num_script_partitions == 0, f'Number DFT counts ({self.DFTcount}) must be divisible by number of partions ({num_script_partitions})'
        num_calcs_per_submission = self.DFTcount / num_script_partitions

        '''
        Since individual DFT calculations can take a while, it is a lot faster to partition calculations over a number of submission scripts. This is done via: num_script_partitions.
        '''
        
        print(f'Each script will carry out {num_calcs_per_submission} DFT caculations, over a time of {time_per_script}')
        
        
        if save_and_remove == True:
            assert database_file != None, "Error you have save data as your calculations are performed but did not specify a database_file name"
            assert pythonenv != None, 'Error you have save data, this requires python, please define a python enviroment.'
            # relies on a lot of formatting across the code
            savedata_cmd = f'{pythonenv}/bin/python -m mlipts.append_to_database $i {database_file}'
            remove_cmd = 'rm -r $i'
            
        else:
            savedata_cmd = ''
            remove_cmd = ''
        
        for i in range(num_script_partitions):
            
            vasp_submit = ''
        
            if hpc == 'Archer2':
            
                vasp_submit += archer2_submission_template(nodes=nodes,ranks=ranks,time=time_per_script)
        
            vasp_submit+='\n'
            
            current_dirs = ''
            for dir in self.electronic_calculation_dirs[int(i*num_calcs_per_submission):int((i+1)*num_calcs_per_submission)]:
                current_dirs += f'{dir} '
            
            vasp_submit += f'directories="{current_dirs}"\n'

            vasp_submit+=f'''for i in $directories; 
do 
cd $i
    {vasp_cmd_line}
cd -
{savedata_cmd}
{remove_cmd}
done\n'''
            
            if show==True:
                print('SCRIPT PREVIEW: ')
                print(vasp_submit)
            
            Path('./QM_submission').mkdir(exist_ok=True)
            
            with open(f'./QM_submission/submit_vasp_#{i}','w') as f:
                f.write(vasp_submit)
                self.QM_scripts.append(f'./QM_submission/submit_vasp_#{i}')
                
                print(f'Submission script saved to: ./QM_submission/submit_vasp_#{i}')



def write_POSCAR_str(config: Atoms) -> str:
    '''
    writes a POSCAR string given an atomic configuration.
    '''
    
    poscar = 'System\n 1.0\n'
    
    cell = np.array(config.cell)
    poscar += f' {cell[0,0]} {cell[0,1]} {cell[0,2]}\n {cell[1,0]} {cell[1,1]} {cell[1,2]}\n {cell[2,0]} {cell[2,1]} {cell[2,2]}\n'
    
    type_labels = config.symbols.species()
    for type in type_labels:
        poscar += f' {type} '
    poscar+='\n'

    for type in type_labels:
        count = config.symbols.count(type)
        poscar += f' {count} '
    poscar+='\nDirect\n'
    

    for pos in config.positions:
        poscar+=f'{pos[0]} {pos[1]} {pos[2]}\n'
 
    return poscar
             