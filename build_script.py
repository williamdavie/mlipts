
'''


Maximum efficiency is to batch every single job individually.

On archer2:

    Limit = 64 jobs and 16 running.
    
    Example dataset size = 10,000 DFT calculation
    
    10,000/64 = 157 submissions (way to many)
    
    So we may want to reduce this to say 16 submissions node repeats 10 DFT calculations.
    
    Therefore we define a partitioning system. 
    
    1 partition = 1 submission to sbatch. 
    
To maintain some kind of simplicity, partitioning is only done for:

main
  |
  |-base_dir
  |-lammps calculations
            |-calc_i
                |-DFT caculations
    
    
Example:

    5 lammps cal: 1 partion
        - 20 DFT calculations per = 100 total: 
            - quickest route: 50 partitions: 4 hours per
    
    

- Run LAMMPS calculations - (partition part 1)

- Run DFT calculations 

'''
