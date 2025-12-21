
# Contributing

Please first contact me via: willdavie2002@gmail.com. There are many areas for improvement, some possible routes are outlined here.

## Version 0.1.X

### mlipt.codes

Main focus on perfecting support for LAMMPS and VASP, but aim to add support for other codes too.

|  job                                 | priority | description |
|--------------------------------------|----------|----------|
| vary LAMMPS parameters across sample space | High      |    need flexability on defining a sample space, e.g. defining a different random seed for velocities   |  
| Native VASP directory construction   | Med     |  Using information in the MD base directory define a simple vasp calculation setup.    |
| Native LAMMPS directory construction | Low     |   Given some basic information, construct a lammps directory    |  

### mlpit.similarity

|  job                                 | priority | description |
|--------------------------------------|----------|----------|
| Add additional descriptors to access similarity of configurations | Med      | E.g. SOAP, High priority if MLIPTS is going to be used for non-periodic solids. |  
| Optimize intensive calculations | Low     |  Unless analysing a huge dataset, this is not a bottleneck in the workflow  |

### mlipts.hpc_submission

|  job                                 | priority | description |
|--------------------------------------|----------|----------|
| Add further support for custom hpcs | Med |  |
 
### mlipts.append_to_database

|  job                                 | priority | description |
|--------------------------------------|----------|----------|
| Catch errors from QM calculations before appending | High | E.g. if calculations fail to converge.  |
| Add options to save additional data from QM calculation | Med | e.g. magnetic and charge states, (+ any useful output parameter of choice with a keyword.) |


### mlipts.data_collection

Support for all changes in the main class. 

## Version 0.2.X

Development of the active learning workflow. Native support for model training.


