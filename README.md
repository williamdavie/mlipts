

# Machine Learned Interatomic Potentials - Training Suite.
<img src="https://img.shields.io/badge/version-1.0.0-blue"> <img src="https://img.shields.io/badge/contributors-welcome-green">


MLIPTS is a python package for training/fine-tuning machine learned interatomic potentials. 

The key idea is to perform the following active learning workflow with as little user input as possible:

<p align="center">
  <img src="https://github.com/williamdavie/mlipts/blob/docs-edit/media/active_learning_flowchart.png" width="50%" height="auto">
</p>

> [!NOTE]
> The scope of a fully-fledged python package to perform seemless MLIP training is significant, with many availible MD and DFT/Quantum Chemistry codes and ways to quantify data quality.
> Contributors are welcome to help make this goal a reality.

## Version 1.0.0

### Installation

MLIPTS can be installed via pip

``` 
pip install mlipts
 ```

### Capability

```1.0.0``` is built to address the creation of an inital data set (Part 1 of the workflow above), however, note part of the current functionality is applicable across the main workflow, ultimately arriving at the following reduced workflow to address:

<p align="center">
  <img src="https://github.com/williamdavie/mlipts/blob/docs-edit/media/workflow_%231.png" width="50%" height="auto">
</p>

The MD code supported is ```LAMMPS``` and DFT code supported is ```VASP```, where the earth movers distance (EMD) has been implemented to filter configurations. The details of this method are found at [1] and [average-minimum-distance](https://github.com/dwiddo/average-minimum-distance) (Copyright (C) 2025 Daniel Widdowson).

### Quickstart

It is highly recommended to follow the availible example at (insert link to example). 

The working directory is set up in the following way:
```
collect_data/
├─ MD_base/
├─ QM_base/
└─ workflow.ipynb
```

Where ```MD_base``` and ```QM_base``` include the input files for molecular dynamics and quantum mechanical simulations respectively. Since the current version only supports lammps and vasp, these directories will have the following format:

```
├─ lammps_base/
    ├─ in.test
    └─ test.dat
├─ vasp_base/
    ├─ INCAR
    ├─ KPOINTS
    └─ POTCAR
```

Noting ```POSCAR``` is intentially missing as this is to be generated. The key to a successful data collection workflow is ensuring all files in the above are formatted correctly, as to collect a full dataset we will be calling each calculation many times. With a directory set up, mlipts allows simply following of the flow chart above: 

1. Run many MD calculations:
```python
workflow.build_MD_calculations('./lammps_base',variables,outdir='./MD_calculations')
workflow.write_MD_submission_scripts(MD_cmd_line,submit=True)
```
2. Filter new configurations from MD:
```python
workflow.
```
3. Run DFT calculations on new configurations
```python
workflow.
```
4. Save to a training data set. 
```python
workflow.
```



The final workflow will then appear as:
```
collect_data/
├─ MD_base/
├─ MD_calculations/
├─ MD_scripts/
├─ QM_base/
├─ QM_calculations/
├─ QM_scripts/
├─ workflow.ipynb
└─ training_data.xyz
```









