

# Machine Learned Interatomic Potentials - Training Suite (MLIPTS).

Seemless data collection for training/fine-tuning machine learned interatomic potentials. 

The key idea behind MLIPTS, is to perform the following active learning workflow with as little user input as possible:

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
