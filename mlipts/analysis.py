'''

Includes a number of options for data set analysis.

'''

import sys,os
import matplotlib.pyplot as plt
import numpy as np

def runFormat(titlesize: float = 16, axeslabelsize: float = 16, legendfontsize: float = 16, xticksize: float = 16, yticksize: float = 16,
              ) -> None:
    
    plt.rcParams['figure.constrained_layout.use'] = True
    plt.rcParams["font.family"] = 'STIXGeneral'
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['axes.linewidth'] = 1
    plt.rc('axes', titlesize=titlesize)     
    plt.rc('axes', labelsize=axeslabelsize)    
    plt.rc('xtick', labelsize=xticksize)    
    plt.rc('ytick', labelsize=yticksize)    
    plt.rc('legend', fontsize=legendfontsize)
    plt.rc('figure', titlesize=16)


def MLIPTS_colour_palette(colourname):
    
    colourHex = {
        "Red": "#D91E1E",
        "PastelRed" : "#F24B4B", 
        "PastelOrange" : "#F28749", 
    
    }
    
   
    return colourHex[colourname]
    

def energy_distribution(dataset: str, plot: bool=False, savefig_str: str=None) -> list[float]:
    '''
    Reads the energies from a dataset file and plots the energy distribution.
    
    Returns energies in a list for personal plotting settings. 
    
    Currently supporting xtb.xyz format only
    '''
    
    energies = []
    
    data = open(dataset,'r').read()

    for line in data.splitlines():
        
        line_split = line.split()
        
        if 'energy=' in line:
            
            location = [i for i, s in enumerate(line_split) if 'energy=' in s][0]
            
            energy = float(line_split[location].split("=")[1])

            energies.append(energy)
            
    if plot==True:
        
        runFormat()
        fig,ax = plt.subplots()
        ax.hist(energies,bins='auto',color=MLIPTS_colour_palette('PastelRed'),edgecolor=MLIPTS_colour_palette('Red'))
        ax.set(title="Dataset Energy Distribution",xlabel='Energy (eV)',ylabel='Count')
        if savefig_str != None:
            plt.savefig(savefig_str)

        
        
    return energies


def force_distribution(dataset: str, return_type: str='max',plot: bool=False, savefig_str: str=None) -> list[float]:
    '''
    Reads forces from data sets, finds the maximum force and plots the distribution
    
    Returns 'max' or 'ave' forces in list for personal plotting style
    '''
        
    max_forces = []
    ave_forces = []
    
    data = open(dataset,'r').read()
    
    data_lines = data.splitlines()

    for i, line in enumerate(data_lines):
        
        line_split = line.split()
        
        if 'Properties=' in line:
            
            num_atoms = int(data_lines[i-1])
         
            forces_current_structure = []
            
            location = [i for i, s in enumerate(line_split) if 'Properties=' in s][0]
            
            properties = (line_split[location].split("=")[1])
            
            properties_split = properties.split(":")
            
            location = [i for i, s in enumerate(properties_split) if 'forces' in s][0]
            
            num_before_forces = sum(int(x) for x in properties_split[0:location] if x.isdigit())
            
            for j, atom_line in enumerate(data_lines[i+1:i+num_atoms+1]):
                
                atom_line_split = atom_line.split()
                
                force_str = atom_line_split[num_before_forces:num_before_forces+3]
                force = np.array([force_str[0],force_str[1],force_str[2]])
                
                forces_current_structure.append(np.linalg.norm(force))
                
            ave_force = sum(forces_current_structure) / num_atoms
            max_force = max(forces_current_structure)
            
            ave_forces.append(ave_force)
            max_forces.append(max_force)
    
    
    

    if return_type == 'ave':
        return_forces = ave_forces
        plot_label = 'Average'
    else:
        return_forces= max_forces
        plot_label = 'Max'
    
    if plot==True:
        
        runFormat()
        fig,ax = plt.subplots()
        ax.hist(return_forces,bins='auto',color=MLIPTS_colour_palette('PastelRed'),edgecolor=MLIPTS_colour_palette('Red'))
        ax.set(title=f"Dataset {plot_label} Force Distribution",xlabel=f'{plot_label} Force (eV/Å)',ylabel='Count')
        if savefig_str != None:
            plt.savefig(savefig_str)

        
        
        
if __name__ == '__main__':
    
    #energy_distribution('solvent_xtb.xyz',plot=True,savefig_str='test_energies.png')
    force_distribution('solvent_xtb.xyz',plot=True,savefig_str='test_forces.png')
    
    
    
    