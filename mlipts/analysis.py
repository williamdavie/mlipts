'''

Includes a number of options for data set analysis.

'''

import sys,os
import matplotlib.pyplot as plt

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
    

    
    


def energy_distribution(dataset: str, plot: bool=False, savefig_str: str=None):
    '''
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
        
        
if __name__ == '__main__':
    
    energy_distribution('solvent_xtb.xyz',plot=True,savefig_str='test_energies.png')
    
    
    
    