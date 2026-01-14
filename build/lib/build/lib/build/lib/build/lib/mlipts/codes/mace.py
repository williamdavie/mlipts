'''
File containing mace specific functionality.

Copyright (c) 2022 ACEsuit/mace
'''

import sys
from mace.cli.eval_configs import main as mace_eval_configs_main

def eval_mace(configs, model, output): # exactly as defined in the MACE tutorials 
    sys.argv = ["program", "--configs", configs, "--model", model, "--output", output]
    mace_eval_configs_main()


def main():
    
    calc_type = sys.argv[1]
    
    if calc_type == 'eval_mace':
        configs = sys.argv[2] 
        model = sys.argv[3]
        output = sys.argv[4]
        eval_mace(configs,model,output)
    
    
if __name__ == "__main__":
    main()
    