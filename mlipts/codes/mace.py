"""
File containing mace specific functionality.

Copyright (c) 2022 ACEsuit/mace
Batatia, I., Kovacs, D.P., Simm, G., Ortner, C. and Csányi, G., 2022. MACE: Higher order equivariant message passing neural networks for fast and accurate force fields. Advances in neural information processing systems, 35, pp.11423-11436.
"""

import sys

from mace.cli.eval_configs import main as mace_eval_configs_main


def eval_mace(configs, model, output):
    """
    exactly as defined in the MACE tutorials.
    https://colab.research.google.com/drive/1ZrTuTvavXiCxTFyjBV4GqlARxgFwYAtX
    """

    sys.argv = ["program", "--configs", configs, "--model", model, "--output", output]
    mace_eval_configs_main()


def main():

    calc_type = sys.argv[1]

    if calc_type == "eval_mace":
        configs = sys.argv[2]
        model = sys.argv[3]
        output = sys.argv[4]
        eval_mace(configs, model, output)


if __name__ == "__main__":
    main()
