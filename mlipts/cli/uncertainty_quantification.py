"""
William Davie.

Uses mlipts with a choice of [MACE] to evaluate a set of samples and quantify uncertainty.
"""

import argparse


def run_uncertainty_quantification():

    return None


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Train a committee of models.")

    parser.add_argument(
        "--samples", type=str, required=True, help="Path to samples file."
    )
    parser.add_argument(
        "--N",
        type=int,
        required=True,
        help="Defines the number of output samples desired after UQ.",
    )
    parser.add_argument("--gpus", type=int, required=True)
    parser.add_argument("--output", type=str, required=True)

    parser.add_argument(
        "--models_dir", type=str, default="models", help="Directory of models"
    )
    parser.add_argument("--architechture", type=str, default="mace")

    args = parser.parse_args()
