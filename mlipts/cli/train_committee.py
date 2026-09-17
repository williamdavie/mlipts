"""
William Davie.

Uses mlipts with a choice of [MACE] to generate and run a committee of mace models.
"""

import argparse
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor

from mlipts.active_learn import ActiveLearn


def run_train_mace(config_file: str, device: str = None):

    from mace.cli.run_train import run as run_train
    from mace.tools import build_default_arg_parser

    parser = build_default_arg_parser()
    cli_args = ["--config", config_file]
    if device is not None:  # overwrites
        cli_args += ["--device", device]

    args = parser.parse_args(cli_args)
    run_train(args)


def _train_worker(task: tuple):
    """
    Runs in its own process (spawned fresh, no inherited CUDA context).
    Pins this process to a single GPU before mace/torch ever get imported.
    """
    gpu_id, config_file = task
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    run_train_mace(config_file, device="cuda")


def run_train_committee(
    N_models: int,
    hpc_config: str,
    base_config: str,
    N_gpus: int,
    architecture: str = "mace",
    device: str = "cuda",
):

    active_learn = ActiveLearn(hpc_config=hpc_config)
    active_learn.define_commitee(base_config, N_models)
    # -> This sets up a directory named "model_configs" with:
    # |- config_#0
    # |- config_#1
    # ...

    if architecture != "mace":
        return None

    config_files = [f"model_configs/config_#{i}" for i in range(N_models)]

    if device == "cpu":
        # Smoke-test path: no parallelism, no CUDA_VISIBLE_DEVICES games,
        # just run every model sequentially on this process's CPU.
        for config_file in config_files:
            run_train_mace(config_file, device="cpu")
        return None

    # GPU path: distribute the committee across N_gpus workers.
    n_gpus = max(N_gpus, 1)
    gpu_ids = [i % n_gpus for i in range(N_models)]
    tasks = list(zip(gpu_ids, config_files))

    # spawn (not fork) so each worker gets a clean CUDA context
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=n_gpus, mp_context=ctx) as executor:
        list(executor.map(_train_worker, tasks))

    return None


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Train a committee of models.")

    parser.add_argument("--N_models", type=int, required=True)
    parser.add_argument("--base_config", type=str, required=True)
    parser.add_argument("--hpc_config", type=str, required=True)
    parser.add_argument(
        "--gpus",
        type=int,
        default=0,
        help="Number of GPUs to distribute training over. Ignored if --device cpu.",
    )
    parser.add_argument("--architecture", type=str, default="mace")
    parser.add_argument(
        "--device",
        type=str,
        choices=["cuda", "cpu"],
        default="cuda",
        help="Use 'cpu' for a serial smoke test with no parallelization.",
    )

    args = parser.parse_args()

    if args.device == "cuda" and args.gpus < 1:
        parser.error("--gpus must be >= 1 when --device cuda")

    run_train_committee(
        N_models=args.N_models,
        hpc_config=args.hpc_config,
        base_config=args.base_config,
        N_gpus=args.gpus,
        architecture=args.architecture,
        device=args.device,
    )
