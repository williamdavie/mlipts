# python -m mlipts.cli.find_best_model.py --best-model models

"""
Returns the name of the model with the lowest RMSE on validation set (final epoch).
"""

import argparse
import json
import math
from pathlib import Path

import yaml


def load_eval_entries(results_file: Path, head: str):
    """Yield parsed eval-mode dicts for the given validation head."""
    with results_file.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                continue
            if d.get("mode") != "eval":
                continue
            if d.get("head") != head:
                continue
            if d.get("epoch") is None:
                # Skip the pre-training baseline eval; not a real epoch.
                continue
            yield d


def final_epoch_entry(results_file: Path, head: str):
    latest = None
    for entry in load_eval_entries(results_file, head):
        if latest is None or entry["epoch"] > latest["epoch"]:
            latest = entry
    return latest


def is_valid_number(value) -> bool:
    if value is None:
        return False
    if isinstance(value, float) and math.isnan(value):
        return False
    return True


def find_best_model(
    results_dir: Path, head: str = "Default", metric: str = "rmse_e_per_atom"
):
    results_files = sorted(results_dir.glob("model_*_run-*.txt"))
    if not results_files:
        return None

    per_model_final = {}
    skipped = []

    for results_file in results_files:
        # Filename looks like model_2_run-994.txt -> model name is model_2
        model_name = results_file.stem.rsplit("_run-", 1)[0]
        entry = final_epoch_entry(results_file, head)
        if entry is None:
            skipped.append(model_name)
            continue

        rmse = entry.get(metric)
        if not is_valid_number(rmse):
            skipped.append(model_name)
            continue

        per_model_final[model_name] = entry

    if not per_model_final:
        return None

    best_model_name, best_entry = min(
        per_model_final.items(), key=lambda kv: kv[1][metric]
    )

    for name, entry in sorted(per_model_final.items()):
        marker = "  <-- best" if name == best_model_name else ""
        value = entry[metric]
        # rmse_e_per_atom / rmse_f are in eV / eV-per-Angstrom; report in meV for readability
        value_mev = value * 1e3

    print(best_model_name)
    return best_model_name


def resolve_results_dir(args, parser) -> Path:
    if args.base_config:
        with open(args.base_config, "r") as f:
            cfg = yaml.safe_load(f)
        results_dir = cfg.get("results_dir")
        if not results_dir:
            parser.error(f"No 'results_dir' key found in {args.base_config}")
        return Path(results_dir)
    return Path(args.results_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find the committee model with the lowest validation RMSE at its final epoch."
    )
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "--results_dir", type=str, help="MACE results_dir to search."
    )
    source_group.add_argument(
        "--base_config",
        type=str,
        help="Path to the MACE base_config.yml; its 'results_dir' field is used.",
    )
    parser.add_argument(
        "--head", type=str, default="Default", help="Validation head name"
    )
    parser.add_argument(
        "--metric",
        type=str,
        default="rmse_e_per_atom",
        choices=["rmse_e_per_atom", "rmse_f", "rmse_stress", "rmse_virials_per_atom"],
        help="RMSE metric to compare models on",
    )

    args = parser.parse_args()
    results_dir = resolve_results_dir(args, parser)
    find_best_model(results_dir, head=args.head, metric=args.metric)
