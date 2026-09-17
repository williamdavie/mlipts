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


def parse_md_logs(models_dir="models"):
    """
    Parse model_*_debug.log files in a directory and extract
    stage 1 and stage 2 error metrics (RMSE_E, RMSE_F) for
    train and valid sets.

    Returns
    -------
    stage1, stage2 : dict
        {model_id: {"train": {"RMSE_E": float, "RMSE_F": float},
                    "valid": {"RMSE_E": float, "RMSE_F": float}}}
    """
    stage1 = {}
    stage2 = {}

    log_files = sorted(glob.glob(os.path.join(models_dir, "model_*_debug.log")))

    # matches a row like: | train_Default |   18.5   |   64.8   | ...
    row_re = re.compile(
        r"\|\s*(train|valid)_Default\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
    )
    model_id_re = re.compile(r"model_(\d+)_")

    for log_path in log_files:
        m = model_id_re.search(os.path.basename(log_path))
        if not m:
            continue
        model_id = int(m.group(1))

        with open(log_path, "r") as f:
            content = f.read()

        table_blocks = content.split("Error-table on TRAIN and VALID:")[1:]

        stage_dicts = [stage1, stage2]

        for stage_idx, block in enumerate(table_blocks[:2]):
            metrics = {}
            for match in row_re.finditer(block):
                split_name, rmse_e, rmse_f = match.groups()
                metrics[split_name] = {
                    "RMSE_E": float(rmse_e),
                    "RMSE_F": float(rmse_f),
                }
            if metrics:
                stage_dicts[stage_idx][model_id] = metrics


def return_best_model(
    models_dir="models",
    stage=2,
    metric="both",
    split="valid",
    weight_e=0.5,
    weight_f=0.5,
):
    """
    Determine the best performing model based on parsed MD log metrics.

    To be used for active learning MD runs.

    Returns
    -------
    str
        e.g. "model_2" for the best performing model.
    """
    stage1, stage2 = parse_md_logs(models_dir)

    if stage == 1:
        stage_dict = stage1
    elif stage == 2:
        stage_dict = stage2
    else:
        raise ValueError("stage must be 1 or 2")

    if not stage_dict:
        return None

    best_model_id = None
    best_score = float("inf")

    for model_id, splits in stage_dict.items():
        if split not in splits:
            continue

        metrics = splits[split]

        if metric == "energy":
            score = metrics["RMSE_E"]
        elif metric == "force":
            score = metrics["RMSE_F"]
        elif metric == "both":
            score = weight_e * metrics["RMSE_E"] + weight_f * metrics["RMSE_F"]
        else:
            raise ValueError("metric must be 'energy', 'force', or 'both'")

        if score < best_score:
            best_score = score
            best_model_id = model_id

    if best_model_id is None:
        return None

    return f"model_{best_model_id}"


if __name__ == "__main__":
    best = return_best_model(models_dir="models", stage=2, metric="both", split="valid")
    print("Best model:", best)
