#!/usr/bin/env python3

from pathlib import Path
import os
import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
)

INPUT_DIR = Path("results/random_forest/input_blocks")
OUTDIR = Path("results/random_forest/ibs_multiblock_python")
OUTDIR.mkdir(parents=True, exist_ok=True)

N_SPLITS = int(os.getenv("RF_N_SPLITS", "5"))
N_REPEATS = int(os.getenv("RF_N_REPEATS", "20"))
N_TREES = int(os.getenv("RF_NUM_TREES", "1000"))
N_JOBS = int(os.getenv("RF_N_JOBS", os.getenv("SLURM_CPUS_PER_TASK", "4")))

META_COLS = ["sample", "subject_id", "metabolomics_id", "ibs_status"]

input_files = sorted(INPUT_DIR.glob("rf_234_*.tsv"))

if not input_files:
    raise FileNotFoundError(f"No rf_234_*.tsv files found in {INPUT_DIR}")

summary_rows = []

print("Random Forest IBS classification")
print(f"Input dir: {INPUT_DIR}")
print(f"Output dir: {OUTDIR}")
print(f"CV: {N_REPEATS} repeats x {N_SPLITS} folds")
print(f"Trees: {N_TREES}")
print(f"Jobs: {N_JOBS}")
print()

for path in input_files:
    block_name = path.stem
    print("=" * 80)
    print(block_name)
    print(path)

    df = pd.read_csv(path, sep="\t")

    if "ibs_status" not in df.columns:
        raise ValueError(f"ibs_status column missing in {path}")

    y = (df["ibs_status"] == "IBS").astype(int).values

    feature_cols = [c for c in df.columns if c not in META_COLS]
    X = df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    # убираем константные признаки
    non_constant = X.var(axis=0) > 0
    X = X.loc[:, non_constant]
    feature_names = X.columns.to_numpy()

    print(f"Samples: {X.shape[0]}")
    print(f"Features after removing constant columns: {X.shape[1]}")
    print("Class counts:")
    print(pd.Series(df["ibs_status"]).value_counts())

    cv = RepeatedStratifiedKFold(
        n_splits=N_SPLITS,
        n_repeats=N_REPEATS,
        random_state=42,
    )

    fold_rows = []
    importance_list = []
    oof_rows = []

    for fold_id, (train_idx, test_idx) in enumerate(cv.split(X, y), start=1):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_train = y[train_idx]
        y_test = y[test_idx]

        model = RandomForestClassifier(
            n_estimators=N_TREES,
            max_features="sqrt",
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=42 + fold_id,
            n_jobs=N_JOBS,
        )

        model.fit(X_train, y_train)

        prob = model.predict_proba(X_test)[:, 1]
        pred = (prob >= 0.5).astype(int)

        auc = roc_auc_score(y_test, prob)
        pr_auc = average_precision_score(y_test, prob)
        acc = accuracy_score(y_test, pred)
        bal_acc = balanced_accuracy_score(y_test, pred)
        f1 = f1_score(y_test, pred)

        fold_rows.append({
            "block": block_name,
            "fold": fold_id,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "n_features": X.shape[1],
            "roc_auc": auc,
            "pr_auc": pr_auc,
            "accuracy": acc,
            "balanced_accuracy": bal_acc,
            "f1": f1,
        })

        importance_list.append(model.feature_importances_)

        tmp = df.iloc[test_idx][["sample", "subject_id", "metabolomics_id", "ibs_status"]].copy()
        tmp["block"] = block_name
        tmp["fold"] = fold_id
        tmp["true_label"] = y_test
        tmp["prob_IBS"] = prob
        tmp["pred_label"] = pred
        oof_rows.append(tmp)

    block_outdir = OUTDIR / block_name
    block_outdir.mkdir(parents=True, exist_ok=True)

    fold_metrics = pd.DataFrame(fold_rows)
    fold_metrics.to_csv(block_outdir / "cv_metrics.tsv", sep="\t", index=False)

    oof = pd.concat(oof_rows, axis=0)
    oof.to_csv(block_outdir / "out_of_fold_predictions.tsv", sep="\t", index=False)

    importances = np.vstack(importance_list)

    imp = pd.DataFrame({
        "feature": feature_names,
        "importance_mean": importances.mean(axis=0),
        "importance_sd": importances.std(axis=0),
    }).sort_values("importance_mean", ascending=False)

    imp.to_csv(block_outdir / "feature_importance.tsv", sep="\t", index=False)
    imp.head(50).to_csv(block_outdir / "top50_feature_importance.tsv", sep="\t", index=False)

    summary = {
        "block": block_name,
        "n_samples": X.shape[0],
        "n_features": X.shape[1],
        "roc_auc_mean": fold_metrics["roc_auc"].mean(),
        "roc_auc_sd": fold_metrics["roc_auc"].std(),
        "pr_auc_mean": fold_metrics["pr_auc"].mean(),
        "pr_auc_sd": fold_metrics["pr_auc"].std(),
        "balanced_accuracy_mean": fold_metrics["balanced_accuracy"].mean(),
        "balanced_accuracy_sd": fold_metrics["balanced_accuracy"].std(),
        "accuracy_mean": fold_metrics["accuracy"].mean(),
        "accuracy_sd": fold_metrics["accuracy"].std(),
        "f1_mean": fold_metrics["f1"].mean(),
        "f1_sd": fold_metrics["f1"].std(),
    }

    summary_rows.append(summary)

    print("ROC-AUC:", round(summary["roc_auc_mean"], 3), "+/-", round(summary["roc_auc_sd"], 3))
    print("PR-AUC:", round(summary["pr_auc_mean"], 3), "+/-", round(summary["pr_auc_sd"], 3))
    print("Balanced accuracy:", round(summary["balanced_accuracy_mean"], 3))
    print("Top features:")
    print(imp.head(10).to_string(index=False))
    print()

summary_df = pd.DataFrame(summary_rows).sort_values("roc_auc_mean", ascending=False)
summary_df.to_csv(OUTDIR / "rf_block_summary.tsv", sep="\t", index=False)

print("=" * 80)
print("FINAL BLOCK SUMMARY")
print(summary_df.to_string(index=False))
print()
print("Written:")
print(OUTDIR / "rf_block_summary.tsv")
print("Done.")
