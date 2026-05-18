#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

BASE = Path("/home/alina_tgrv/beegfs/IBS_SQ")

metaphlan_file = BASE / "results/metaphlan_metatranscriptome/joined/metaphlan_merged.tsv"
meta_file = BASE / "metadata/metadata_326_clean_v2.tsv"

species_matrix_out = BASE / "results/metaphlan_metatranscriptome/joined/species_sample_matrix_with_groups.tsv"
stats_out = BASE / "results/metaphlan_metatranscriptome/joined/species_stats_mwu.tsv"

# ---------- load ----------
mpa = pd.read_csv(metaphlan_file, sep="\t", comment="#")
meta = pd.read_csv(meta_file, sep="\t")

print("MetaPhlAn merged shape:", mpa.shape)
print("Metadata shape:", meta.shape)

feature_col = mpa.columns[0]
mpa[feature_col] = mpa[feature_col].astype(str)

# ---------- keep species only ----------
species = mpa[mpa[feature_col].str.contains(r"\|s__")].copy()
print("Species rows before collapse:", species.shape[0])

species["species"] = species[feature_col].str.split("|").str[-1]
sample_cols = [c for c in species.columns if c not in [feature_col, "species"]]

species[sample_cols] = species[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)

dup_n = species["species"].duplicated().sum()
print("Duplicated species labels:", dup_n)

species_collapsed = species.groupby("species", as_index=False)[sample_cols].sum()
print("Species rows after collapse:", species_collapsed.shape[0])

# ---------- transpose to samples x species ----------
species_mat = species_collapsed.set_index("species").T.reset_index()
species_mat = species_mat.rename(columns={"index": "sample"})
species_mat["sample"] = species_mat["sample"].astype(str)

print("Species matrix shape:", species_mat.shape)

# ---------- metadata ----------
print("Metadata columns:", meta.columns.tolist())

sample_col_meta = None
for c in ["sample", "Sample", "run", "Run", "SRR"]:
    if c in meta.columns:
        sample_col_meta = c
        break

if sample_col_meta is None:
    raise ValueError("Не найдена колонка sample в metadata.")

group_col_meta = None
for c in ["ibs_status", "IBS_status", "group", "_group", "status"]:
    if c in meta.columns:
        group_col_meta = c
        break

if group_col_meta is None:
    raise ValueError("Не найдена колонка группы в metadata.")

meta2 = meta[[sample_col_meta, group_col_meta]].copy()
meta2 = meta2.rename(columns={sample_col_meta: "sample", group_col_meta: "_group"})
meta2["sample"] = meta2["sample"].astype(str)
meta2["_group"] = meta2["_group"].astype(str).str.strip()

meta2["_group"] = meta2["_group"].replace({
    "HC": "Control",
    "control": "Control",
    "Control ": "Control",
    "ibs": "IBS",
    "IBS ": "IBS"
})

joined = meta2.merge(species_mat, on="sample", how="inner")
print("Joined shape:", joined.shape)
print("\nGroup counts:")
print(joined["_group"].value_counts(dropna=False))

joined.to_csv(species_matrix_out, sep="\t", index=False)
print("\nSaved species matrix:")
print(species_matrix_out)

# ---------- statistics ----------
species_cols = [c for c in joined.columns if c.startswith("s__")]

ibs_mask = joined["_group"].astype(str).str.lower().eq("ibs")
ctrl_mask = joined["_group"].astype(str).str.lower().eq("control")

print("\nIBS samples:", ibs_mask.sum())
print("Control samples:", ctrl_mask.sum())

if ibs_mask.sum() == 0 or ctrl_mask.sum() == 0:
    raise ValueError("После merge одна из групп пуста.")

results = []

for sp in species_cols:
    x = pd.to_numeric(joined.loc[ibs_mask, sp], errors="coerce").fillna(0.0).values
    y = pd.to_numeric(joined.loc[ctrl_mask, sp], errors="coerce").fillna(0.0).values

    prevalence_ibs = (x > 0).mean()
    prevalence_ctrl = (y > 0).mean()
    mean_ibs = x.mean()
    mean_ctrl = y.mean()
    median_ibs = np.median(x)
    median_ctrl = np.median(y)
    log2fc = np.log2((mean_ibs + 1e-9) / (mean_ctrl + 1e-9))

    try:
        stat = mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
        u = stat.statistic
        p = stat.pvalue
        effect_rbc = (2 * u) / (len(x) * len(y)) - 1
    except Exception:
        u = np.nan
        p = np.nan
        effect_rbc = np.nan

    results.append({
        "species": sp,
        "prevalence_IBS": prevalence_ibs,
        "prevalence_Control": prevalence_ctrl,
        "delta_prevalence": prevalence_ibs - prevalence_ctrl,
        "mean_IBS": mean_ibs,
        "mean_Control": mean_ctrl,
        "median_IBS": median_ibs,
        "median_Control": median_ctrl,
        "log2FC_IBS_vs_Control": log2fc,
        "U_statistic": u,
        "p_value": p,
        "effect_size_rbc": effect_rbc
    })

res = pd.DataFrame(results)

res["q_value"] = np.nan
ok = res["p_value"].notna()
if ok.sum() > 0:
    res.loc[ok, "q_value"] = multipletests(res.loc[ok, "p_value"], method="fdr_bh")[1]

res = res.sort_values(["q_value", "p_value", "delta_prevalence"], ascending=[True, True, False])
res.to_csv(stats_out, sep="\t", index=False)

print("\nSaved stats:")
print(stats_out)

print("\nTop 20 results:")
print(res.head(20).to_string(index=False))