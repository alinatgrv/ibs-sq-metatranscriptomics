import csv
import os
import math
from collections import Counter

import matplotlib.pyplot as plt


QC_SUMMARY = "results/humann_metatranscriptome/joined_all/sample_qc_summary.tsv"
QC_FLAGGED = "results/humann_metatranscriptome/joined_all/sample_qc_flagged.tsv"
OUTDIR = "results/humann_metatranscriptome/qc_report_2026-04-13"
FIGDIR = os.path.join(OUTDIR, "figures")
TABDIR = os.path.join(OUTDIR, "tables")


def ensure_dirs():
    os.makedirs(FIGDIR, exist_ok=True)
    os.makedirs(TABDIR, exist_ok=True)


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def to_float(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


def to_int(x):
    try:
        return int(float(x))
    except Exception:
        return None


def quartiles(vals):
    vals = sorted(vals)
    n = len(vals)

    def median(arr):
        m = len(arr)
        if m == 0:
            return None
        if m % 2 == 1:
            return arr[m // 2]
        return (arr[m // 2 - 1] + arr[m // 2]) / 2

    med = median(vals)
    lower = vals[: n // 2]
    upper = vals[(n + 1) // 2 :]
    q1 = median(lower)
    q3 = median(upper)
    return q1, med, q3


def iqr_upper(vals):
    q1, med, q3 = quartiles(vals)
    iqr = q3 - q1
    upper = q3 + 1.5 * iqr
    return q1, med, q3, iqr, upper


def save_stats(rows):
    detected = [to_int(r["detected_pathways"]) for r in rows]
    gf_unmapped = [to_float(r["genefamilies_UNMAPPED"]) for r in rows]
    pa_unmapped = [to_float(r["pathabundance_UNMAPPED"]) for r in rows]
    pa_unintegrated = [to_float(r["pathabundance_UNINTEGRATED"]) for r in rows]

    with open(os.path.join(TABDIR, "qc_metric_stats.tsv"), "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["metric", "n", "min", "q1", "median", "q3", "iqr", "upper_outlier_threshold", "max"])

        for name, vals in [
            ("detected_pathways", detected),
            ("genefamilies_UNMAPPED", gf_unmapped),
            ("pathabundance_UNMAPPED", pa_unmapped),
            ("pathabundance_UNINTEGRATED", pa_unintegrated),
        ]:
            vals = [v for v in vals if v is not None and not math.isnan(v)]
            q1, med, q3, iqr, upper = iqr_upper(vals)
            writer.writerow([name, len(vals), min(vals), q1, med, q3, iqr, upper, max(vals)])


def write_flag_tables(rows):
    hard = []
    borderline = []
    other_flagged = []

    for r in rows:
        flags = r["qc_flags"]
        if flags == "ok":
            continue
        if "hard_low_pathways" in flags:
            hard.append(r)
        elif "borderline_low_pathways" in flags:
            borderline.append(r)
        else:
            other_flagged.append(r)

    for name, subset in [
        ("hard_low_pathways.tsv", hard),
        ("borderline_low_pathways.tsv", borderline),
        ("other_flagged_samples.tsv", other_flagged),
    ]:
        out = os.path.join(TABDIR, name)
        with open(out, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(subset)


def plot_hist(vals, title, xlabel, outfile, threshold=None):
    plt.figure(figsize=(8, 5))
    plt.hist(vals, bins=30, edgecolor="black")
    if threshold is not None:
        plt.axvline(threshold, linestyle="--", linewidth=2, label=f"IQR upper = {threshold:.2f}")
        plt.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Number of samples")
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()


def plot_rank(vals, labels, title, ylabel, outfile):
    pairs = sorted(zip(labels, vals), key=lambda x: x[1])
    x = list(range(1, len(pairs) + 1))
    y = [p[1] for p in pairs]

    plt.figure(figsize=(9, 5))
    plt.plot(x, y)
    plt.title(title)
    plt.xlabel("Samples ranked from low to high")
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()


def plot_scatter(rows, xkey, ykey, title, xlabel, ylabel, outfile):
    colors = {
        "ok": "#4daf4a",
        "borderline": "#ffb000",
        "hard": "#d62728",
        "other_flagged": "#7f7f7f",
    }

    groups = {
        "ok": [],
        "borderline": [],
        "hard": [],
        "other_flagged": [],
    }

    for r in rows:
        x = to_float(r[xkey])
        y = to_float(r[ykey])
        flags = r["qc_flags"]

        if flags == "ok":
            cat = "ok"
        elif "hard_low_pathways" in flags:
            cat = "hard"
        elif "borderline_low_pathways" in flags:
            cat = "borderline"
        else:
            cat = "other_flagged"

        groups[cat].append((x, y, r["sample"]))

    plt.figure(figsize=(8, 6))
    for cat in ["ok", "borderline", "hard", "other_flagged"]:
        pts = groups[cat]
        if not pts:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        plt.scatter(xs, ys, s=24, alpha=0.75, label=cat, color=colors[cat])

    # annotate hard samples only
    for x, y, s in groups["hard"]:
        plt.annotate(s, (x, y), fontsize=7, alpha=0.85)

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outfile, dpi=220)
    plt.close()


def plot_flag_counts(rows):
    cnt = Counter(r["qc_flags"] for r in rows)
    items = cnt.most_common()

    labels = [k for k, _ in items]
    values = [v for _, v in items]

    plt.figure(figsize=(10, 6))
    plt.bar(range(len(labels)), values)
    plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
    plt.ylabel("Number of samples")
    plt.title("Counts of QC flag categories")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "09_qc_flag_counts.png"), dpi=220)
    plt.close()


def plot_pipeline_overview():
    stages = [
        ("QC sample dirs", 1184),
        ("HUMAnN sample dirs", 326),
        ("HUMAnN success", 326),
        ("MetaPhlAn profiles", 326),
        ("MetaPhlAn success", 325),
    ]

    labels = [x[0] for x in stages]
    values = [x[1] for x in stages]

    plt.figure(figsize=(8, 5))
    plt.bar(range(len(labels)), values)
    plt.xticks(range(len(labels)), labels, rotation=20, ha="right")
    plt.ylabel("Count")
    plt.title("Pipeline overview")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "01_pipeline_overview.png"), dpi=220)
    plt.close()


def write_overview_text(rows):
    cnt = Counter(r["qc_flags"] for r in rows)
    with open(os.path.join(TABDIR, "qc_overview.txt"), "w") as f:
        f.write("QC overview for HUMAnN metatranscriptome outputs\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total samples in sample_qc_flagged.tsv: {len(rows)}\n\n")
        f.write("QC flag counts:\n")
        for k, v in cnt.most_common():
            f.write(f"  {v}\t{k}\n")

        f.write("\nInterpretation notes:\n")
        f.write("- UNMAPPED (genefamilies): signal not assigned to known gene families by HUMAnN.\n")
        f.write("- UNMAPPED (pathabundance): abundance not assigned to reconstructed pathways.\n")
        f.write("- UNINTEGRATED: annotated functions present, but not integrated into pathway reconstruction.\n")
        f.write("- hard_low_pathways: samples with extremely low number of detected pathways.\n")
        f.write("- borderline_low_pathways: weak but not catastrophic pathway detection.\n")
        f.write("- high_* flags: high values by IQR-based outlier threshold.\n")


def main():
    ensure_dirs()

    rows = read_tsv(QC_FLAGGED)
    save_stats(rows)
    write_flag_tables(rows)
    write_overview_text(rows)

    samples = [r["sample"] for r in rows]
    detected = [to_int(r["detected_pathways"]) for r in rows]
    gf_unmapped = [to_float(r["genefamilies_UNMAPPED"]) for r in rows]
    pa_unmapped = [to_float(r["pathabundance_UNMAPPED"]) for r in rows]
    pa_unintegrated = [to_float(r["pathabundance_UNINTEGRATED"]) for r in rows]

    _, _, _, _, detected_thr = iqr_upper([v for v in detected if v is not None])
    _, _, _, _, gf_thr = iqr_upper([v for v in gf_unmapped if not math.isnan(v)])
    _, _, _, _, pau_thr = iqr_upper([v for v in pa_unmapped if not math.isnan(v)])
    _, _, _, _, pai_thr = iqr_upper([v for v in pa_unintegrated if not math.isnan(v)])

    plot_pipeline_overview()

    plot_hist(
        detected,
        "Detected pathways per sample",
        "Number of detected pathways",
        os.path.join(FIGDIR, "02_detected_pathways_hist.png"),
        threshold=detected_thr,
    )

    plot_hist(
        gf_unmapped,
        "Genefamilies UNMAPPED",
        "UNMAPPED abundance",
        os.path.join(FIGDIR, "03_gf_unmapped_hist.png"),
        threshold=gf_thr,
    )

    plot_hist(
        pa_unmapped,
        "Pathabundance UNMAPPED",
        "UNMAPPED abundance",
        os.path.join(FIGDIR, "04_pa_unmapped_hist.png"),
        threshold=pau_thr,
    )

    plot_hist(
        pa_unintegrated,
        "Pathabundance UNINTEGRATED",
        "UNINTEGRATED abundance",
        os.path.join(FIGDIR, "05_pa_unintegrated_hist.png"),
        threshold=pai_thr,
    )

    plot_scatter(
        rows,
        "detected_pathways",
        "genefamilies_UNMAPPED",
        "Detected pathways vs genefamilies UNMAPPED",
        "Detected pathways",
        "Genefamilies UNMAPPED",
        os.path.join(FIGDIR, "06_scatter_detected_vs_gf_unmapped.png"),
    )

    plot_scatter(
        rows,
        "detected_pathways",
        "pathabundance_UNMAPPED",
        "Detected pathways vs pathabundance UNMAPPED",
        "Detected pathways",
        "Pathabundance UNMAPPED",
        os.path.join(FIGDIR, "07_scatter_detected_vs_pa_unmapped.png"),
    )

    plot_scatter(
        rows,
        "detected_pathways",
        "pathabundance_UNINTEGRATED",
        "Detected pathways vs pathabundance UNINTEGRATED",
        "Detected pathways",
        "Pathabundance UNINTEGRATED",
        os.path.join(FIGDIR, "08_scatter_detected_vs_pa_unintegrated.png"),
    )

    plot_flag_counts(rows)

    plot_rank(
        detected,
        samples,
        "Ranked samples by detected pathways",
        "Detected pathways",
        os.path.join(FIGDIR, "10_rank_detected_pathways.png"),
    )

    plot_rank(
        gf_unmapped,
        samples,
        "Ranked samples by genefamilies UNMAPPED",
        "Genefamilies UNMAPPED",
        os.path.join(FIGDIR, "11_rank_gf_unmapped.png"),
    )

    plot_rank(
        pa_unintegrated,
        samples,
        "Ranked samples by pathabundance UNINTEGRATED",
        "Pathabundance UNINTEGRATED",
        os.path.join(FIGDIR, "12_rank_pa_unintegrated.png"),
    )

    print("QC figures and tables written to:", OUTDIR)


if __name__ == "__main__":
    main()