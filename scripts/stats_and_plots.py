import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


analysis_type = snakemake.wildcards.type

df = pd.read_csv(snakemake.input["summary"], sep="\t")

# Convert methylation columns to numeric
df["Maternal_Methylation"] = pd.to_numeric(df["Maternal_Methylation"], errors="coerce")
df["Paternal_Methylation"] = pd.to_numeric(df["Paternal_Methylation"], errors="coerce")

df["Difference"] = df["Maternal_Methylation"] - df["Paternal_Methylation"]

# Remove rows where either is missing
paired = df.dropna(subset=["Maternal_Methylation", "Paternal_Methylation"]).copy()

x = paired["Maternal_Methylation"]
y = paired["Paternal_Methylation"]

# Compute difference
paired["Difference"] = x - y

# ---------- Correct Wilcoxon test ----------
if len(paired) == 0:
    wilcoxon_stat, wilcoxon_p = None, None
elif (paired["Difference"] == 0).all():
    wilcoxon_stat, wilcoxon_p = 0.0, 1.0   # all ties → uninformative
else:
    res = wilcoxon(x, y)
    wilcoxon_stat, wilcoxon_p = res.statistic, res.pvalue

# Save results table
paired.to_csv(snakemake.output["results_table"], sep="\t", index=False)
paired.to_excel(snakemake.output["results_xlsx"], index=False)

# ===============================
#            PLOTS
# ===============================

wilcoxon_label = f"Wilcoxon p = {wilcoxon_p:.3e}" if wilcoxon_p is not None else "Wilcoxon NA"

pdf_path = snakemake.output["pdf"]   # add to your rule outputs

with PdfPages(pdf_path) as pdf:

    # -------- Scatter plot --------
    plt.figure(figsize=(7,7))
    sns.scatterplot(x=y, y=x, edgecolor=None)
    plt.xlabel("Paternal Methylation")
    plt.ylabel("Maternal Methylation")
    plt.title(f"Maternal vs Paternal {analysis_type}\n{wilcoxon_label}")
    plt.savefig(snakemake.output["scatter"], dpi=300)
    plt.close()

    # -------- Boxplot --------
    plt.figure(figsize=(6,6))
    sns.boxplot(data=paired[["Maternal_Methylation", "Paternal_Methylation"]])
    plt.ylabel("Methylation Score")
    plt.title(f"Distribution of {analysis_type} Levels\n{wilcoxon_label}")
    plt.savefig(snakemake.output["boxplot"], dpi=300)
    plt.close()
    
    # ----- Histogram of differences -----
    plt.figure(figsize=(6,6))
    sns.histplot(df["Difference"].dropna(), kde=True, bins=50)
    plt.xlabel("Maternal - Paternal Methylation")
    plt.title(f"Distribution of Differences in {analysis_type}\n{wilcoxon_label}")
    plt.savefig(snakemake.output["histogram"], dpi=300)
    plt.close()
    
    # -------- Bland–Altman plot --------
    mean_vals = (x + y) / 2
    diff_vals = x - y

    plt.figure(figsize=(7,6))
    sns.scatterplot(x=mean_vals, y=diff_vals, alpha=0.6)
    plt.axhline(0, color="red", linestyle="--")
    plt.xlabel("Mean Methylation")
    plt.ylabel("Maternal − Paternal")
    plt.title(f"Bland–Altman Plot of {analysis_type}\n{wilcoxon_label}")
    pdf.savefig(); plt.close()

print("Wilcoxon statistic:", wilcoxon_stat)
print("Wilcoxon p-value:", wilcoxon_p)
print("Plots saved.")
