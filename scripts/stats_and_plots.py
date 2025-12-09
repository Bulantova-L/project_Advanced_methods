import pandas as pd
import numpy as np
from scipy.stats import wilcoxon, ttest_rel
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv(snakemake.input["summary"], sep="\t")

# ================================
# Step 3: Methylation Quantification
# ================================
maternal = df["Maternal_Methylation"]
paternal = df["Paternal_Methylation"]

# Difference
df["Difference"] = df["Maternal_Methylation"] - df["Paternal_Methylation"]

# ================================
# Step 4: Statistical Analysis
# ================================

# Wilcoxon Signed-Rank Test
wilcoxon_pvals = []
t_test_pvals   = []

for m, p in zip(maternal, paternal):
    # Wilcoxon (paired nonparametric)
    try:
        w_p = wilcoxon([m], [p]).pvalue
    except:
        w_p = np.nan

    # Paired t-test
    try:
        t_p = ttest_rel([m], [p]).pvalue
    except:
        t_p = np.nan

    wilcoxon_pvals.append(w_p)
    t_test_pvals.append(t_p)

df["Wilcoxon_p"] = wilcoxon_pvals
df["Ttest_p"] = t_test_pvals

# Multiple testing correction
df["Wilcoxon_fdr"] = multipletests(df["Wilcoxon_p"], method="fdr_bh")[1]
df["Ttest_fdr"] = multipletests(df["Ttest_p"], method="fdr_bh")[1]

df.to_csv(snakemake.output["results_table"], sep="\t", index=False)

print("Saved statistics table:", snakemake.output["results_table"])

# ================================
# Step 5: Visualization
# ================================

# -------- Scatter plot --------
plt.figure(figsize=(7,7))
sns.scatterplot(
    x=df["Paternal_Methylation"],
    y=df["Maternal_Methylation"],
    edgecolor=None
)
plt.xlabel("Paternal methylation")
plt.ylabel("Maternal methylation")
plt.title("Maternal vs Paternal Methylation per Gene")
plt.savefig(snakemake.output["scatter"], dpi=300)
plt.close()

# -------- Volcano plot --------
plt.figure(figsize=(8,7))
sns.scatterplot(
    x=df["Difference"],
    y=-np.log10(df["Wilcoxon_p"]),
    edgecolor=None
)
plt.xlabel("Maternal - Paternal methylation")
plt.ylabel("-log10(Wilcoxon p-value)")
plt.title("Volcano Plot of Methylation Differences")
plt.savefig(snakemake.output["volcano"], dpi=300)
plt.close()

# -------- Boxplot --------
plt.figure(figsize=(6,6))
data = pd.DataFrame({
    "Maternal": maternal,
    "Paternal": paternal
})
sns.boxplot(data=data)
plt.title("Distribution of Methylation Levels")
plt.ylabel("Methylation score")
plt.savefig(snakemake.output["boxplot"], dpi=300)
plt.close()

print("Plots saved.")
