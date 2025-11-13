import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Compute mean gene methylation per haplotype.")
parser.add_argument("--maternal", required=True, help="BEDTools intersection file for maternal haplotype")
parser.add_argument("--paternal", required=True, help="BEDTools intersection file for paternal haplotype")
parser.add_argument("--output", required=True, help="Output CSV file")
args = parser.parse_args()

# Load inputs
maternal = pd.read_csv(args.maternal, sep="\t", header=None)
paternal = pd.read_csv(args.paternal, sep="\t", header=None)

# Example: gene ID column = 3, methylation column = 17 (adjust if needed)
maternal_summary = maternal.groupby(3)[17].mean().reset_index()
paternal_summary = paternal.groupby(3)[17].mean().reset_index()

maternal_summary.columns = ["GeneID", "Maternal_Methylation"]
paternal_summary.columns = ["GeneID", "Paternal_Methylation"]

merged = pd.merge(maternal_summary, paternal_summary, on="GeneID", how="inner")
merged["Difference"] = merged["Maternal_Methylation"] - merged["Paternal_Methylation"]

merged.to_csv(args.output, index=False)
print("Saved:", args.output)