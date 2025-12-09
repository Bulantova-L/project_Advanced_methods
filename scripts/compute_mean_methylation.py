import pandas as pd

maternal_file = snakemake.input["maternal_overlap"]
paternal_file = snakemake.input["paternal_overlap"]
output_file   = snakemake.output[0]

# Load inputs
maternal = pd.read_csv(maternal_file, sep="\t", header=None)
paternal = pd.read_csv(paternal_file, sep="\t", header=None)

# Example: gene ID column = 3, methylation column = 17 (adjust if needed)
maternal_summary = maternal.groupby(3)[17].mean().reset_index()
paternal_summary = paternal.groupby(3)[17].mean().reset_index()

maternal_summary.columns = ["GeneID", "Maternal_Methylation"]
paternal_summary.columns = ["GeneID", "Paternal_Methylation"]

merged = pd.merge(maternal_summary, paternal_summary, on="GeneID", how="inner")
merged["Difference"] = merged["Maternal_Methylation"] - merged["Paternal_Methylation"]

merged.to_csv(output_file, index=False)
print("Saved:", output_file)