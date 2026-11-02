import pandas as pd
import sys

# Input and output files
input_file = sys.argv[1]
output_file = sys.argv[2]

# Read PAV matrix
df = pd.read_csv(input_file, sep="\t")

# Rename coordinate columns
df = df.rename(columns={
    "chrom": "#Chromosome",
    "start": "FeatureStart",
    "end": "FeatureStop"
})

# Remove the original name column
df = df.drop(columns=["name"], errors="ignore")

# Add required Panache annotation columns
df.insert(3, "Sequence_IUPAC_Plus", ".")
df.insert(4, "SimilarBlocks", ".")
df.insert(5, "Function", ".")

# Save formatted matrix
df.to_csv(output_file, sep="\t", index=False)

print(f"Formatted Panache matrix saved to {output_file}")
