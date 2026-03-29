# %%
from pathlib import Path

# Base path
BASE = Path("/home/ggiordano/snap/main")

# Project structure
DIRS = {
    "scripts": BASE / "scripts/trpv project",
    "export": BASE / "export/trpv_project",
    "data": BASE / "data/trpv_project",
    "import": BASE / "import/trpv_project"
}

# Create dirs if they don't exist (safe)
for name, path in DIRS.items():
    path.mkdir(parents=True, exist_ok=True)

# Optional: print for sanity check
for k, v in DIRS.items():
    print(f"{k}: {v}")

# %%
import pandas as pd
import numpy as np
from pathlib import Path

IMPORT_DIR = DIRS["import"]

def load_expression_tsv(file_path):
    # Read file, skipping comment lines
    with open(file_path) as f:
        lines = [l for l in f if not l.startswith("#")]

    # Build DataFrame
    df = pd.read_csv(
        pd.io.common.StringIO("".join(lines)),
        sep="\t"
    )

    # First column = experiment names
    df = df.rename(columns={df.columns[0]: "experiment"})
    df = df.set_index("experiment")

    # Convert to numeric
    df = df.replace("NA", np.nan)
    df = df.apply(pd.to_numeric, errors="coerce")

    return df
# Output directory
SAVE_DIR = DIRS["data"] / "expression_matrices_tsv"
SAVE_DIR.mkdir(parents=True, exist_ok=True)

# Load all TSV files into expression_data
expression_data = {}

for file in IMPORT_DIR.glob("*.tsv"):
    gene_name = file.stem.lower().replace(" ", "_")
    expression_data[gene_name] = load_expression_tsv(file)

print(f"Loaded {len(expression_data)} genes")
print(list(expression_data.keys()))

# Save each gene matrix as TSV
for gene, df in expression_data.items():
    out_file = SAVE_DIR / f"{gene}.tsv"
    
    df.to_csv(
        out_file,
        sep="\t",
        na_rep="NA"  # keep consistency with original files
    )

print(f"Saved {len(expression_data)} TSV files to:")
print(SAVE_DIR)

# %%
#Get col names 
from pathlib import Path

LOAD_DIR = DIRS["data"] / "expression_matrices_tsv"

columns_per_gene = {}
all_columns = set()

# Collect columns
for file in LOAD_DIR.glob("*.tsv"):
    gene = file.stem
    df = pd.read_csv(file, sep="\t", nrows=0)  # only header
    
    cols = set(df.columns) - {"experiment"}  # exclude index column
    columns_per_gene[gene] = cols
    all_columns.update(cols)

# Compare
differences = {}

for gene, cols in columns_per_gene.items():
    missing = all_columns - cols
    extra = cols - all_columns  # usually empty
    
    if missing:
        differences[gene] = {
            "missing_cols": sorted(missing)
        }

# Output
print(f"Total unique tissues: {len(all_columns)}")
print(f"Genes with differences: {len(differences)}\n")

for gene, diff in differences.items():
    print(f"{gene}: missing {len(diff['missing_cols'])} columns")

# %%
# Extract GTEx-only matrices
gtex_data = {}

for gene, df in expression_data.items():
    # Select GTEx row
    gtex_df = df[df.index.str.contains("GTEX", case=False, na=False)]
    
    if len(gtex_df) == 0:
        continue
    
    # Keep only first GTEx row (there is usually one)
    gtex_row = gtex_df.iloc[[0]]
    
    gtex_data[gene] = gtex_row

print(f"Extracted GTEx for {len(gtex_data)} genes")
# Create output directory
SAVE_DIR = DIRS["data"] / "gtex_only_raw"
SAVE_DIR.mkdir(parents=True, exist_ok=True)

# Save each gene GTEx row
for gene, df in gtex_data.items():
    out_file = SAVE_DIR / f"{gene}_gtex.tsv"
    
    df.to_csv(
        out_file,
        sep="\t",
        na_rep="NA"
    )

print(f"Saved {len(gtex_data)} GTEx-only files to:")
print(SAVE_DIR)

from pathlib import Path

# Output directory
SAVE_DIR_CLEAN = DIRS["data"] / "gtex_only_noNA"
SAVE_DIR_CLEAN.mkdir(parents=True, exist_ok=True)

clean_gtex_data = {}
columns_per_gene = {}

# Function to remove NA columns
def drop_na_columns(df):
    return df.dropna(axis=1)

# Process each gene
for gene, df in gtex_data.items():
    df_clean = drop_na_columns(df)
    
    clean_gtex_data[gene] = df_clean
    columns_per_gene[gene] = set(df_clean.columns)
    
    # Save cleaned file
    out_file = SAVE_DIR_CLEAN / f"{gene}_gtex_noNA.tsv"
    df_clean.to_csv(out_file, sep="\t", na_rep="NA")

print(f"Saved cleaned GTEx (no NA columns) for {len(clean_gtex_data)} genes")
print(SAVE_DIR_CLEAN)
# Get all column sets
col_sets = list(columns_per_gene.values())

# Check equality
all_equal = all(col_sets[0] == s for s in col_sets)

print("All genes have identical columns:", all_equal)
if not all_equal:
    all_cols = set.union(*columns_per_gene.values())
    
    for gene, cols in columns_per_gene.items():
        missing = all_cols - cols
        if missing:
            print(f"{gene}: missing {len(missing)} columns")

# Combine into one matrix (genes × tissues)
combined_df = pd.DataFrame({
    gene: df.iloc[0] for gene, df in clean_gtex_data.items()
}).T

# Optional: sort genes (trpv1 → trpv6)
combined_df = combined_df.sort_index()

# Output directory
SAVE_DIR_COMBINED = DIRS["data"] / "gtex_only_combined"
SAVE_DIR_COMBINED.mkdir(parents=True, exist_ok=True)

# Save file
out_file = SAVE_DIR_COMBINED / "trpv_gtex_combined.tsv"

combined_df.to_csv(
    out_file,
    sep="\t",
    na_rep="NA"
)

print("Saved combined GTEx matrix:")
print(out_file)
# Columns to remove (exact names)
REMOVE_COLS = [
    "adrenal gland",
    "cervical vertebra",
    "coronary artery",
    "ectocervix",
    "endocervix",
    "esophagogastric junction",
    "esophagus mucosa",
    "fallopian tube",
    "mammary gland",
    "muscle layer of esophagus",
    "omental fat pad",
    "ovary",
    "pituitary gland",
    "skeletal muscle tissue",
    "subcutaneous adipose tissue",
    "thyroid gland",
    "tibial artery",
    "tibial nerve",
    "uterus",
    "vagina", 
    "aorta"
]

# Filter columns (safe: ignore missing just in case)
combined_filtered = combined_df.drop(
    columns=[c for c in REMOVE_COLS if c in combined_df.columns]
)

print(f"Original shape: {combined_df.shape}")
print(f"Filtered shape: {combined_filtered.shape}")
SAVE_DIR_FILTERED = DIRS["data"] / "gtex_only_combined_filtered"
SAVE_DIR_FILTERED.mkdir(parents=True, exist_ok=True)

out_file = SAVE_DIR_FILTERED / "trpv_gtex_combined_filtered.tsv"

combined_filtered.to_csv(
    out_file,
    sep="\t",
    na_rep="NA"
)

print("Saved filtered combined matrix:")
print(out_file)

# Define pooling groups
POOLING = {
    "Brain_pooled": [
        "amygdala",
        "anterior cingulate cortex",
        "caudate nucleus",
        "cerebellar hemisphere",
        "cerebellum",
        "cerebral cortex",
        "frontal cortex",
        "hippocampus",
        "hypothalamus",
        "nucleus accumbens",
        "putamen",
        "substantia nigra"
    ],
    "Heart_pooled": [
        "heart left ventricle",
        "atrium auricular region"
    ],
    "Kidney_pooled": [
        "cortex of kidney",
        "renal medulla"
    ],
    "Intestine_pooled": [
        "ileum",
        "sigmoid colon",
        "transverse colon"
    ],
    "Skin_pooled": [
        "skin",
        "lower leg skin",
        "suprapubic skin"
    ]
}
# Copy original matrix
combined_pooled = combined_filtered.copy()

# Keep track of columns to remove
cols_to_remove = set()

# Create pooled columns
for new_col, cols in POOLING.items():
    valid_cols = [c for c in cols if c in combined_pooled.columns]
    
    if len(valid_cols) == 0:
        print(f"Warning: no valid columns found for {new_col}")
        continue
    
    # Compute mean
    combined_pooled[new_col] = combined_pooled[valid_cols].mean(axis=1)
    
    # Mark original columns for removal
    cols_to_remove.update(valid_cols)

# Remove pooled source columns
combined_pooled = combined_pooled.drop(columns=list(cols_to_remove))

print("Added pooled columns:")
print([c for c in combined_pooled.columns if "pooled" in c])

print(f"Final shape after pooling: {combined_pooled.shape}")
SAVE_DIR_POOLED = DIRS["data"] / "gtex_only_combined_pooled"
SAVE_DIR_POOLED.mkdir(parents=True, exist_ok=True)

out_file = SAVE_DIR_POOLED / "trpv_gtex_combined_pooled.tsv"

combined_pooled.to_csv(
    out_file,
    sep="\t",
    na_rep="NA"
)

print("Saved pooled matrix:")
print(out_file)



# %%
# =========================================================
# TRPV heatmap (your data, publication style)
# =========================================================

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.font_manager as fm
from pathlib import Path

# ---------------------------------------------------------
# Font (Arial)
# ---------------------------------------------------------

font_dir = Path.home() / ".local/share/fonts"

for font in font_dir.glob("Arial*.TTF"):
    fm.fontManager.addfont(str(font))

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.dpi": 300
})

# ---------------------------------------------------------
# Figure size (9 cm width panel)
# ---------------------------------------------------------

cm_to_inch = 1 / 2.54
fig_width = 9 * cm_to_inch
fig_height = 4 * cm_to_inch

# ---------------------------------------------------------
# Prepare data (FIXED)
# ---------------------------------------------------------

df = combined_pooled.copy()

# Standardize column names (important!)
df.columns = df.columns.str.lower().str.strip()

# Rename pooled + key columns to clean labels
df = df.rename(columns={
    "brain_pooled": "Brain",
    "heart_pooled": "Heart",
    "kidney_pooled": "Kidney",
    "intestine_pooled": "Intestine",
    "skin_pooled": "Skin",
    "blood": "Blood",
    "liver": "Liver",
    "lung": "Lung",
    "pancreas": "Pancreas",
    "spleen": "Spleen",
    "testis": "Testis",
    "prostate gland": "Prostate",
    "minor salivary gland": "Salivary gland",
    "urinary bladder": "Urinary bladder",
    "stomach": "Stomach"   
})

# Uppercase gene names
df.index = [g.upper() for g in df.index]

# Desired order
column_order = [
    "Blood",
    "Brain",
    "Heart",
    "Intestine",
    "Kidney",
    "Liver",
    "Lung",
    "Pancreas",
    "Prostate",
    "Salivary gland",
    "Skin",
    "Spleen",
    "Stomach",
    "Testis",
    "Urinary bladder"
]

# Keep only available columns (NOW this will work correctly)
df = df[[c for c in column_order if c in df.columns]]

# ---------------------------------------------------------
# Log transform
# ---------------------------------------------------------

df = df.apply(pd.to_numeric, errors="coerce")
df_log = np.log2(df + 1)

# Convert 0 values to NaN so they appear as grey
df_log = df_log.mask(df_log == 0)
# ---------------------------------------------------------
# Colormap with NA = light grey
# ---------------------------------------------------------

cmap = plt.cm.Blues.copy()
cmap.set_bad(color="#d9d9d9")

# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

sns.heatmap(
    df_log,
    cmap=cmap,
    linewidths=0.15,
    linecolor="#eeeeee",
    square=True,
    cbar=True,
    cbar_kws={
        "shrink": 0.5,
        "aspect": 15,
        "pad": 0.02
    },
    ax=ax
)

# ---------------------------------------------------------
# Labels
# ---------------------------------------------------------

ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

ax.tick_params(axis="y", which="both", length=2, width=0.3, direction="out")
ax.tick_params(axis="x", which="both", length=2, width=0.3, direction="out")

# ---------------------------------------------------------
# Light frame (panel style)
# ---------------------------------------------------------

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.3)
    spine.set_color("#888888")

# ---------------------------------------------------------
# Colorbar styling
# ---------------------------------------------------------

cbar = ax.collections[0].colorbar
cbar.set_label("log2(TPM + 1)", fontsize=6)
cbar.ax.tick_params(labelsize=6)
cbar.outline.set_linewidth(0.3)
cbar.outline.set_edgecolor("#888888")

# ---------------------------------------------------------
# Layout
# ---------------------------------------------------------

plt.tight_layout(pad=0.2)

# ---------------------------------------------------------
# Save
# ---------------------------------------------------------

out_dir = DIRS["export"]
out_dir.mkdir(parents=True, exist_ok=True)

plt.savefig(out_dir / "TRPV_heatmap_from_GTEX.pdf", bbox_inches="tight")
plt.savefig(out_dir / "TRPV_heatmap_from_GTEX.jpeg", dpi=300, bbox_inches="tight")
plt.savefig(out_dir / "TRPV_heatmap_from_GTEX.svg", dpi=300, bbox_inches="tight")

plt.show()

# %%



