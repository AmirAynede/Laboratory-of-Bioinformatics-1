import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Domain ID to check in InterProScan
DOMAIN_ID = "PF00014"  # Pfam ID for Kunitz domain

# ---------- 1. Load .class files ----------
def load_class_file(path, label):
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["Sequence_ID", "Label", "E-value", "Bit_score"])
    # Extract UniProt accession only (e.g., from sp|A0A6P8HC43|XYZ to A0A6P8HC43)
    df["Sequence_ID"] = df["Sequence_ID"].apply(lambda x: x.split("|")[1] if "|" in x else x)
    df["HMM_Prediction"] = df["Label"].map({1: "Kunitz", 0: "Non-Kunitz"})
    df["Set"] = label
    return df[["Sequence_ID", "HMM_Prediction", "Set"]]

set1 = load_class_file("set_1.class", "Set 1")
set2 = load_class_file("set_2.class", "Set 2")
hmm_df = pd.concat([set1, set2], ignore_index=True)

# ---------- 2. Load InterProScan TSVs ----------
def parse_interpro(path):
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    df = df[[0, 4]]  # Sequence ID and Pfam domain accession
    df.columns = ["Sequence_ID", "Domain"]
    df["InterPro_Annotation"] = df["Domain"].apply(lambda x: "Kunitz" if DOMAIN_ID in str(x) else "Other")
    grouped = df.groupby("Sequence_ID")["InterPro_Annotation"].apply(
        lambda x: "Kunitz" if "Kunitz" in x.values else "Non-Kunitz"
    ).reset_index()
    return grouped

interpro_pos = parse_interpro("interproscan_pos.tsv")
interpro_neg = parse_interpro("interproscan_neg.tsv")
interpro_all = pd.concat([interpro_pos, interpro_neg], ignore_index=True)

# ---------- 3. Merge HMM + InterPro ----------
merged = pd.merge(hmm_df, interpro_all, on="Sequence_ID", how="inner")
merged.to_csv("HMM_vs_InterPro_Merged.csv", index=False)

# ---------- 4. Plot Heatmap ----------
plt.figure(figsize=(8, 5))
conf_matrix = pd.crosstab(
    [merged["HMM_Prediction"], merged["Set"]],
    merged["InterPro_Annotation"]
)

sns.set(style="whitegrid")
sns.heatmap(conf_matrix, annot=True, fmt="d", cmap="Blues", cbar=False)
plt.title("HMM vs InterProScan: Kunitz Domain Detection")
plt.ylabel("HMM Prediction by Set")
plt.xlabel("InterPro Annotation")
plt.tight_layout()
plt.savefig("HMM_vs_InterPro_Heatmap.png", dpi=300)
plt.show()