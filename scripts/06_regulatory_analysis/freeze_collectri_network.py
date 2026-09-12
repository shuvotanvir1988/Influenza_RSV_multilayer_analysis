from pathlib import Path
import hashlib
import json
import pandas as pd
import decoupler as dc
import importlib.metadata as md

ROOT = Path.home() / "Influenza_RSV_Project"

OUTDIR = ROOT / "results/regulatory_driver_analysis/design"
QCDIR = ROOT / "results/regulatory_driver_analysis/qc"
LOGDIR = ROOT / "logs/regulatory_driver_analysis/provenance"

OUTDIR.mkdir(parents=True, exist_ok=True)
QCDIR.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

network_file = OUTDIR / "COLLECTRI_HUMAN_FROZEN_v1.0.tsv.gz"
summary_file = QCDIR / "COLLECTRI_HUMAN_NETWORK_SUMMARY_v1.0.tsv"
schema_file = QCDIR / "COLLECTRI_HUMAN_SCHEMA_v1.0.txt"
metadata_file = OUTDIR / "COLLECTRI_HUMAN_FROZEN_METADATA_v1.0.json"
sha_file = LOGDIR / "COLLECTRI_HUMAN_FROZEN_v1.0.sha256"

print("Retrieving human CollecTRI...")
net = dc.op.collectri(
    organism="human",
    remove_complexes=False,
    license="academic",
    verbose=True
)

if not isinstance(net, pd.DataFrame):
    raise TypeError(f"Expected pandas DataFrame, received {type(net)}")

print("\n=== NETWORK SHAPE ===")
print(net.shape)

print("\n=== COLUMNS ===")
print(list(net.columns))

print("\n=== FIRST ROWS ===")
print(net.head(10).to_string(index=False))

required = {"source", "target", "weight"}
missing = required - set(net.columns)

if missing:
    raise ValueError(
        f"Required CollecTRI columns missing: {sorted(missing)}"
    )

# Standardize only row ordering; preserve all supplied columns.
sort_cols = ["source", "target"]
if "weight" in net.columns:
    sort_cols.append("weight")

net = net.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)

# Basic integrity checks
if net["source"].isna().any():
    raise ValueError("Missing TF/source identifiers found.")

if net["target"].isna().any():
    raise ValueError("Missing target identifiers found.")

if net["weight"].isna().any():
    raise ValueError("Missing regulatory weights found.")

n_edges = len(net)
n_sources = net["source"].nunique()
n_targets = net["target"].nunique()
n_positive = int((net["weight"] > 0).sum())
n_negative = int((net["weight"] < 0).sum())
n_zero = int((net["weight"] == 0).sum())
n_duplicate_pairs = int(
    net.duplicated(subset=["source", "target"]).sum()
)

summary = pd.DataFrame([{
    "resource": "CollecTRI",
    "organism": "human",
    "decoupler_version": md.version("decoupler"),
    "remove_complexes": False,
    "license": "academic",
    "edges": n_edges,
    "unique_sources": n_sources,
    "unique_targets": n_targets,
    "positive_weight_edges": n_positive,
    "negative_weight_edges": n_negative,
    "zero_weight_edges": n_zero,
    "duplicate_source_target_pairs": n_duplicate_pairs
}])

print("\n=== SUMMARY ===")
print(summary.to_string(index=False))

# Save exact frozen network
net.to_csv(
    network_file,
    sep="\t",
    index=False,
    compression="gzip"
)

summary.to_csv(
    summary_file,
    sep="\t",
    index=False
)

with open(schema_file, "w") as f:
    f.write("COLLECTRI HUMAN NETWORK SCHEMA v1.0\n\n")
    f.write(f"Shape: {net.shape}\n\n")
    f.write("Columns:\n")
    for c in net.columns:
        f.write(f"- {c}: {net[c].dtype}\n")

metadata = {
    "resource": "CollecTRI",
    "organism": "human",
    "retrieval_function": "decoupler.op.collectri",
    "decoupler_version": md.version("decoupler"),
    "remove_complexes": False,
    "license": "academic",
    "edges": n_edges,
    "unique_sources": n_sources,
    "unique_targets": n_targets,
    "positive_weight_edges": n_positive,
    "negative_weight_edges": n_negative,
    "zero_weight_edges": n_zero,
    "duplicate_source_target_pairs": n_duplicate_pairs,
    "frozen_network_path": str(network_file.relative_to(ROOT))
}

with open(metadata_file, "w") as f:
    json.dump(metadata, f, indent=2)

# SHA256 of frozen compressed file
h = hashlib.sha256()
with open(network_file, "rb") as f:
    for block in iter(lambda: f.read(1024 * 1024), b""):
        h.update(block)

digest = h.hexdigest()

with open(sha_file, "w") as f:
    f.write(f"{digest}  {network_file.relative_to(ROOT)}\n")

print("\n=== FROZEN OUTPUTS ===")
print(network_file.relative_to(ROOT))
print(summary_file.relative_to(ROOT))
print(schema_file.relative_to(ROOT))
print(metadata_file.relative_to(ROOT))

print("\nSHA256:")
print(digest)

print("\nCOLLECTRI NETWORK SUCCESSFULLY FROZEN.")
