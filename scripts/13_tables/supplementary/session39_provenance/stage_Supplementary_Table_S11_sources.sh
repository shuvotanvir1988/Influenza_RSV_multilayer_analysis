#!/usr/bin/env bash
set -euo pipefail

ROOT="$HOME/Influenza_RSV_Project"
SRC="$ROOT/results/session38_additional_analysis/frozen/Session38A"
OUT="$ROOT/results/session39_additional_supplementary_tables/staging"
DESKTOP="${DESKTOP:-$HOME/Desktop}"

mkdir -p "$OUT"

FILES=(
  "tables/Session38A_EVIDENCE_WEIGHTED_GENE_PRIORITY_MASTER_v1.0.tsv"
  "tables/Session38A_EVIDENCE_WEIGHTED_SUMMARY_v1.0.tsv"
  "tables/Session38A_RESPONSE_DEPENDENCY_BRIDGE_GENES_v1.0.tsv"
  "tables/Session38A_TOP30_EVIDENCE_WEIGHTED_GENES_v1.0.tsv"
  "tables/Session38A_WEIGHT_SENSITIVITY_CORRELATIONS_v1.0.tsv"
  "tables/Session38A_WEIGHT_SENSITIVITY_GENE_RANKS_v1.0.tsv"
)

for f in "${FILES[@]}"; do
  [[ -f "$SRC/$f" ]] || { echo "ERROR: missing $SRC/$f" >&2; exit 1; }
done

TAR="$OUT/Session39_SuppTable_S11_SOURCE_BUNDLE_v1.0.tar.gz"
MAN="$OUT/Session39_SuppTable_S11_SOURCE_SHA256_v1.0.txt"

(
  cd "$SRC"
  sha256sum "${FILES[@]}"
) > "$MAN"

tar -czf "$TAR" -C "$SRC" "${FILES[@]}"

cp -p "$TAR" "$MAN" "$DESKTOP/"

echo "=== SESSION 39 SUPPLEMENTARY TABLE S11 SOURCE BUNDLE READY ==="
echo "$TAR"
echo "$MAN"
echo
echo "Copied to Desktop:"
echo "$DESKTOP/$(basename "$TAR")"
echo "$DESKTOP/$(basename "$MAN")"
echo
echo "Bundle contents:"
tar -tzf "$TAR"
