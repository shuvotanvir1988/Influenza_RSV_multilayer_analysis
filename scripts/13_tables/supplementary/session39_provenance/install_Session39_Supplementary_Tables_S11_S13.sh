#!/usr/bin/env bash
set -euo pipefail

ROOT="$HOME/Influenza_RSV_Project"
DOWNLOADS="${DOWNLOADS:-$HOME/Downloads}"

BASE="$ROOT/results/session39_additional_supplementary_tables"
FROZEN="$BASE/frozen"
PROV="$BASE/provenance"
TMP="$BASE/.install_tmp"

mkdir -p "$FROZEN" "$PROV"
rm -rf "$TMP"
mkdir -p "$TMP"

declare -A ARCHIVES=(
  [S11]="Session39_Supplementary_Table_S11_v1.0_FROZEN.tar.gz"
  [S12]="Session39_Supplementary_Table_S12_v1.0_FROZEN.tar.gz"
  [S13]="Session39_Supplementary_Table_S13_v1.0_FROZEN.tar.gz"
)

declare -A TOPDIRS=(
  [S11]="Session39_Supplementary_Table_S11_v1.0_FROZEN"
  [S12]="Session39_Supplementary_Table_S12_v1.0_FROZEN"
  [S13]="Session39_Supplementary_Table_S13_v1.0_FROZEN"
)

echo "=== SESSION 39 SUPPLEMENTARY TABLE INSTALL ==="

for S in S11 S12 S13; do
  ARC="$DOWNLOADS/${ARCHIVES[$S]}"
  TOP="${TOPDIRS[$S]}"
  DEST="$FROZEN/$S"

  echo
  echo "[$S] Checking archive..."
  [[ -f "$ARC" ]] || {
    echo "ERROR: missing archive:"
    echo "  $ARC"
    exit 1
  }

  rm -rf "$TMP/$S"
  mkdir -p "$TMP/$S"
  tar -xzf "$ARC" -C "$TMP/$S"

  SRC="$TMP/$S/$TOP"
  [[ -d "$SRC" ]] || {
    echo "ERROR: expected extracted directory not found:"
    echo "  $SRC"
    exit 1
  }

  # Verify the internal freeze manifest before installation.
  FREEZE_MANIFEST=$(find "$SRC/provenance" -maxdepth 1 -type f \
    -name "*FREEZE_SHA256.txt" | head -n 1)

  [[ -n "${FREEZE_MANIFEST:-}" && -f "$FREEZE_MANIFEST" ]] || {
    echo "ERROR: freeze manifest not found for $S"
    exit 1
  }

  echo "[$S] Verifying frozen package before install..."
  (
    cd "$SRC"
    sha256sum -c "provenance/$(basename "$FREEZE_MANIFEST")"
  )

  echo "[$S] Installing..."
  rm -rf "$DEST"
  mkdir -p "$DEST"
  cp -a "$SRC/." "$DEST/"

  # Verify again after installation.
  INSTALLED_MANIFEST=$(find "$DEST/provenance" -maxdepth 1 -type f \
    -name "*FREEZE_SHA256.txt" | head -n 1)

  echo "[$S] Verifying installed copy..."
  (
    cd "$DEST"
    sha256sum -c "provenance/$(basename "$INSTALLED_MANIFEST")"
  )

  echo "[$S] INSTALLED + VERIFIED"
done

# ------------------------------------------------------------------
# Build project-level master manifest across all installed S11-S13 files
# ------------------------------------------------------------------
MASTER="$PROV/Session39_SUPPLEMENTARY_TABLES_S11_S13_MASTER_SHA256_v1.0.txt"

(
  cd "$ROOT"
  find \
    "results/session39_additional_supplementary_tables/frozen/S11" \
    "results/session39_additional_supplementary_tables/frozen/S12" \
    "results/session39_additional_supplementary_tables/frozen/S13" \
    -type f -print0 \
    | sort -z \
    | xargs -0 sha256sum
) > "$MASTER"

# Freeze inventory
INVENTORY="$PROV/Session39_SUPPLEMENTARY_TABLES_S11_S13_INVENTORY_v1.0.tsv"
{
  printf "supplementary_table\tcanonical_directory\tworkbook\tfreeze_manifest\n"
  for S in S11 S12 S13; do
    DEST="$FROZEN/$S"
    WB=$(find "$DEST/workbook" -maxdepth 1 -type f -name "*.xlsx" -printf "%f\n" | head -n 1)
    FM=$(find "$DEST/provenance" -maxdepth 1 -type f -name "*FREEZE_SHA256.txt" -printf "%f\n" | head -n 1)
    printf "%s\t%s\t%s\t%s\n" \
      "$S" \
      "results/session39_additional_supplementary_tables/frozen/$S" \
      "$WB" \
      "$FM"
  done
} > "$INVENTORY"

# Master-manifest checksum
sha256sum "$MASTER" "$INVENTORY" > \
  "$PROV/Session39_SUPPLEMENTARY_TABLES_S11_S13_CONTROL_SHA256_v1.0.txt"

rm -rf "$TMP"

echo
echo "=== INSTALLATION COMPLETE ==="
echo
echo "Canonical locations:"
echo "  $FROZEN/S11"
echo "  $FROZEN/S12"
echo "  $FROZEN/S13"
echo
echo "Master manifest:"
echo "  $MASTER"
echo
echo "Inventory:"
echo "  $INVENTORY"
echo
echo "Control hashes:"
cat "$PROV/Session39_SUPPLEMENTARY_TABLES_S11_S13_CONTROL_SHA256_v1.0.txt"
