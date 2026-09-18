#!/usr/bin/env bash
set -euo pipefail

# Downloads archives from a Zenodo record and unpacks them into a target directory.
#
#
# Usage:
#   zenodo_download_data.sh <zenodo_record_id> [data_dir]
#
# <zenodo_record_id> is the numeric id from the record's DOI
# (e.g. for 10.5281/zenodo.1234567 the record id is 1234567).
# See Data/README.md for both record ids once published.
#
# Requires: curl, jq, unzip, and md5sum (or md5 on macOS).

RECORD_ID="${1:?Usage: $0 <zenodo_record_id> [data_dir]}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${2:-$REPO_ROOT/Data}"
API_URL="https://zenodo.org/api/records/${RECORD_ID}"

mkdir -p "$DATA_DIR"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

md5_of() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | awk '{print $1}'
  else
    md5 -q "$1"
  fi
}

echo "Fetching record metadata from $API_URL"
curl -fsSL "$API_URL" -o "$TMP_DIR/record.json"

n_files=$(jq '.files | length' "$TMP_DIR/record.json")
if [ "$n_files" -eq 0 ]; then
  echo "Error: record $RECORD_ID has no files (or doesn't exist / isn't public)." >&2
  exit 1
fi
echo "Record has $n_files file(s)"

for i in $(seq 0 $((n_files - 1))); do
  fname=$(jq -r ".files[$i].key" "$TMP_DIR/record.json")
  url=$(jq -r ".files[$i].links.self" "$TMP_DIR/record.json")
  expected_md5=$(jq -r ".files[$i].checksum" "$TMP_DIR/record.json" | sed 's/^md5://')

  echo "Downloading $fname ..."
  curl -fsSL "$url" -o "$TMP_DIR/$fname"

  actual_md5="$(md5_of "$TMP_DIR/$fname")"
  if [ "$actual_md5" != "$expected_md5" ]; then
    echo "Checksum mismatch for $fname (expected $expected_md5, got $actual_md5)" >&2
    exit 1
  fi

  case "$fname" in
    *.zip)
      echo "Unpacking $fname -> $DATA_DIR/"
      unzip -qo "$TMP_DIR/$fname" -d "$DATA_DIR"
      ;;
    *)
      cp "$TMP_DIR/$fname" "$DATA_DIR/"
      ;;
  esac
done

echo "Done. Data downloaded into $DATA_DIR"
