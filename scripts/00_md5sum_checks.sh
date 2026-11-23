#!/usr/bin/env bash
set -euo pipefail

# Run from the repository root, regardless of where this script is invoked.
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE:-$0}")" && pwd)
ROOT_DIR=$(cd -- "$SCRIPT_DIR/.." && pwd)
cd "$ROOT_DIR"

mkdir -p data/metadata
CHECKSUM_FILE=data/metadata/data-checksums.txt

# Find recursively and only include files, avoiding unmatched-shell-glob errors.
find data/source-data -type f ! -name '.DS_Store' -exec md5sum {} + > "$CHECKSUM_FILE"
find data/processed_data -type f -name '*.csv' -exec md5sum {} + >> "$CHECKSUM_FILE"
find data/xrefdata -type f -name '*.csv' -exec md5sum {} + >> "$CHECKSUM_FILE"
find results -type f -exec md5sum {} + >> "$CHECKSUM_FILE"

## Check all md5sums with:
# md5sum -c data/metadata/data-checksums.txt
