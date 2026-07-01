#!/bin/bash
set -e

if [[ "$#" -lt 3 ]]; then
  echo "Usage: run.sh <input_h5ad_path> <azimuth_reference> <output_csv_path>"
  exit 1
fi

INPUT_H5AD=$1
AZIMUTH_REF=$2
OUTPUT_CSV=$3
MTX_DIR="${INPUT_H5AD}_10x"

echo "Converting h5ad to 10x format..."
/opt/py311/bin/python /scripts/convert_h5ad.py "$INPUT_H5AD" "$MTX_DIR"

echo "Running Azimuth annotation..."
Rscript /scripts/run_azimuth.R "$MTX_DIR" "$AZIMUTH_REF" "$OUTPUT_CSV"

echo "Cleaning up intermediate 10x files..."
rm -rf "$MTX_DIR"