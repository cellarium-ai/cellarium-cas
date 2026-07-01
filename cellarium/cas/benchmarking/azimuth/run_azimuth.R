#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Azimuth)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript run_azimuth.R <mtx_dir> <azimuth_reference> <output_csv_path>")
}

mtx_dir <- args[1]
azimuth_reference <- args[2]
output_csv_path <- args[3]

cat("Reading 10x matrix...\n")
counts <- Read10X(data.dir = mtx_dir)
query <- CreateSeuratObject(counts = counts)

cat("Running Azimuth annotation...\n")
query <- RunAzimuth(query = query, reference = azimuth_reference)

cat("Writing metadata to CSV...\n")
write.csv(query@meta.data, output_csv_path, row.names = TRUE)

cat("Done.\n")
