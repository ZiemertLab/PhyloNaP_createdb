#!/usr/bin/env Rscript
# Step 13 helper: MAD rooting for a single tree.
#
# Usage:  Rscript 13_root_trees_mad.R <input.nwk> <output_rooted.nwk>
#
# Requires: ape, phangorn, inline, Rcpp
# The madRoot.R source and Phylib headers must be available.
# Set MAD_ROOT_FILE and PHYLIB_DIR env vars, or place madRoot.R alongside this script.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript 13_root_trees_mad.R <input.nwk> <output_rooted.nwk>\n")
  quit(status = 1)
}

input_file  <- args[1]
output_file <- args[2]

library("inline")
library("Rcpp")

# ── Locate madRoot.R ──────────────────────────────────────────────────────
SCRIPT_DIR <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
  grep("--file=", commandArgs(trailingOnly = FALSE))
] |> sub("--file=", "", x = _), mustWork = FALSE))

# Check environment variable, then fallback to resources/ relative to pipeline root
MAD_ROOT_FILE <- Sys.getenv("MAD_ROOT_FILE", "")
if (MAD_ROOT_FILE == "" || !file.exists(MAD_ROOT_FILE)) {
  # Try alongside this script
  MAD_ROOT_FILE <- file.path(SCRIPT_DIR, "madRoot.R")
}
if (!file.exists(MAD_ROOT_FILE)) {
  # Try resources/ in pipeline root
  pipeline_root <- dirname(dirname(SCRIPT_DIR))
  MAD_ROOT_FILE <- file.path(pipeline_root, "resources", "madRoot.R")
}

PHYLIB_DIR <- Sys.getenv("PHYLIB_DIR", "")
if (PHYLIB_DIR == "" || !dir.exists(PHYLIB_DIR)) {
  PHYLIB_DIR <- dirname(MAD_ROOT_FILE)
}

if (!file.exists(MAD_ROOT_FILE)) {
  cat("ERROR: madRoot.R not found. Set MAD_ROOT_FILE env var.\n")
  quit(status = 1)
}

# ── Source madRoot (needs working dir for Phylib includes) ────────────────
old_wd <- getwd()
setwd(PHYLIB_DIR)
source(MAD_ROOT_FILE)
setwd(old_wd)

# ── Root the tree ─────────────────────────────────────────────────────────
tryCatch({
  tree_content  <- readLines(input_file, warn = FALSE)
  rooted_trees  <- madRoot(tree_content)
  writeLines(rooted_trees, output_file)
  cat("OK\n")
  quit(status = 0)
}, error = function(e) {
  cat(sprintf("ERROR: %s\n", e$message))
  quit(status = 1)
})
