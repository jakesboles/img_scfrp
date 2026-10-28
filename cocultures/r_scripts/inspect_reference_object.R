# One-off diagnostic, NOT a numbered pipeline stage. The coculture arm's
# original preprocessing scripts were lost, but a finished/cleaned Seurat
# object from that lost pipeline still exists. This reverse-engineers the
# pipeline's parameters from that object rather than from a human's
# memory of it -- Seurat logs every function that touches an object
# (method, every argument) in obj@commands, which is close to a
# parameter-level transcript of what was actually run (NormalizeData,
# FindVariableFeatures, RunPCA, IntegrateLayers, FindClusters, etc.).
# Combined with the meta.data schema/value tables and variable features,
# that's enough to reconstruct a 01_obj_creation.R-through-0N-equivalent
# pipeline for this arm without transferring the (likely multi-GB) object
# itself anywhere.
#
# Usage: edit OBJ_PATH below, run on Quest (module load R/4.4.0 first),
# then share the full contents of the resulting
# cocultures/r_scripts/object_diagnostic.txt.

suppressMessages({
  library(Seurat)
})

OBJ_PATH <- "CHANGE_ME.rds"

obj <- readRDS(OBJ_PATH)

sink("object_diagnostic.txt")

cat("=== Basic structure ===\n")
print(class(obj))
cat("dim: ")
print(dim(obj))

cat("\nAssays:\n")
print(Assays(obj))
for (a in Assays(obj)) {
  cat(sprintf("\nAssay '%s' layers:\n", a))
  print(Layers(obj[[a]]))
}

cat("\nReductions:\n")
print(Reductions(obj))
for (r in Reductions(obj)) {
  cat(sprintf("\nReduction '%s': dims = %s\n", r,
              paste(dim(Embeddings(obj, r)), collapse = " x ")))
}

cat("\nGraphs:\n")
print(names(obj@graphs))

cat("\nmisc slot contents (names only):\n")
print(names(obj@misc))

cat("\n=== meta.data schema ===\n")
str(obj@meta.data)

cat("\n=== meta.data column value tables/summaries ===\n")
for (col in colnames(obj@meta.data)) {
  cat(sprintf("\n--- %s ---\n", col))
  x <- obj@meta.data[[col]]
  if (is.factor(x) || is.character(x) || is.logical(x)) {
    print(table(x, useNA = "ifany"))
  } else {
    print(summary(x))
  }
}

cat("\n=== Variable features ===\n")
cat(sprintf("Total: %d\n", length(VariableFeatures(obj))))
print(head(VariableFeatures(obj), 50))

cat("\n=== Command history (obj@commands) ===\n")
cat("This records every Seurat function run on this object, in order,\n")
cat("with every argument it was called with -- the main thing this\n")
cat("script exists to surface.\n")
print(names(obj@commands))
for (cmd_name in names(obj@commands)) {
  cat(sprintf("\n--- %s ---\n", cmd_name))
  print(obj@commands[[cmd_name]])
}

sink()

cat("Wrote object_diagnostic.txt in the current directory -- share its full contents.\n")
