# Computes each coculture pool's sequencing depth (reads, mean reads/cell,
# cells) from Cell Ranger's qc_library_metrics.csv, to see how much depth
# varies across pools before deciding whether a seqtk downsampling step
# (like monoculture's) is actually needed for this arm. Coculture-arm
# counterpart of monocultures/seqtk/get_seqtk_proportions.R -- same core
# computation (reads/rpc/cells per pool; target = the lowest-depth pool;
# target_fraction = target/rpc for everything else), trimmed to only the
# packages this actually uses -- the monoculture version's ~40-library
# import (Seurat, DESeq2, MAST, ComplexHeatmap, ...) is mostly unrelated
# to this computation and just slows startup; not fixed there since this
# is a new file, but worth porting that fix back if that script is ever
# revisited.
#
# This is diagnostic only for now -- the user has decided coculture does
# not need a downsampling step yet -- so this writes both the
# target_fraction table (in case that decision changes later, same format
# monoculture's seqtk.sh would consume) AND a plot, so the actual
# magnitude of the depth differences is visible at a glance rather than
# needing interpretation from raw fractions.
#
# Caveat worth remembering (see CLAUDE.md's coculture section): depth
# here is necessarily per-pool, not per-batch -- qc_library_metrics.csv
# is a whole-GEM-well metric, and JSB134-3 holds samples from two
# different batches (_2 and _3) sequenced together in that one well, so
# those two batches' depth can't be distinguished from this file. The
# other 3 pools each map to exactly one batch, so per-pool depth doubles
# as per-batch depth for those.

suppressMessages({
  library(tidyverse)
})

setwd("/projects/b1169/boles/img_scfrp/cocultures")

pools <- list.dirs("cellranger",
                    recursive = F,
                    full.names = F)

pools <- pools[str_detect(pools, "test", negate = T)]

numbers <- list()

for (pool in pools) {

  qc <- read.csv(paste0("cellranger/", pool, "/outs/qc_library_metrics.csv"))

  reads <- qc %>%
    filter(Metric.Name == "Number of reads" &
             Grouped.By == "Fastq ID") %>%
    pull(Metric.Value) %>%
    as.numeric() %>%
    sum()

  rpc <- qc %>%
    filter(Metric.Name == "Mean reads per cell") %>%
    pull(Metric.Value) %>%
    as.numeric()

  cells <- qc %>%
    filter(Metric.Name == "Cells") %>%
    pull(Metric.Value) %>%
    as.numeric()

  numbers[[pool]] <- tibble(
    reads = reads,
    rpc = rpc,
    cells = cells,
    pool = pool
  )

}

numbers <- list_rbind(numbers)

numbers$target <- min(numbers$rpc)

numbers <- numbers %>%
  mutate(target_fraction = target / rpc)

write.csv(numbers, "seqtk/depth_by_pool.csv",
          row.names = F, quote = F)

# Same table format monoculture's seqtk.sh consumes, in case this arm
# ends up needing downsampling after all -- pools already at the target
# depth (target_fraction == 1) are dropped, same as monoculture's
# version.
targets <- numbers %>%
  filter(target_fraction != 1) %>%
  dplyr::select(c(pool, target_fraction))

write.table(targets, "seqtk/downsample_targets.txt",
            row.names = F,
            col.names = F,
            quote = F,
            sep = ",")

# The actual "how bad is it" view -- monoculture's script never needed
# one of these (it only ever fed the raw fractions to seqtk), but a bare
# target_fraction table doesn't convey magnitude at a glance the way a
# bar chart does.
p <- ggplot(numbers, aes(x = pool, y = rpc)) +
  geom_col(fill = "steelblue") +
  geom_hline(aes(yintercept = target), linetype = "dashed", color = "red") +
  geom_text(aes(label = round(rpc)), vjust = -0.5) +
  labs(
    x = "Pool",
    y = "Mean reads per cell",
    title = "Coculture sequencing depth by pool",
    subtitle = "Dashed line = lowest-depth pool (the seqtk downsampling target, if used)"
  ) +
  theme_minimal()

ggsave(p, filename = "seqtk/reads_per_cell_by_pool.png",
       units = "in", dpi = 300, width = 7, height = 5)
