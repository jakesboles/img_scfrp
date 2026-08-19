library(tidyverse)
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(msigdbr)

setwd("/projects/b1169/boles/img_scfrp")

plots_dir <- "plots/wgcna/"
tab_dir <- "tab_data/wgcna/"

# Prep gene sets ----------------------------------------------------------

m_t2g <- msigdbr(species = "Homo sapiens",
                 category = "C2")
m_t2g <- m_t2g %>%
  filter(gs_subcollection %in% c("CP:BIOCARTA", "CP:KEGG_LEGACY", "CP:REACTOME", "CP:WIKIPATHWAYS", "CP:PID")) %>%
  dplyr::select(c(gs_name, gene_symbol))
# unique(m_t2g$gs_name)

go_t2g <- msigdbr(species = "Homo sapiens",
                  category = "C5")
go_t2g <- go_t2g %>%
  # filter(gs_subcollection %in% c("GO:BP", "GO:CC") %>%
  dplyr::select(c(gs_name, gene_symbol))

t2g <- rbind(go_t2g, m_t2g)

# Get module members & run analyses ------------------------------------------------------

df <- read.csv("tab_data/wgcna/module_members_consensus.csv")

background <- df$gene_name %>% unique()

modules <- unique(df$module)
modules <- modules[modules != "grey"]

for (j in seq_along(modules)){
  
  message(paste0(str_to_title(modules[j]), " module"))
  
  genes <- df %>%
    filter(module == modules[j]) %>%
    pull(gene_name)
  
  em <- enricher(genes,
                 TERM2GENE = t2g,
                 universe = background)
  
  if (em@result$p.adjust %>% min() < 0.05){
    p <- dotplot(em,
                 showCategory = 20,
                 x = "FoldEnrichment",
                 color = "p.adjust",
                 size = "GeneRatio")
    ggsave(p,
           filename = paste0(plots_dir, "path_enrich_consensus_", modules[j], ".png"),
           units = "in", dpi = 600,
           height = 10, width = 8)
  }
  
  write.csv(as.data.frame(em@result),
            file = paste0(tab_dir, "path_enrich_consensus_", modules[j], ".csv"))
}
