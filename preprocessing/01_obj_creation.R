library(Seurat)
library(scCustomize)
library(tidyverse)

proj_dir <- "/projects/b1169/boles/img_scfrp/"

data_dir <- paste0(proj_dir, "data/01_obj_creation/")
dir.create(data_dir,
           recursive = T,
           showWarnings = F)

create_object <- function(file, id){
  counts <- Read_CellBender_h5_Mat(file)
  return(CreateSeuratObject(counts = counts, 
                            project = id))
}

# Pull only probes included in kit
# CellBender finds other genes somehow? 
probes <- read.csv("/projects/b1042/Gate_Lab/boles/cellranger_references/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv",
                   skip = 5)
genes <- probes$probe_id %>%
  str_split_i(pattern = "[|]",
              i = 2) %>% 
  unique()

# Get CellBender directories and filtered .h5 files 
in_dir <- "/projects/b1042/Gate_Lab/boles/img_scfrp/cellbender"
cellbender_dirs <- list.dirs(in_dir,
                             recursive = F)[str_detect(list.dirs(in_dir, 
                                                                 recursive = F),
                                                       "out")]

h5 <- list()

for (i in seq_along(cellbender_dirs)){
  h5[[i]] <- list.files(cellbender_dirs[i],
                        full.names = T)[str_detect(list.files(cellbender_dirs[i]), "_filtered.h5")]
}

h5 <- list_c(h5)
id <- str_split_i(h5, "/", 9) %>%
  str_split_i("[.]", 1) %>%
  str_remove_all("_filtered")

# Assemble objects

obj_list <- list()

for (i in seq_along(id)){
  message(paste0("Making Seurat object for ", id[i]))
  
  obj_list[[i]] <- create_object(h5[i],
                                 id[i])
}

for (i in seq_along(id)){
  message(paste0("Removing non-probe list genes from ", id[i]))
  idx <- rownames(obj_list[[i]]) %in% genes
  print(table(idx))
  obj_list[[i]] <- obj_list[[i]][idx, ]
}

message("Merging objects into one")
obj <- Merge_Seurat_List(obj_list,
                         add.cell.ids = id)
obj

message("Saving object")
saveRDS(obj,
        file = paste0(data_dir, "obj.rds"))