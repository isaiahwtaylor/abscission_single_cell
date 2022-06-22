library(Seurat)
library(here)
library(future)
#library(BiocParallel)
options(future.globals.maxSize = 100000 * 1024^2)

rice_seu = readRDS(here("../data/", "Rice_Atlas_8WT_ref_tz2_seu4.rds"))
Idents(object = rice_seu) = "consensus.anno"
ident_vec = unique(rice_seu@meta.data$consensus.anno)

cluster_markers = list()
count = 1
for a in ident_vec {
    cluster_marker[count] <- FindMarkers(object = pbmc, ident.1 = a, logfc.threshold = 2, max.cells.per.ident = 1000)
    saveRDS(Markers, here("../data/", paste(a, "_markers.rds")))
}

# markers = FindAllMarkers(rice_seu, max.cells.per.ident = 1000, logfc.threshold = 2)

# saveRDS(Markers, here("../data/", "markers.rds"))


# ident_vec = unique(rice_seu@meta.data$consensus.anno)
# M=length(ident_vec)
# N=M-1

# FindMarker.wrapper <- function(x){
#   FindMarkers(rice_seu, ident.1=ident_vec[x+1], only.pos = TRUE, min.pct=1)
# }

# Markers <- bplapply(0:N, FindMarker.wrapper, BPPARAM=MulticoreParam(16))

