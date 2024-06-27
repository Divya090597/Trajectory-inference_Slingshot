###################### PSEUDO BULK ######################

######### 1.Done some Data Manupulation

View(Bulk_agg@meta.data)

# Aggregate counts by cell type for pseudobulk analysis
Bulk_agg <- AggregateExpression(data,
                                group.by =c("Cell_labels","Sample_name"),
                                assays = "RNA",
                                slot = "counts",
                                return.seurat = F)
View(data@meta.data)
DefaultAssay(data)

Bulk_agg <- Bulk_agg$RNA

## transpose
Bulk_agg.t <- t(Bulk_agg)

# convert to data.frame
Bulk_agg.t <- as.data.frame(Bulk_agg.t)
Bulk_agg.t[1:10,1:10]

# get values where to split
spltRows <- gsub('_.*','',rownames(Bulk_agg.t))
View(spltRows)
str(spltRows)

# split data.frame
cts_splt <- split.data.frame(Bulk_agg.t,
                             f = factor(spltRows))

cts_splt$`CMP CD41`[1:10,1:10]
cts_splt$`Cebpa control`[1:10,1:10]
cts_splt$`CMP Irf8-GFP+ MHCII+`[1:10,1:10]

# fix colnames and transpose
cts_mod <- lapply(cts_splt, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1',rownames(x))
  t(x)
})

gsub('.*_(.*)', '\\1', 'Cebpa control_GSM1873474')

cts_mod$`Cebpa KO`[1:10,1:10]

library(DESeq2)
counts_Cebpa <- cts_mod$`Cebpa control`



################## 2.Actual Pseudo Bulk ###################



Bulk_agg <- AggregateExpression(data,
                                group.by =c("Cell_labels","Sample_name"),
                                assays = "RNA",
                                slot = "counts",
                                return.seurat = T)

### View and Extract Data

View(Bulk_agg@assays$RNA$data)
col = colnames(Bulk_agg@assays$RNA$data)

###Parse Column Names and Create Metadata

meta = do.call(rbind, strsplit(col, "_"))
print(head(meta))
m_text = cbind(col, meta)
colnames(m_text) = c("ID", "CellType","Sample")
head(m_text)

### Create Sample_CellType Column

Sample_celltype = paste(m_text[, "Sample"], m_text[,"CellType"], sep="_")
m_text = cbind(m_text, Sample_celltype);
colnames(m_text)
head(m_text)

### Define Output Folder and File Paths
output_folder <- "~/DataAnalysis/R"
Project_Name <- "Bulk_D"
matrix_file = paste(output_folder, Project_Name, "Bulk_D.txt", sep="");

### Write Data to Files
write.table(2^ agg@assays$RNA$data, file = matrix_file, sep="\t", row.names = TRUE,col.names=NA)
head(Bulk_agg@assays$RNA$data)
metaData_file = paste(output_folder, Project_Name, "_meta.txt", sep="");
write.table(m_text, file = metaData_file, sep="\t", row.names=F)

### Write out the matrix ###
outputmatrix = matrix_file
outputmeta = metaData_file
