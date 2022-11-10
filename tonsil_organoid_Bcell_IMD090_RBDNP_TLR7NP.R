library(Seurat)
library(dplyr)
library(ggplot2)
library(xlsx)
library(ggrepel)
library("org.Hs.eg.db")
hs <- org.Hs.eg.db
library("readxl")
setwd('/Volumes/GoogleDrive/My\ Drive/Stanford/RNA-seq/data_analysis/I53_NP_tonsil_Rhapsody_120d')

######################### data input ##############################################
tonsil_Rhapsody_table <- read.table("tonsil_organoid_I53_TLR7NP_processed_matrix.tsv", sep = '\t', row.names = 1, header = TRUE,check.names = FALSE)
tonsil_Rhapsody_tag_table <- read.table("tonsil_organoid_I53_TLR7NP_index_label.tsv", sep = '\t', row.names = 1, header = TRUE,check.names = FALSE)
p <- match(rownames(tonsil_Rhapsody_table),rownames(tonsil_Rhapsody_tag_table))
tonsil_Rhapsody_tag_table <- tonsil_Rhapsody_tag_table[p,,drop = FALSE]

all.markers <- colnames(tonsil_Rhapsody_table)
gene_list <- all.markers[!grepl('^ab-',all.markers)]
Abs_list <- all.markers[grepl('^ab-',all.markers)]

tonsil_Rhapsody_gene_table <- tonsil_Rhapsody_table[,gene_list]
tonsil_Rhapsody_Abs_table <- tonsil_Rhapsody_table[,Abs_list]
###### Seurat #######################################
# Here the seurat has columns as cell index and rows as genes
Bcell.subset <- CreateSeuratObject(counts = t(tonsil_Rhapsody_table), project = "tonsil_Rhapsody")
Bcell.subset[["antibody"]] <- CreateAssayObject(counts = t(tonsil_Rhapsody_Abs_table))
Bcell.subset[["gene"]] <- CreateAssayObject(counts = t(tonsil_Rhapsody_gene_table))

Bcell.subset@meta.data <- cbind(Bcell.subset@meta.data,tonsil_Rhapsody_tag_table)

ordered_condition <- c("day4 unstim","day4 RBD-NP" ,"day4 RBD-NP+TLR7NP","day9 unstim","day9 RBD-NP","day9 RBD-NP+TLR7NP")
Bcell.subset$condition <- factor(Bcell.subset$condition, levels = ordered_condition)

# plot the # genes & # counts to double check the data quality control (QC)
VlnPlot(Bcell.subset, features = "nFeature_RNA", group.by = 'days',split.by = 'treatment',pt.size = 0.1) + theme(axis.text = element_text(size = 20))
dev.print(png, 'tonsil_Rhapsody_IMD090_cleaned_QCmetric_ngene_RNA.png',width = 300, height = 500)
VlnPlot(Bcell.subset, features = "nCount_RNA", group.by = 'days',split.by = 'treatment',pt.size = 0.1,log = TRUE) + NoLegend() + theme(axis.text = element_text(size = 20))
dev.print(png, 'tonsil_Rhapsody_IMD090_cleaned_QCmetric_ncount_RNA.png',width = 300, height = 500)
VlnPlot(Bcell.subset, features = "nFeature_antibody", group.by = 'days',split.by = 'treatment',pt.size = 0.1) + NoLegend() + theme(axis.text = element_text(size = 20))
dev.print(png, 'tonsil_Rhapsody_IMD090_cleaned_QCmetric_ngene_antibody.png',width = 300, height = 500)
VlnPlot(Bcell.subset, features = "nCount_antibody", group.by = 'days',split.by = 'treatment',pt.size = 0.1,log = TRUE) + NoLegend() + theme(axis.text = element_text(size = 20))
dev.print(png, 'tonsil_Rhapsody_IMD090_cleaned_QCmetric_ncount_antibody.png',width = 300, height = 500)
VlnPlot(Bcell.subset, features = "nFeature_gene", group.by = 'days',split.by = 'treatment',pt.size = 0.1) + NoLegend() + theme(axis.text = element_text(size = 20))
dev.print(png, 'tonsil_Rhapsody_IMD090_cleaned_QCmetric_ngene_gene.png',width = 300, height = 500)
VlnPlot(Bcell.subset, features = "nCount_gene", group.by = 'days',split.by = 'treatment',pt.size = 0.1,log = TRUE) + NoLegend() + theme(axis.text = element_text(size = 20))
dev.print(png, 'tonsil_Rhapsody_IMD090_cleaned_QCmetric_ncount_gene.png',width = 300, height = 500)

log1p_assay <- CreateAssayObject(counts = log1p(t(tonsil_Rhapsody_table)))
Bcell.subset[["log1p"]] <- log1p_assay
DefaultAssay(Bcell.subset) <- "log1p"


####### B cell clustering #####################################################
Bcell.subset <- ScaleData(Bcell.subset)
Bcell.subset <- RunPCA(Bcell.subset, features = all.markers, npcs = 50)
# ElbowPlot(object = Bcell.subset,ndims = 50) + theme(axis.text = element_text(size = 20))
# dev.print(png, 'tonsil_Rhapsody_IMD090_Bcell_log1p_scale_PC_elbow.png',width = 700, height = 434)

n_pca_selected <- 30
Bcell.subset <- RunUMAP(Bcell.subset, reduction = "pca", dims = 1:n_pca_selected)
Bcell.subset <- FindNeighbors(Bcell.subset, reduction = "pca", dims = 1:n_pca_selected)
Bcell.subset <- FindClusters(Bcell.subset, resolution = 0.5)
# Bcell.subset$subclusters <- Bcell.subset$seurat_clusters
Bcell.subset$clusters <- Bcell.subset$seurat_clusters
DimPlot(Bcell.subset, label = TRUE, reduction = "umap")
dev.print(pdf, 'tonsil_Rhapsody_IMD090_Bcell_log1p_scale_umap_clusters.pdf',width = 5, height = 4)

DimPlot(Bcell.subset, label = TRUE, reduction = "umap",split.by = 'condition',ncol = 3)
dev.print(pdf, 'tonsil_Rhapsody_IMD090_Bcell_log1p_scale_umap_clusters_condition.pdf',width = 6, height = 4)

####### Define PB population ##############################
Bcell.markers <- FindAllMarkers(Bcell.subset)
Bcell.markers.top10 <- Bcell.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Bcell.subset, features = c(Bcell.markers.top10$gene)) +
  theme(text = element_text(size=20))
dev.print(png, 'tonsil_Rhapsody_IMD090_Bcell_log1p_scale_cluster_gene_heatmap.png',width = 1200, height = 2000)
write.xlsx(Bcell.markers, "tonsil_Rhapsody_IMD090_Bcell_log1p_scale_cluster_gene.xlsx")

# gene_name <- 'XBP1'
# View(Bcell.markers %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))
# c((Bcell.markers %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))['cluster'])
# View(nonPB.markers %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))
# c((nonPB.markers %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))['cluster'])
# FeaturePlot(Bcell.subset, features = c(gene_name),pt.size = 0.2, ncol = 1, sort.cell = TRUE,min.cutoff = 0,max.cutoff = 3)
# dev.print(pdf, paste('tonsil_Rhapsody_IMD090_Bcell_log1p_scale_umap_',gene_name,'.pdf',sep = ''),width = 5, height = 4)
# # BCL6: 4
# # ab-CD38: 4 1
# # ab-CD27: 2 0
# 
# Bcell.subset$if_PB <- 'nonPB'
# Bcell.subset$if_PB[Bcell.subset$clusters %in% c(5)] <- 'PB'
# DimPlot(Bcell.subset, label = TRUE, reduction = "umap",group.by = 'if_PB')
# 
# nonPB.subset <- subset(Bcell.subset,if_PB == 'nonPB')
# nonPB.markers <- FindAllMarkers(nonPB.subset)
# nonPB.markers.top10 <- nonPB.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# DoHeatmap(nonPB.subset, features = c(nonPB.markers.top10$gene)) +
#   theme(text = element_text(size=20))
# dev.print(png, 'tonsil_Rhapsody_IMD090_Bcell_log1p_scale_nonPB_subcluster_gene_heatmap.png',width = 1200, height = 2000)
# write.xlsx(nonPB.markers, "tonsil_Rhapsody_IMD090_Bcell_log1p_scale_nonPB_subcluster_gene.xlsx")

Bcell.subset$cluster_manual_label <- 'Naive B cells'
Bcell.subset$cluster_manual_label[Bcell.subset$clusters %in% c(4)] <- 'GC B cells'
Bcell.subset$cluster_manual_label[Bcell.subset$clusters %in% c(1)] <- 'Pre-GC'
Bcell.subset$cluster_manual_label[Bcell.subset$clusters %in% c(0)] <- 'Memory B cells'
Bcell.subset$cluster_manual_label[Bcell.subset$clusters %in% c(5)] <- 'PB'
Bcell.subset$cluster_manual_label[Bcell.subset$clusters %in% c(2)] <- 'CD11c+ memory B cells'
Bcell.subset$cluster_manual_label[Bcell.subset$clusters %in% c(6)] <- 'noise'

DimPlot(Bcell.subset, label = TRUE, group.by = 'cluster_manual_label',repel = T, reduction = "umap")
dev.print(pdf, 'tonsil_Rhapsody_IMD090_Bcell_log1p_scale_umap_cluster_labeled.pdf',width = 6.5, height = 4)

############### B cell subset visualization ##########################################
Bcell_cluster_condition_count <- as.data.frame.matrix(table(Bcell.subset$condition,Bcell.subset$cluster_manual_label))
Bcell_condition_count <- data.frame(table(Bcell.subset$condition),row.names = 1)
Bcell_cluster_condition_ratio <- floor((Bcell_cluster_condition_count/Bcell_condition_count$Freq)*1000)/10
# write.xlsx(Bcell_cluster_condition_ratio, "tonsil_Rhapsody_Bcell_log1p_scale_cluster_ratio.xlsx")
# write.xlsx(Bcell_cluster_condition_count, "tonsil_Rhapsody_Bcell_log1p_scale_cluster_count.xlsx")

Bcell_cluster_condition_ratio <- Bcell_cluster_condition_ratio[ordered_condition,]
Bcell_cluster_condition_ratio$condition <- factor(rownames(Bcell_cluster_condition_ratio),levels = ordered_condition)
Bcell_cluster_condition_ratio$treatment <- 'RBD-NP'
Bcell_cluster_condition_ratio$treatment[grepl('unstim',Bcell_cluster_condition_ratio$condition)] <- 'unstim'
Bcell_cluster_condition_ratio$treatment[grepl('TLR7NP',Bcell_cluster_condition_ratio$condition)] <- 'RBD-NP+TLR7NP'
Bcell_cluster_condition_ratio$treatment <- factor(Bcell_cluster_condition_ratio$treatment,levels = c('unstim','RBD-NP','RBD-NP+TLR7NP'))

for (cluster_name in unique(Bcell.subset$cluster_manual_label))
{
  print(cluster_name)
  graphics.off()
  plot <- ggplot(Bcell_cluster_condition_ratio,aes(x=condition, y=!!sym(cluster_name),fill = treatment)) + ggtitle(paste(cluster_name)) + ylab(paste(cluster_name,'(%)')) + geom_bar(stat="identity", width=0.5) + theme(text = element_text(size = 20),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis() + scale_fill_manual(values=c('grey','blue','red'))
  print(plot)
  dev.print(pdf, paste('tonsil_Rhapsody_IMD090_Bcell_log1p_scale_',cluster_name,'_condition_ratio.pdf',sep = ''),width = 5, height = 5)
}

library("readxl")
raw_cell_count <- read_excel("./qianyin1 12092020-sorting raw file/FACS_sorting_table.xlsx")
raw_cell_count <- raw_cell_count[raw_cell_count$sample %in% c('2','3','4','5','6','8'),]
raw_cell_count$condition <- ordered_condition
raw_cell_count <- cbind(raw_cell_count,Bcell_cluster_condition_ratio)
for (cluster_name in sort(as.character(unique(Bcell.subset$cluster_manual_label))))
{
  raw_cell_count[paste(cluster_name,'count')] <- raw_cell_count[cluster_name]*raw_cell_count$`B cell estimated`/100
}
write.xlsx(raw_cell_count, "tonsil_Rhapsody_Bcell_log1p_scale_cluster_count_estimated_FACS.xlsx")
raw_cell_count$condition <- factor(raw_cell_count$condition,levels = ordered_condition)
# plot <- ggplot(raw_cell_count,aes(x=condition, y=!!sym('B cell estimated'),fill = treatment)) + ggtitle(paste('#B cells estimated',sep = '')) + geom_bar(stat="identity", width=0.5) + theme(text = element_text(size = 20),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis() + scale_fill_manual(values=c('blue','red'))
# print(plot)
# dev.print(pdf, paste('tonsil_Rhapsody_IMD090_Bcell_log1p_scale_condition_count_estimated_FACS.pdf',sep = ''),width = 5, height = 5)
for (cluster_name in colnames(raw_cell_count)[grepl('count',colnames(raw_cell_count))])
{
  print(cluster_name)
  graphics.off()
  plot <- ggplot(raw_cell_count,aes(x=condition, y=!!sym(cluster_name),fill = treatment)) + ggtitle(paste('#',gsub(' count','',cluster_name))) + ylab('# cells') + geom_bar(stat="identity", width=0.5) + theme(text = element_text(size = 20),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis() + scale_fill_manual(values=c('grey','blue','red'))
  # plot <- ggplot(raw_cell_count,aes(x=condition, y=!!sym(cluster_name),fill = treatment)) + ggtitle(paste('#cluster',gsub(' count','',cluster_name),'B cells')) + ylab(paste('#cluster',cluster_name,'B cells')) + geom_bar(stat="identity", width=0.5) + theme(text = element_text(size = 20),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis() + scale_fill_manual(values=c('blue','red'))
  print(plot)
  dev.print(pdf, paste('tonsil_Rhapsody_IMD090_Bcell_log1p_scale_',cluster_name,'_condition_count.pdf',sep = ''),width = 7, height = 7)
}


#################PB#################################
PB.subset <- subset(Bcell.subset,cluster_manual_label == 'PB')
PB_secreted <- data.frame(t(as.matrix(PB.subset@assays$RNA@data)),check.names = F)[,all.markers[grepl('-secreted',all.markers)]]
PB_secreted$max <- colnames(PB_secreted)[apply(PB_secreted,1,which.max)]
PB_secreted$isotype <- PB_secreted$max
PB_secreted$isotype[grepl('IGHG',PB_secreted$max)] <- 'IgG'
PB_secreted$isotype[grepl('IGHM',PB_secreted$max)] <- 'IgM'
PB_secreted$isotype[grepl('IGHA',PB_secreted$max)] <- 'IgA'
PB.subset$isotype <- PB_secreted$isotype


# within B cells
PB_cluster_condition_count <- as.data.frame.matrix(table(PB.subset$condition,PB.subset$isotype))
Bcell_condition_count <- data.frame(table(Bcell.subset$condition),row.names = 1)
PB_cluster_condition_ratio <- floor((PB_cluster_condition_count/Bcell_condition_count$Freq)*1000)/10

PB_cluster_condition_ratio$condition <- factor(rownames(PB_cluster_condition_ratio),levels = ordered_condition)
PB_cluster_condition_ratio$treatment <- 'RBD-NP'
PB_cluster_condition_ratio$treatment[grepl('unstim',PB_cluster_condition_ratio$condition)] <- 'unstim'
PB_cluster_condition_ratio$treatment[grepl('TLR7NP',PB_cluster_condition_ratio$condition)] <- 'RBD-NP+TLR7NP'
PB_cluster_condition_ratio$treatment <- factor(PB_cluster_condition_ratio$treatment,levels = c('unstim','RBD-NP','RBD-NP+TLR7NP'))

for (cluster_name in unique(PB.subset$isotype))
{
  print(cluster_name)
  graphics.off()
  plot <- ggplot(PB_cluster_condition_ratio,aes(x=condition, y=!!sym(cluster_name),fill = treatment)) + ggtitle(paste(cluster_name,' PB (%)',sep = '')) + geom_bar(stat="identity", width=0.5) + ylab('% in B cells') + theme(text = element_text(size = 20),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis() + scale_fill_manual(values=c('grey','blue','red'))
  print(plot)
  dev.print(pdf, paste('tonsil_Rhapsody_IMD090_Bcell_log1p_scale_',cluster_name,'_condition_ratio.pdf',sep = ''),width = 7, height = 7)
}


################ DE analysis ##############################################
library(ggrepel)
logFC_cutoff <- 0.25
p_cutoff <- 0.05
#Seurat findmarker volcano

DefaultAssay(Bcell.subset) <- "log1p"

for (day_name in unique(Bcell.subset$days))
{
print(day_name)
cluster_subset <- subset(Bcell.subset, days == day_name)
  if ('RBD-NP+TLR7NP' %in% unique(cluster_subset$treatment) & 'RBD-NP' %in% unique(cluster_subset$treatment)){
    Bcell_LAIV_ns.markers <- FindMarkers(cluster_subset, group.by = "treatment",ident.1 = 'RBD-NP+TLR7NP',ident.2 = 'RBD-NP',logfc.threshold = 0,min.pct = 0)
    Bcell_LAIV_ns.markers$significant <- 'no'
    Bcell_LAIV_ns.markers$significant[((Bcell_LAIV_ns.markers$avg_log2FC <= -logFC_cutoff) | (Bcell_LAIV_ns.markers$avg_log2FC >= logFC_cutoff)) & (Bcell_LAIV_ns.markers$p_val_adj <= p_cutoff)] <- 'yes'
    Bcell_LAIV_ns.markers$gene <- rownames(Bcell_LAIV_ns.markers)
    geneID <- select(hs, 
                     keys = toupper(Bcell_LAIV_ns.markers$gene),#selected_gene_set1_org$gene[!grepl('ab',selected_gene_set1_org$gene)],#
                     columns = c("ENTREZID", "SYMBOL"),
                     keytype = "SYMBOL")
    Bcell_LAIV_ns.markers$ENTREZID <- geneID$ENTREZID
    Bcell_LAIV_ns.markers$labels <- rownames(Bcell_LAIV_ns.markers)
    Bcell_LAIV_ns.markers$labels[Bcell_LAIV_ns.markers$significant == 'no'] <- ''
    Bcell_LAIV_ns.markers$abs_pct_diff <- abs(Bcell_LAIV_ns.markers$pct.1 - Bcell_LAIV_ns.markers$pct.2)
    Bcell_LAIV_ns.markers$abs_avg_log2FC <- abs(Bcell_LAIV_ns.markers$avg_log2FC)
    Bcell_LAIV_ns.markers <- Bcell_LAIV_ns.markers %>% arrange(desc(significant),desc(avg_log2FC))
    write.xlsx(Bcell_LAIV_ns.markers,paste('tonsil_Rhapsody_IMD090_Bcell_RBDNP+TLR7NP_vs_RBDNP_',day_name,'_log1p.xlsx',sep = ''))
    graphics.off()
    plot<-ggplot(Bcell_LAIV_ns.markers, aes(x=avg_log2FC, y=-log10(p_val_adj),color = significant,legend = significant)) +
      geom_point() +
      scale_colour_manual(values = c("yes" = "red", "no" = "darkgrey")) +
      ggtitle(paste('DE RBD-NP+TLR7NP vs RBD-NP on',day_name)) +
      xlab("log2 of Fold Change") + ylab("-log10 adjusted p-value") +
      theme(text = element_text(size = 20)) +
      geom_text_repel(aes(x=avg_log2FC, y=-log10(p_val_adj),label = labels), color = "black", size = 5) +
      theme(legend.position="none") +
      geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff), linetype="dotted") +
      geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")
    print(plot)
    dev.print(png, paste('tonsil_Rhapsody_IMD090_Bcell_RBDNP+TLR7NP_vs_RBD-NP_',day_name,'_log1p.png',sep = ''),width = 600, height = 434*1.2)
  }
}


##### DAVID ###############
# 1. download .txt file
# 2. open in excel and save as excel
# 3. change the sheet name to Sheet1
# 4. double check in David about the total upregualted gene number and put this in GeneRatio command

library("readxl")
table1 <- 'tonsil_organoid_120d_Bcell_RBD-NP+TLR7NP_vs_RBD-NP_day4_log1p_GO_BP_FAT.xlsx'
table1_DE_data <- read_excel(table1, sheet = "Sheet1")
table1_DE_data$GeneRatio <- paste(table1_DE_data$Count,122,sep = '/')
n_present <- 12#sum(table1_DE_data$FDR <= 0.05)
# sapply(table1_DE_data, class)
table1_DE_data$GO_name <- gsub('.*~','',table1_DE_data$Term)
ego@result$Description[1:n_present] <- table1_DE_data$GO_name[1:n_present]
ego@result$GeneRatio[1:n_present] <- table1_DE_data$GeneRatio[1:n_present]
ego@result$p.adjust[1:n_present] <- table1_DE_data$FDR[1:n_present]
ego@result$Count[1:n_present] <- table1_DE_data$Count[1:n_present]
ego@result <- ego@result[1:n_present,]
dotplot(ego, showCategory=n_present) + ggtitle('GO of Bcell RBD-NP+TLR7NP vs RBD-NP\nupregulated on day4') + 
  theme(axis.text = element_text(size = 20),plot.title = element_text(size = 25, face = "bold"))
dev.print(pdf, paste('tonsil_Rhapsody_IMD090_Bcell_RBD-NP+TLR7NP_vs_RBD-NP_day4_log1p_upregulation_DAVID_GO_BP_FAT_top12.pdf',sep = ''),width = 8, height = 5)


