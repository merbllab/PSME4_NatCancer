library(stringr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(hexbin)
library(viridis)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(wesanderson)
library(corrplot)
library(Hmisc)
library(ggrepel)
# path =  "D:\\Dropbox (Weizmann Institute)\\Aaron Javitt\\PSME4_PaperDraft\\PSME4 Figures\\Figure2_MAPPQC\\"
path =  "C:\\Users\\ajavi\\Dropbox (Weizmann Institute)\\PSME4\\PSME4 Figures\\Figure2_MAPPQC\\"

Proteins_noImp = read.csv(paste0(path,"060219_1_Analysis\\060219_1_DataSnapshots\\060219_1_D_ProteinsAfterFilteringNonValidProteins.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)

diffProteinTable = Proteins_sub1[unlist(diffProteins),c(Tumor_col,Adj_col)]
rownames(diffProteinTable) = unlist(diffProteins)
rownames(names_frame) = names_frame$original
colnames(diffProteinTable) = names_frame[c(Tumor_col,Adj_col),"cleanColname"]
col_fun = colorRamp2(c(0,10), c("white", "darkred"))
diffProteinTable[is.na(diffProteinTable)] = 0
abun_col = colnames(Proteins_noImp)[198:213]
cleanColname = str_sub(abun_col,start = 39,end = -9)
names_frame  <-  as.data.frame(cleanColname)
names_frame$original = abun_col
names_frame$Sample = factor(str_sub(names_frame$cleanColname ,start = 1,end = 1))
names_frame$Sample = recode(names_frame$Sample,N = "Adjacent",T = "Tumor")
names_frame$ID = str_sub(names_frame$cleanColname,start = 3)
Tumor_col = names_frame[names_frame$Sample== "Tumor","original"]
Adj_col= names_frame[names_frame$Sample== "Adjacent","original"]
Proteins_sub = Proteins_noImp[,c("First.Gene.Name",Tumor_col,Adj_col)]
Proteins_sub$T_count = apply(Proteins_noImp[,Tumor_col], 1 ,function(x) 8 - sum(is.na(x)))
Proteins_sub$A_count = apply(Proteins_noImp[,Adj_col], 1 ,function(x) 8 - sum(is.na(x)))
Proteins_sub$dif_count = Proteins_sub$T_count - Proteins_sub$A_count
Proteins_sub1 = Proteins_sub[!str_detect(Proteins_sub$First.Gene.Name,"PSM") & !str_detect(Proteins_sub$First.Gene.Name,"IGH") & !str_detect(Proteins_sub$First.Gene.Name,"IGK"),]
Proteins_sub1 = Proteins_sub1[order(Proteins_sub1$dif_count),]
Proteins_sub1$First.Gene.Name = factor(Proteins_sub1$First.Gene.Name,levels = (Proteins_sub1[order(Proteins_sub1$dif_count),"First.Gene.Name"]))
Proteins_sub1$ind = 1:nrow(Proteins_sub1)

diffProteins = Proteins_sub1 %>% filter(abs(dif_count) >=5) %>% select(First.Gene.Name)
#TCGA
TCGA = read.csv(paste0(path,"TCGA_MAPPdifferentials.csv",sep = ""),row.names = 1)
TCGA_sub = TCGA %>% filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal"))
TCGA_annot = TCGA_sub %>% select(sample_type)
TCGA_numeric = TCGA_sub %>% select(CDH5:DENR) %>% t() %>% data.frame()
TCGA_numeric_norm = data.frame(t(scale(t(TCGA_numeric))))
TCGA_numeric_norm = TCGA_numeric_norm[as.character(diffProteins$First.Gene.Name),]
rownames(TCGA_numeric_norm) = as.character(diffProteins$First.Gene.Name)
TCGA_annot_map = HeatmapAnnotation(df = data.frame(Sample = TCGA_annot$sample_type),
                                col = list(Sample = c("Solid Tissue Normal" = "skyblue3","Primary Tumor" = "darkorange3")),
                                annotation_name_gp = element_blank())
DIFF_TCGA = Heatmap(as.matrix(TCGA_numeric_norm),
                     show_row_names = TRUE,
                     show_column_names = FALSE,
                     show_row_dend =  FALSE,
                     show_column_dend = FALSE,
                     cluster_columns = TRUE,
                     cluster_rows = TRUE,
                     name = "Normalized RPKM",
                     na_col = "lightgrey",
                     top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("skyblue3","darkorange3")))),
                     column_split  = TCGA_annot$sample_type,
                     row_names_gp = gpar(fontsize = 8),
                     column_title = "TCGA LUAD",
                     height = 8,
                     width = 10)
#CPTAC
ProtData = read.csv(paste(path,"alldata_publicproteomeLC.csv",sep = ""),row.names = 1)
SampData = read.csv(paste(path,"sampleData.csv",sep = ""))
select_col = str_detect(colnames(ProtData),"Unshared")
Prot_data_sub = ProtData[-3:-1,select_col]
select_col1 = (str_detect(colnames(Prot_data_sub),"CPT"))
Prot_data_sub1 = Prot_data_sub[,select_col1]
new_names = str_split(colnames(Prot_data_sub1),"\\.",simplify = TRUE)[,1]
colnames(Prot_data_sub1) = new_names
Prot_data_sub1_DIFF = Prot_data_sub1[as.character(diffProteins$First.Gene.Name),]
rownames(Prot_data_sub1_DIFF) = as.character(diffProteins$First.Gene.Name)
ColType_cptac = SampData %>% filter(Aliquot..Specimen.Label. %in% colnames(Prot_data_sub1_DIFF)) %>% select(Aliquot..Specimen.Label.,Type) %>% tibble::deframe()

CPTAC_annot = HeatmapAnnotation(df = data.frame(Sample = ColType_cptac),
                                 col = list(Sample = c("Normal" = "skyblue3","Tumor" = "darkorange3")),
                                annotation_name_gp = element_blank())
DIFF_CPTAC = Heatmap(as.matrix(Prot_data_sub1_DIFF),
                        show_row_names = TRUE,
                        show_column_names = FALSE,
                        show_row_dend =  FALSE,
                        show_column_dend = FALSE,
                        cluster_columns = TRUE,
                        cluster_rows = TRUE,
                        name = "CPTAC Protein Abundance",
                        na_col = "lightgrey",
                     top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("skyblue3","darkorange3")))),
                     column_split  = ColType_cptac,
                        row_names_gp = gpar(fontsize = 8),
                     column_title = "CPTAC Proteomics",
                     height = 8,
                     width = 10)
draw(DIFF_CPTAC)


#MAPP

MAPPdata = read.csv(paste0(path,"060219_1_Analysis\\060219_1_FilteredData\\060219_1_Lung-cancer-midgam_proteins.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
MAPPdata_sub = MAPPdata %>%select(First.Gene.Name,N.107.Log10.median:N.76.Log10.median,T.107.Log10.median:T.76.Log10.median)
rownames(MAPPdata_sub) = MAPPdata_sub$First.Gene.Name
MAPPdata_sub = MAPPdata_sub[,-1]
MAPPdata_sub_diff = MAPPdata_sub[as.character(diffProteins$First.Gene.Name),]
MAPPdata_sub_diff_norm = data.frame(t(scale(t(MAPPdata_sub_diff))))
colnames(MAPPdata_sub_diff_norm) = str_sub(colnames(MAPPdata_sub_diff_norm),end = -14)

sample_annot = HeatmapAnnotation(df = data.frame(Sample = c(rep("Adjacent",8),rep("Tumor",8))),
                                 col = list(Sample = c("Adjacent" = "skyblue3","Tumor" = "darkorange3")))
MAPP_diff = Heatmap(MAPPdata_sub_diff_norm,
                        show_row_names = TRUE,
                        show_column_names = FALSE,
                        show_row_dend =  FALSE,
                        show_column_dend = FALSE,
                        cluster_columns = TRUE,
                        cluster_rows = TRUE,
                        name = "MAPP intensity",
                        na_col = "lightgrey",
                    top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("skyblue3","darkorange3")))),
                    row_names_gp = gpar(fontsize = 8),
                        column_split  = c(rep("Adjacent",8),rep("Tumor",8)),
                    column_title = "MAPP",
                    height = 16,
                    width = 8)
## 1D


Protdata = read.csv(paste0(path,"110319_1_Analysis-1D\\110319_1_FilteredData\\110319_1_Lung-Midgam-1D_proteins.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
Protdata_sub = Protdata %>%select(First.Gene.Name,N.107.Log10.median:N.76.Log10.median,T.107.Log10.median:T.76.Log10.median) %>% filter(First.Gene.Name %in% as.character(diffProteins$First.Gene.Name))
rownames(Protdata_sub) = Protdata_sub$First.Gene.Name
Protdata_sub = Protdata_sub[,-1]
Protdata_sub_diff = Protdata_sub[as.character(diffProteins$First.Gene.Name),]
Protdata_sub_diff_norm = data.frame(t(scale(t(Protdata_sub_diff))))
colnames(Protdata_sub_diff_norm) = str_sub(colnames(Protdata_sub_diff_norm),end = -14)

Prot_annot = HeatmapAnnotation(df = data.frame(Sample = c(rep("Adjacent",8),rep("Tumor",8))),
                                 col = list(Sample = c("Adjacent" = "skyblue3","Tumor" = "darkorange3")),annotation_name_gp = element_blank())
Prot_diff = Heatmap(Protdata_sub_diff_norm,
                    show_row_names = TRUE,
                    show_column_names = FALSE,
                    show_row_dend =  FALSE,
                    show_column_dend = FALSE,
                    cluster_columns = TRUE,
                    cluster_rows = TRUE,
                    name = "WCE Proteomics Protein Abundance",
                    na_col = "lightgrey",
                    top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("skyblue3","darkorange3")))),
                    row_names_gp = gpar(fontsize = 8),
                    column_split  = c(rep("Adjacent",8),rep("Tumor",8)),
                    column_title = "WCE Proteomics",
                    height = 10,
                    width = 8)
pdf(paste0(path,"heatmapDifferentials_combined_TCGA.pdf"),width = 12) 
draw(MAPP_diff + Prot_diff + DIFF_CPTAC + DIFF_TCGA)
dev.off() 

#ttest_RANK
Protdata_ALL = Protdata %>%select(First.Gene.Name,N.107.Log10.median:N.76.Log10.median,T.107.Log10.median:T.76.Log10.median) 
Protdata_ALL[,2:17] = 10^Protdata_ALL[,2:17] 
rownames(Protdata_ALL) = make.unique(Protdata_ALL$First.Gene.Name)
TTEST_RESULTS = as.data.frame(matrix(nrow = nrow(Protdata_ALL),ncol = 4))
colnames(TTEST_RESULTS) = c("Gene","Mean-Adj","Mean-Tum","pval")
rownames(TTEST_RESULTS) = rownames(Protdata_ALL)
for (gene in rownames(Protdata_ALL)){
  TTEST_RESULTS[gene,"Gene"] = Protdata_ALL[gene,"First.Gene.Name"]
  TTEST_RESULTS[gene,"Mean-Adj"] = rowMeans(Protdata_ALL[gene,c(2:9)]) 
  TTEST_RESULTS[gene,"Mean-Tum"] = rowMeans(Protdata_ALL[gene,c(10:17)]) 
    test = t.test(x = Protdata_ALL[gene,c(2:9)],y = Protdata_ALL[gene,c(10:17)]) 
    TTEST_RESULTS[gene,"pval"] = test$p.value
}
TTEST_RESULTS$adjP = p.adjust(TTEST_RESULTS$pval,method = "BH")
TTEST_RESULTS$FC = log2(TTEST_RESULTS$`Mean-Tum`/TTEST_RESULTS$`Mean-Adj`)
TTEST_RESULTS$DIFF = "NO"
TTEST_RESULTS[TTEST_RESULTS$Gene %in% as.character(diffProteins$First.Gene.Name),"DIFF"] = "Differential MAPP"
TTEST_RESULTS$SIZE = ifelse(TTEST_RESULTS$DIFF == "Differential MAPP",.02,.01)
TTEST_RESULTS %>%
  arrange(SIZE) %>%
  ggplot() + 
  geom_point(aes(x = FC,y = -log10(adjP),color = DIFF))+
  geom_hline(yintercept = 1.3)+
  geom_vline(xintercept = -1) +
  geom_vline(xintercept =  1)+
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=11,face = "bold"),
        axis.title.x  = element_text(size=11,face = "bold"),
        legend.text = element_text(size=8),
        legend.position = "right",
        legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(colour = "black", fill=NA, size=.5))+
  scale_color_manual(values = c("darkred","lightgrey")) +
  xlab(paste("Log2 Fold Change Tumor/Adjacent")) + ylab("-Log10(Pvalue)")
ggsave(filename = paste0(path,"diffVolcano_MAPP.pdf"),height = 3)


write.csv(diffProteins,file = paste0(path,"MAPP_differentials.csv"))

