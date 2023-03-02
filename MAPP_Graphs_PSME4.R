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
library(ggbeeswarm)
# path =  "D:\\Dropbox (Weizmann Institute)\\PSME4\\PSME4 Figures\\Figure2_MAPPQC\\"
path =  "C:\\Users\\ajavi\\Dropbox (Weizmann Institute)\\PSME4\\PSME4 Figures\\Figure2_MAPPQC\\"
SampRename =  read.csv(paste0(path,"SampleRename.csv",sep = ""))
rownames(SampRename) = SampRename$OldName
#Peptide Number
Count =  read.csv(paste0(path,"060219_1_Analysis\\060219_1_Statistics_Summary\\060219_1_C_PeptidesDatasetGeneralStatistics.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
Count$Name = str_sub(Count$X,start = 35,end = -9)
Count$Sample = factor(str_sub(Count$Name,start = 1,end = 1))
Count$Sample = recode(Count$Sample,N = "Adjacent",T = "Tumor")
Count$ID = str_sub(Count$Name,start = 3,end = -1)
Count %>%
 filter(ID != "81") %>% #filter Sample81
  ggplot(aes(x = Sample,y = ValidValues)) +
  geom_dotplot(aes(fill = Sample,color = Sample),binaxis='y', stackdir='center', stackratio=1.5, dotsize=1,lwd = .1) + 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.4,lwd = .1) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom="errorbar", width=0.2,lwd = .1) +
   theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12,face = "bold"),
        axis.title.x  = element_text(size=12,face = "bold"),
        legend.text = element_text(size=8),
        legend.position = "right",
        legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
        panel.background = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.05),
        aspect.ratio = 1.5,
        panel.border = element_rect(colour = "black", fill=NA, size=.05),
        strip.text.x = element_text(size = 3)) +
  ylim(0,4000) +
  scale_fill_manual(values = c("lightskyblue3","darkorange3"),aesthetics = c("color","fill")) +
  stat_compare_means(paired = TRUE,label.x.npc = "middle",label = "p.signif") +
  xlab("") + ylab("Peptide Count")
ggsave(file = paste0(path,"peptideCount_001.pdf"),useDingbats=FALSE,height = 3)
#Protein Number
Count =  read.csv(paste0(path,"060219_1_Analysis\\060219_1_Statistics_Summary\\060219_1_C_ProteinsDatasetGeneralStatistics.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
Count$Name = str_sub(Count$X,start = 39,end = -9)
Count$Sample = factor(str_sub(Count$Name,start = 1,end = 1))
Count$Sample = recode(Count$Sample,N = "Adjacent",T = "Tumor")
Count$ID = str_sub(Count$Name,start = 3,end = -1)
Count %>%
  filter(ID != "81") %>% #filter Sample81
  ggplot(aes(x = Sample,y = ValidValues)) +
  geom_dotplot(aes(fill = Sample,color = Sample),binaxis='y', stackdir='center', stackratio=1.5, dotsize=1,lwd = .1) + 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.4,lwd = .1) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom="errorbar", width=0.2,lwd = .1) +
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12,face = "bold"),
        axis.title.x  = element_text(size=12,face = "bold"),
        legend.text = element_text(size=8),
        legend.position = "right",
        legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
        panel.background = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.05),
        aspect.ratio = 1.5,
        panel.border = element_rect(colour = "black", fill=NA, size=.05),
        strip.text.x = element_text(size = 3)) +
  scale_fill_manual(values = c("lightskyblue3","darkorange3"),aesthetics = c("color","fill")) +
  stat_compare_means(paired = TRUE,label.x.npc = "middle",label = "p.signif") +
  xlab("") + ylab("Protein Count")
ggsave(file = paste0(path,"proteinCount_001.pdf"),useDingbats=FALSE,height = 3)
#Protein Number_ proteome
ProCount =  read.csv(paste0(path,"110319_1_Analysis-1D\\110319_1_Statistics_Summary\\110319_1_C_ProteinsDatasetGeneralStatistics.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
ProCount$Name = str_sub(ProCount$X,start = 35,end = -9)
ProCount$Sample = factor(str_sub(ProCount$Name,start = 1,end = 1))
ProCount$Sample = recode(ProCount$Sample,N = "Adjacent",T = "Tumor")
ProCount$ID = str_sub(ProCount$Name,start = 3,end = -1)
ProCount %>%
  filter(ID != "81") %>% #filter Sample81
  ggplot(aes(x = Sample,y = Mean)) +
  geom_dotplot(aes(fill = Sample),binaxis='y', stackdir='center', stackratio=1.5, dotsize=1,lwd = 1.5) + 
  #geom_boxplot(aes(color = Sample)) +
  #geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2, position=position_dodge(.9)) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.4) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom="errorbar", width=0.2) +
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12,face = "bold"),
        axis.title.x  = element_text(size=12,face = "bold"),
        legend.text = element_text(size=8),
        legend.position = "right",
        legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
        panel.background = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.05),
        aspect.ratio = 1.5,
        panel.border = element_rect(colour = "black", fill=NA, size=.05),
        strip.text.x = element_text(size = 3)) +
  scale_fill_manual(values = c("lightskyblue3","darkorange3")) +
  stat_compare_means(paired = TRUE) +
  xlab("") + ylab("Mean Intensity Proteome")
ggsave(file = paste0(path,"ProteomeIntensity.pdf"),useDingbats=FALSE)

#PCA
PCA = read.csv(paste0(path,"060219_1_Analysis\\060219_1_Statistics_Summary\\060219_1_ProteinsPCAValues.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
PCA$X.1 = str_replace(str_sub(PCA$X.1,end = -3),"-",".")
PCA$NewName = SampRename[PCA$X.1,"NewName"]
PCA$Sample = factor(str_sub(PCA$NewName  ,start = 1,end = 1))
PCA$Sample = recode(PCA$Sample,A = "Adjacent",T = "Tumor")
PCA$ID = PCA$NewName
PCA %>%
  filter(ID != "T.9") %>% filter(ID != "A.9") %>% #filter Sample81
  ggplot(aes(x = X,y = Y)) +
  geom_point(aes(color = Sample),size = 4) +
  geom_text_repel(aes(label = ID),hjust = 0, nudge = 0.5) +
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12,face = "bold"),
        axis.title.x  = element_text(size=12,face = "bold"),
        legend.text = element_text(size=10),
        legend.position = "right",
        legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
        panel.background = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.05),
        aspect.ratio = 1,
        panel.border = element_rect(colour = "black", fill=NA, size=.05),
        strip.text.x = element_text(size = 3)) +
  scale_color_manual(values = c("lightskyblue3","darkorange3")) +
  xlab("PC1") + ylab("PC2") + ggtitle("MAPP Degradaome")
ggsave(file = paste0(path,"PCA_MAPP_001.pdf"),useDingbats = FALSE, height = 4)

#PCA 1D
PCA = read.csv(paste0(path,"110319_1_Analysis-1D\\110319_1_Statistics_Summary\\110319_1_ProteinsPCAValues.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
PCA$X.1 = str_replace(str_sub(PCA$X.1,end = -3),"-",".")
PCA$NewName = SampRename[PCA$X.1,"NewName"]
PCA$Sample = factor(str_sub(PCA$NewName  ,start = 1,end = 1))
PCA$Sample = recode(PCA$Sample,A = "Adjacent",T = "Tumor")
PCA$ID = PCA$NewName
PCA %>%
  filter(ID != "T.9") %>% filter(ID != "A.9") %>% #filter Sample81
  ggplot(aes(x = X,y = Y)) +
  geom_point(aes(color = Sample),size = 4) +
  geom_text_repel(aes(label = ID),hjust = 0, nudge = 0.5) +
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12,face = "bold"),
        axis.title.x  = element_text(size=12,face = "bold"),
        legend.text = element_text(size=10),
        legend.position = "right",
        legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
        panel.background = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.05),
        aspect.ratio = 1,
        panel.border = element_rect(colour = "black", fill=NA, size=.05),
        strip.text.x = element_text(size = 3)) +
  scale_color_manual(values = c("lightskyblue3","darkorange3")) +
  xlab("PC1") + ylab("PC2") + ggtitle("Proteome")
ggsave(file = paste0(path,"PCA_1D.pdf"), useDingbats = FALSE,height = 4)

#differential heatmap

Proteins_noImp = read.csv(paste0(path,"060219_1_Analysis\\060219_1_DataSnapshots\\060219_1_D_ProteinsAfterFilteringNonValidProteins.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
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
#Proteins_sub1 = Proteins_sub
Proteins_sub1 = Proteins_sub1[order(Proteins_sub1$dif_count),]
Proteins_sub1$First.Gene.Name = factor(Proteins_sub1$First.Gene.Name,levels = (Proteins_sub1[order(Proteins_sub1$dif_count),"First.Gene.Name"]))
Proteins_sub1$ind = 1:nrow(Proteins_sub1)
write.csv(Proteins_sub1,paste0(path,"MAPPprot_nocont.csv",sep = ""))
Proteins_sub1 %>%
  ggplot(aes(x = First.Gene.Name)) +
  geom_bar(aes(y = T_count),stat="identity", fill="darkorange3") +
  geom_bar(aes(y = -A_count),stat="identity", fill="lightskyblue3") +
  geom_point(aes(y = dif_count),color = "black") +
  geom_hline(yintercept  = 0) +
  theme(axis.text.y   = element_blank(),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12,face = "bold"),
        axis.title.x  = element_text(size=12,face = "bold"),
        legend.text = element_text(size=10),
        legend.position = "right",
        legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
        panel.background = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.05),
        aspect.ratio = 1,
        panel.border = element_rect(colour = "black", fill=NA, size=.05),
        strip.text.x = element_text(size = 3)) +
  xlab("Proteins") + ylab("Patient Count") + coord_flip() 
ggsave(file = paste0(path,"Differential_proteincounts.pdf"))

diffProteins = Proteins_sub1 %>% filter(abs(dif_count) >=5) %>% select(First.Gene.Name)

#differential_heatmap

diffProteinTable = Proteins_sub1[unlist(diffProteins),c(Tumor_col,Adj_col)]
rownames(diffProteinTable) = unlist(diffProteins)
rownames(names_frame) = names_frame$original
colnames(diffProteinTable) = names_frame[c(Tumor_col,Adj_col),"cleanColname"]
col_fun = colorRamp2(c(0,10), c("white", "darkred"))
diffProteinTable[is.na(diffProteinTable)] = 0

sample_annot = HeatmapAnnotation(df = data.frame(Sample = names_frame[c(Tumor_col,Adj_col),"Sample"]),
                                 col = list(Sample = c("Adjacent" = "skyblue3","Tumor" = "darkorange3")))
colnames(diffProteinTable) = SampRename[colnames(diffProteinTable) ,"NewName"]
Prot_diff_MAP = Heatmap(diffProteinTable,
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        show_row_dend =  FALSE,
                        show_column_dend = FALSE,
                        cluster_columns = TRUE,
                        cluster_rows = TRUE,
                        name = "MAPP intensity",
                        col = col_fun,na_col = "white",
                        top_annotation = sample_annot,
                        column_split = c(rep("Tumor",8),rep("Adjacent",8)),
                        row_names_gp = gpar(fontsize = 8),
                        height = 10,)
pdf(paste0(path,"heatmapDifferentials.pdf")) 
draw(Prot_diff_MAP)
dev.off() 
write.csv(diffProteinTable,paste0(path,"diffMAPPproteins.csv",sep = ""))
#Imp and Ratio calculate
impVAL = log10((10^(min(Proteins_sub1[,c(Tumor_col,Adj_col)],na.rm = TRUE)))/2)
MAPP_sub_imp = Proteins_sub1[,c("First.Gene.Name",Tumor_col,Adj_col)]
MAPP_sub_imp[is.na(MAPP_sub_imp)] = impVAL
names_frame$ID = as.factor(names_frame$ID)

for (IDX in levels(names_frame$ID)){
  variable = paste0(IDX,"_logFC_T-A")
  A_x = names_frame %>% filter(ID == IDX) %>% filter(Sample == "Adjacent") 
  T_x = (names_frame %>% filter(ID == IDX) %>% filter(Sample == "Tumor"))
  MAPP_sub_imp[,variable] = log2((10^MAPP_sub_imp[,(T_x$original)])/ (10^MAPP_sub_imp[,(A_x$original)]))
  
}
MAPP_sub_imp$Mean_Rat = rowMeans(MAPP_sub_imp[,c(17:24)])
MAPP_sub_imp$Med_Rat = apply(MAPP_sub_imp[,c(17:24)],1,median)
MAPP_sub_imp = MAPP_sub_imp %>% filter(Mean_Rat !=0)
rownames(MAPP_sub_imp) = MAPP_sub_imp$First.Gene.Name
for (Gene in MAPP_sub_imp$First.Gene.Name){
 MAPP_sub_imp[Gene,"pval"] =t.test(MAPP_sub_imp[Gene,Tumor_col],MAPP_sub_imp[Gene,Adj_col])$p.value

}
#localization
library(biomaRt)
library(hpar)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
hgnc_swissprot <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),mart = ensembl)
MAPP_ID = merge(MAPP_sub_imp,hgnc_swissprot,by.x = "First.Gene.Name",by.y = "hgnc_symbol")
MAPP_ID = unique(MAPP_ID)




location = getHpa(MAPP_ID$ensembl_gene_id, hpadata = "hpaSubcellularLoc")
location_MAPP = unique(merge(location,MAPP_ID,by.x = "Gene",by.y = "ensembl_gene_id" ))
location_MAPP[,"location"] = as.character(str_c(location_MAPP$Enhanced,location_MAPP$Supported,location_MAPP$Approved,sep = ";"))
location_MAPP$Nuclear = as.factor(ifelse(str_detect(location_MAPP$location,"Nuc"),"Nuclear","Other"))
location_MAPP$Mitochondria = as.factor(ifelse(str_detect(location_MAPP$location,"Mito"),"Mito","Not"))
location_MAPP$Cytosol = as.factor(ifelse(str_detect(location_MAPP$location,"Cyto") &! str_detect(location_MAPP$location,"Cytokin"),"Cyto","Not"))
location_MAPP$Membrane = as.factor(ifelse(str_detect(location_MAPP$location,"Plasma membrane")|str_detect(location_MAPP$location,"Junction"),"Mem","Not"))
location_MAPP$ERGOLG = as.factor(ifelse(str_detect(location_MAPP$location,"Endoplasmic")|str_detect(location_MAPP$location,"Golgi"),"ER","Not"))
location_MAPP$Cytoskeleton = as.factor(ifelse(str_detect(location_MAPP$location,"Actin")|str_detect(location_MAPP$location,"Microtubule") |str_detect(location_MAPP$location,"filament"),"Skeleton","Not"))
location_MAPP_sub = location_MAPP[location_MAPP$Reliability != "Uncertain",]
location_MAPP_sub = location_MAPP_sub[location_MAPP_sub$Reliability != "Approved",]
colors = brewer.pal(8,"Set2")[c(1,8)]
MAPP_sub_long_nuc = location_MAPP_sub %>% dplyr::select(Cytosol,Nuclear,`107_logFC_T-A`:`76_logFC_T-A`) %>% pivot_longer(cols = `107_logFC_T-A`:`76_logFC_T-A`, names_to = "Sample",values_to = "Ratio")
location_MAPP_sub %>% 
  filter(Mean_Rat != 0) %>%
  ggplot(aes(x=Nuclear,y=Mean_Rat)) +
  geom_quasirandom(aes(color = Nuclear),size = 0.5) +
  geom_boxplot(aes(), notch = TRUE,notchwidth = .7, width=0.1,size = 1,
               coef = 0,outlier.size = 0,outlier.stroke = 0,show.legend = FALSE) +
  theme(axis.text.y   = element_text(size=10,face = "bold"),
        axis.text.x   = element_text(size=10,face = "bold"),
        axis.title.y  = element_text(size=12,face = "bold"),
        axis.title.x  = element_text(size=12,face = "bold"),
        legend.text = element_text(size=11),
        legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_color_manual(values = c("#D37371","darkgrey"),aesthetics = c("fill","color")) +
  stat_compare_means() +
  ylim(-5.5,10) +
  xlab("") + ylab("FC Tumor/Normal") + labs(color = "Localization")
ggsave(paste0(path,"nucdeg_updated.pdf"),useDingbats=FALSE,height = 3)

location_MAPP_sub %>% 
#filter(Mean_Rat != 0) %>%
  ggplot(aes(x=Cytosol,y=Mean_Rat)) +
  geom_quasirandom(aes(color = Cytosol),size = 0.5) +
  geom_boxplot(aes(), notch = TRUE,notchwidth = .7, width=0.1,size = 1,
               coef = 0,outlier.size = 0,outlier.stroke = 0,show.legend = FALSE) +
  theme(axis.text.y   = element_text(size=10,face = "bold"),
        axis.text.x   = element_text(size=10,face = "bold"),
        axis.title.y  = element_text(size=12,face = "bold"),
        axis.title.x  = element_text(size=12,face = "bold"),
        legend.text = element_text(size=11),
        legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_color_manual(values = c("#f0ce7f","darkgrey"),aesthetics = c("fill","color")) +
  stat_compare_means() +
  ylim(-5.5,10) +
  xlab("") + ylab("FC Tumor/Normal") + labs(color = "Localization")
ggsave(paste0(path,"cytdeg_updated.pdf"),useDingbats=FALSE,height = 3)

#cleavage analysis
AA_properties = read.csv(paste0(path,"AA_properties.csv"), header = 1,stringsAsFactors = FALSE) 
peptides = read.csv(paste0(path,"060219_1_Analysis\\060219_1_FilteredData\\060219_1_Lung-cancer-midgam_peptides.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
#peptides_sub = peptides[!str_detect(peptides$First.Gene.Name,"PSM") & !str_detect(peptides$First.Gene.Name,"IGH") & !str_detect(peptides$First.Gene.Name,"IGK"),]
peptides_sub = peptides
peptide_IDcol = colnames(peptides)[64:79]
cleanColname_pep = str_sub(peptide_IDcol,start = 36,-3)
names_framepep  <-  as.data.frame(cleanColname_pep)
names_framepep$peptide = peptide_IDcol
names_frame = merge(names_frame,names_framepep,by.x = "cleanColname",by.y = "cleanColname_pep")

intensity_norm = peptides_sub[,names_frame$peptide]
intensity_norm[is.na(intensity_norm)] = 0
intensity_norm[intensity_norm >=1] = 1
intensity_norm = intensity_norm[,names_frame$peptide]
colnames(intensity_norm) = unlist(names_frame$cleanColname)
sums = colSums(intensity_norm,na.rm = TRUE)
intensity_norm$N = as.factor(peptides_sub$Amino.acid.before)
intensity_norm$C = as.factor(peptides_sub$Last.amino.acid)
c_term_agg = aggregate(. ~ C,intensity_norm[,-17], FUN=sum)
n_term_agg = aggregate(. ~ N,intensity_norm[,-18], FUN=sum)
AA_list_names = c_term_agg$C
rownames(c_term_agg) = AA_list_names
c_term_agg = c_term_agg[,-1]
rownames(n_term_agg) = n_term_agg$N
n_term_agg = n_term_agg[,-1]
c_term_agg_norm = sweep(c_term_agg, 2, sums, "/")
n_term_agg_norm = sweep(n_term_agg, 2, sums, "/")
c_term_agg_scale = t(scale(t(c_term_agg_norm)))
n_term_agg_scale = t(scale(t(n_term_agg_norm)))
n_term_agg_norm1 = log2(sweep(n_term_agg_norm, 1, rowMeans(n_term_agg_norm), "/"))
c_term_agg_norm1 = log2(sweep(c_term_agg_norm, 1, rowMeans(c_term_agg_norm), "/"))
rownames(c_term_agg_norm) = AA_list_names
col_fun = colorRamp2(c(-.7, 0, .7), c("darkblue", "white", "darkorange"))
col_fun_annot = colorRamp2(c(-4.5, 4.5), c("darkred",  "azure"))
AA_properties$Biochem = as.factor(AA_properties$Biochem)
AA_properties$Cleavage= as.factor(AA_properties$Cleavage)

rownames(AA_properties) = AA_properties$Letter
biochem_colors = brewer.pal(n = 5, name = 'Set3')
names(biochem_colors) = levels(AA_properties$Biochem)
cleavage_colors = c(brewer.pal(3,"Dark2")[c(1,3)],"white",brewer.pal(3,"Dark2")[c(2)])
names(cleavage_colors)= levels(AA_properties$Cleavage)
annotation_row_hmap = data.frame(Biochem = AA_properties[rownames(c_term_agg_scale),"Biochem"],
                                 Gravy_score =  AA_properties[rownames(c_term_agg_scale),"gravy"],
                                 Cleavage =  AA_properties[rownames(c_term_agg_scale),"Cleavage"])
annot_colors = list(Biochem = biochem_colors,
                    Gravy_score = col_fun_annot,
                    Cleavage = cleavage_colors)

MAP_annot = rowAnnotation(df = annotation_row_hmap,
                          col = annot_colors,
                          gp = gpar(col = "white", lwd = 0.1))
rownames(names_frame) = names_frame$cleanColname
sample_annot = HeatmapAnnotation(df = data.frame(Sample = names_frame[colnames(c_term_agg_norm1),"Sample"]),
                                 col = list(Sample = c("Adjacent" = "skyblue3","Tumor" = "darkorange3")))
colnames(c_term_agg_norm1) = SampRename[colnames(c_term_agg_norm1) ,"NewName"]
cterm = Heatmap(c_term_agg_norm1,
                show_row_names = TRUE,
                show_column_names = TRUE,
                show_row_dend =  FALSE,
                show_column_dend = FALSE,
                cluster_columns = FALSE,
                col = col_fun,
                rect_gp = gpar(col = "white", lwd = 0.5),
                height = unit(10,"cm"),
                name = "normalized percentage",
                column_title = "C Terminus Cleavage",
                right_annotation = MAP_annot,
                top_annotation =sample_annot)

pdf(paste0(path,"cTermCleavage_cleave.pdf")) 
draw(cterm)
dev.off() 

caspase = c("E","D")
tryptic = c("R","K")
chymotryptic = c("W","F","I","L","M","V","Y") 
c_term_agg_cleavage = as.data.frame(t(c_term_agg_norm1))
colnames(c_term_agg_cleavage) = AA_list_names
c_term_agg_cleavage$caspase = rowMeans(c_term_agg_cleavage[,caspase])
c_term_agg_cleavage$tryptic = rowMeans(c_term_agg_cleavage[,tryptic])
c_term_agg_cleavage$chymotryptic = rowMeans(c_term_agg_cleavage[,chymotryptic])

Proteome_proteins = read.csv(paste0(path,"110319_1_Analysis-1D\\110319_1_FilteredData\\110319_1_Lung-Midgam-1D_proteins.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
PSME4_1D = (Proteome_proteins[Proteome_proteins$First.Gene.Name == "PSME4",c(204:219)])
PSME4_1D_ratio = log2(PSME4_1D / rowMeans(PSME4_1D))
colnames(PSME4_1D_ratio) = str_sub(colnames(PSME4_1D_ratio),end = -8)
c_term_agg_cleavage$PSME4 = unlist(t(PSME4_1D_ratio))
# Proteasomeabun = read.csv(paste0(path,"220519_1_LC-MIDGAM-80per\\220519_1_Analysis\\220519_1_FilteredData\\220519_1_LC-MIDGAM-80per_proteins.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
# PSME4_Proteasome = (Proteasomeabun[Proteasomeabun$First.Gene.Name == "PSME4",c(214:221,223:230)])
# PSME4_Proteasome_ratio = log2(PSME4_Proteasome / rowMeans(PSME4_Proteasome))
# colnames(PSME4_Proteasome_ratio) = str_sub(colnames(PSME4_Proteasome_ratio),end = -8)
# c_term_agg_cleavage[,"PSME4_proteasome"] = PSME4_Proteasome_ratio %>% dplyr::select((rownames(c_term_agg_cleavage))) %>% t()





colnames(c_term_agg_cleavage)
enzyme = "chymotryptic"
c_term_agg_cleavage %>%
  ggplot(aes(x = PSME4,y = get(enzyme))) +
  geom_smooth(method=lm, color = brewer.pal(3,"Dark2")[3], se=F) +
  geom_point(aes(),size = 4,color = brewer.pal(3,"Dark2")[3],) +
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12,face = "bold"),
        axis.title.x  = element_text(size=12,face = "bold"),
        legend.text = element_text(size=10),
        legend.position = "right",
        legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
        panel.background = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.05),
        aspect.ratio = .9,
        panel.border = element_rect(colour = "black", fill=NA, size=.05),
        strip.text.x = element_text(size = 3)) +
  xlab("PSME4 Abundance") + ylab(paste0(enzyme,"-like cleavage"))
ggsave(paste0(path,enzyme,"correlation.pdf"),useDingbats=FALSE,height = 3)
c_term_agg_cleavage %>% select(PSME4,chymotryptic,tryptic,caspase) %>% cor(method = 'spearman')


corMAT = data.frame()
i = 1
for (AAX in AA_properties$Letter){
  corX = cor.test(c_term_agg_cleavage[,AAX],c_term_agg_cleavage$PSME4,method = "spearman")
  corMAT[i,"Rho"] = corX$estimate
  corMAT[i,"PVAL"]= corX$p.value
  corMAT[i,"AA"] = AAX
  i = i+1
}
corMAT = cor(c_term_agg_cleavage)

Cterm_aavPSME4 = c_term_agg_cleavage %>% select(PSME4,AA_properties$Letter) %>% pivot_longer(cols = AA_properties$Letter,values_to = "Percent of AA",names_to = "Amino.Acid") %>% data.frame()
colnames(Cterm_aavPSME4) = c("PSME4","Amino.Acid","Percent.Norm")
Cterm_aavPSME4 = merge(Cterm_aavPSME4,AA_properties,by.x = "Amino.Acid",by.y = "Letter")
Cterm_aavPSME4 %>%
  ggplot(aes(x = PSME4,y = Percent.Norm,color = Cleavage)) +
  facet_wrap(.~Amino.Acid, scales = "free") +
  geom_smooth(method=lm ,se=F) +
  geom_point(aes(),size = 4) +
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12,face = "bold"),
        axis.title.x  = element_text(size=12,face = "bold"),
        legend.text = element_text(size=10),
        legend.position = "right",
        legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
        panel.background = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.05),
        aspect.ratio = .9,
        panel.border = element_rect(colour = "black", fill=NA, size=.05),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12)) + scale_color_manual(values = c(brewer.pal(3,"Dark2")[c(1,3)],"grey",brewer.pal(3,"Dark2")[2])) +
  xlab("PSME4 Abundance") + stat_cor(color = "black",label.y.npc = .9,method = "spearman")
ggsave(paste0(path,enzyme,"corAAwPSME4.pdf"),useDingbats=FALSE,height = 10)
ggsave(paste0(path,enzyme,"corAAwPSME4.tiff"),useDingbats=FALSE,height = 10)


#Peptide Heatmap
# peptides = read.csv(paste0(path,"060219_1_Analysis\\060219_1_FilteredData\\060219_1_Lung-cancer-midgam_peptides.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
# Pep_sub = peptides[,c(128:144)]
# 
# Pep_sub = Pep_sub[,-9]
# Pep_sub = t(scale(t(Pep_sub)))
# colnames(Pep_sub) = str_sub(colnames(Pep_sub),end = -14)
# sample_annot = HeatmapAnnotation(df = data.frame(Sample = c(rep("Adjacent",8),rep("Tumor",8))),
#                                  col = list(Sample = c("Adjacent" = "skyblue3","Tumor" = "darkorange3")))
peptides = read.csv(paste0(path,"060219_1_Analysis\\060219_1_DataSnapshots\\060219_1_D_PeptidesAfterFilteringNonValidPeptides.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
Pep_sub = peptides[,c(111:126)]
Pep_sub = Pep_sub[rowSums(is.na(Pep_sub))<15,]
Pep_sub[is.na(Pep_sub)] = log10(((10^min(Pep_sub,na.rm = TRUE))/2))
Pep_sub = t(scale(t(Pep_sub)))
colnames(Pep_sub) = str_sub(colnames(Pep_sub),end = -9,start = 35)

sample_annot = HeatmapAnnotation(df = data.frame(Sample = rep(c("Tumor","Adjacent"),8)),
                                 col = list(Sample = c("Adjacent" = "skyblue3","Tumor" = "darkorange3")))
colnames(Pep_sub) = SampRename[colnames(Pep_sub) ,"NewName"]
Prot_diff_MAP = Heatmap(Pep_sub,
                        show_row_names = FALSE,
                        show_column_names = TRUE,
                        show_row_dend =  FALSE,
                        show_column_dend = TRUE,
                        cluster_columns = TRUE,
                        clustering_distance_columns = "spearman",
                        clustering_distance_rows = "spearman",
                        cluster_rows = TRUE,
                        name = "MAPP intensity",
                        top_annotation = sample_annot,
                        row_names_gp = gpar(fontsize = 8),
                        height = 10,
                        row_km = 3,
                        row_km_repeats = 100,
                        width = 5,
                        row_title = " ")
pdf(paste0(path,"heatmapPeptides.pdf")) 
draw(Prot_diff_MAP)
dev.off() 

#BarplotCleavages
c_term_percentchange = ((c_term_agg_norm[,9:16]) - (c_term_agg_norm[,1:8]))
c_term_percentchange.1 = sweep(c_term_percentchange,1,rowMeans(c_term_agg_norm[,1:16]),"/")
c_term_percentchange.1$Cleavage =  AA_properties[rownames(c_term_percentchange.1),"Cleavage"]
c_term_percentchange.1$AA = rownames(c_term_percentchange.1)
leter_order = c("R","K","H","P","S","G","A","T","L","V","Y","M","F","I","W","C","Q","N","D","E")  
c_term_percentchange.1$AA = factor(c_term_percentchange.1$AA,levels = leter_order)
c_term_percentchange.long = c_term_percentchange.1 %>% pivot_longer(cols = T.107:T.76,names_to = "Sample",values_to = "PercentChange") 

c_term_percentchange.long %>%
  ggplot(aes(y = PercentChange,x = AA,fill = Cleavage)) +
  geom_boxplot() + geom_hline(yintercept = 0) +
  theme(axis.text.y   = element_text(size=10,face = "bold"),
        axis.text.x   = element_text(size=10,face = "bold"),
        axis.title.y  = element_text(size=12,face = "bold"),
        axis.title.x  = element_text(size=12,face = "bold"),
        legend.text = element_text(size=11),
        legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = .5,
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_color_manual(values = cleavage_colors,aesthetics = c("fill","color")) 
ggsave(paste0(path,"AAchange_MAPP_bar_order.pdf"),useDingbats=FALSE,height = 4)

cterm = Heatmap(c_term_agg_norm1[leter_order,],
                show_row_names = TRUE,
                show_column_names = TRUE,
                show_row_dend =  FALSE,
                show_column_dend = FALSE,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                col = col_fun,
                rect_gp = gpar(col = "white", lwd = 0.5),
                height = unit(10,"cm"),
                name = "normalized percentage",
                column_title = "C Terminus Cleavage")

pdf(paste0(path,"cTermCleavage_cleavereorder_noannot.pdf")) 
draw(cterm)
dev.off() 

AA.count = colnames(peptides_sub)[10:29]
AA_counts_raw = peptides_sub %>% select(A.Count:V.Count)
intensity_norm[,AA.count] = (AA_counts_raw / rowSums(AA_counts_raw))*100
AA_comp_MAT = data.frame(matrix(nrow = 20,ncol = 16))
rownames(AA_comp_MAT) = AA.count
colnames(AA_comp_MAT) = cleanColname_pep
for (AA in AA.count){
  for (samp in cleanColname_pep){
    AA_comp_MAT[AA,samp] =   sum(intensity_norm[,AA]*intensity_norm[,samp]) / sum(intensity_norm[,samp])
  }
}

tumors = c(1,3,5,7,9,11,13,15)
normals= c(2,4,6,8,10,12,14,16)
AA_comp_MAT$T = rowMeans(AA_comp_MAT[,tumors])
AA_comp_MAT$N = rowMeans(AA_comp_MAT[,normals])
AA_comp_MAT$AA = str_split(rownames(AA_comp_MAT),"\\.",simplify = TRUE)[,1]
AA_comp_MAT %>%
  ggplot(aes(x = N,y = T)) +
  geom_point(aes(),size = 2) +
  geom_text_repel(aes(label = AA)) +
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12,face = "bold"),
        axis.title.x  = element_text(size=12,face = "bold"),
        legend.text = element_text(size=10),
        legend.position = "right",
        legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
        panel.background = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.05),
        aspect.ratio = .9,
        panel.border = element_rect(colour = "black", fill=NA, size=.05),
        strip.text.x = element_text(size = 3)) +
  xlim(0,10) + ylim(0,10)+ geom_abline()+
  xlab("Adjacent AA composition (%)") + ylab("Tumor AA composition (%)")
ggsave(paste0(path,"AAcomp_tumorvnormal.pdf"),useDingbats=FALSE,height = 3)




#agg by protein
AA_properties = read.csv(paste0(path,"AA_properties.csv"), header = 1,stringsAsFactors = FALSE) 
peptides = read.csv(paste0(path,"060219_1_Analysis\\060219_1_FilteredData\\060219_1_Lung-cancer-midgam_peptides.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
#peptides_sub = peptides[!str_detect(peptides$First.Gene.Name,"PSM") & !str_detect(peptides$First.Gene.Name,"IGH") & !str_detect(peptides$First.Gene.Name,"IGK"),]
peptides_sub = peptides
peptide_IDcol = colnames(peptides)[64:79]
cleanColname_pep = str_sub(peptide_IDcol,start = 36,-3)
names_framepep  <-  as.data.frame(cleanColname_pep)
names_framepep$peptide = peptide_IDcol
names_frame = merge(names_frame,names_framepep,by.x = "cleanColname",by.y = "cleanColname_pep")

intensity_norm = peptides_sub[,names_frame$peptide]
intensity_norm[is.na(intensity_norm)] = 0
intensity_norm[intensity_norm >=1] = 1
intensity_norm = intensity_norm[,names_frame$peptide]
colnames(intensity_norm) = unlist(names_frame$cleanColname)
sums = colSums(intensity_norm,na.rm = TRUE)
intensity_norm$N = as.factor(peptides_sub$Amino.acid.before)
intensity_norm$C = as.factor(peptides_sub$Last.amino.acid)
intensity_norm$Gene = as.factor(peptides_sub$First.Gene.Name)
c_term_agg = aggregate(. ~ C + Gene,intensity_norm[,-c(17)], FUN=sum)
T_cols = names_frame[names_frame$Sample == "Tumor","cleanColname"]
N_cols = names_frame[names_frame$Sample == "Adjacent","cleanColname"]
gene_sums_c = aggregate(. ~ Gene,intensity_norm[,-c(17,18)], FUN=sum)
rownames(gene_sums_c) = gene_sums_c$Gene
c_term_agg.norm = c_term_agg
c_term_agg.norm[,3:18] = (c_term_agg[,3:18] / gene_sums_c[c_term_agg$Gene,2:17])*100
pepcount = apply(gene_sums_c[,2:17],1,median)
Genes_covered = names(pepcount[pepcount>5])
c_term_agg.norm$T = apply(c_term_agg.norm[,T_cols],1,median,na.rm = TRUE)
c_term_agg.norm$N = apply(c_term_agg.norm[,N_cols],1,median,na.rm = TRUE)
c_term_agg.norm$total = apply(c_term_agg.norm[,c(N_cols,T_cols)],1,median,na.rm = TRUE)
c_term_agg.norm$percentchange = (c_term_agg.norm$T-c_term_agg.norm$N)/ c_term_agg.norm$total
c_term_table = c_term_agg.norm %>% select(Gene,percentchange,C) %>% pivot_wider(names_from = C,names_prefix = "count.",values_from = percentchange) %>% data.frame()
rownames(c_term_table) = c_term_table$Gene
c_term_table.ordered = c_term_table[names(pepcount[rev(order(pepcount))]),]
write.csv(c_term_table.ordered,file = paste0(path,"countsperprotein.csv"))


