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
library(patchwork)
library("ggExtra")
themeX = theme(axis.text.y   = element_text(size=12),
               axis.text.x   = element_text(size=12,hjust = .5),
               axis.title.y  = element_text(size=14,face = "bold"),
               axis.title.x  = element_text(size=14,face = "bold"),
               legend.text = element_text(size=10),
               legend.position = "right",
               legend.background = element_rect(fill = rgb(1, 1, 1, 0)),
               panel.background = element_blank(),
               plot.margin = margin(0, 0, 0, 0, "cm"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black",size = 0.05),
               aspect.ratio = 1,
               panel.border = element_rect(colour = "black", fill=NA, size=.05),
               strip.background = element_blank(),
               strip.text.x = element_text(size = 10))
path =  "D:\\Dropbox (Weizmann Institute)\\PSME4\\additionalData\\Immunopeptidomics_KD\\Ctrl4Exc\\HLApsme4KD_withoutCtrl4KD4\\"

# path =  "C:\\Users\\ajavi\\Dropbox (Weizmann Institute)\\PSME4\\additionalData\\Immunopeptidomics_KD\\Ctrl4Exc\\HLApsme4KD_withoutCtrl4KD4\\"
# Peptides = read.csv(paste0(path,"Analysis\\FilteredData\\HLApsme4KD_withoutCtrl4KD4_peptides.txt"),sep =  "\t",header = 1,stringsAsFactors = FALSE)


path =  "C:\\Users\\ajavi\\Dropbox (Weizmann Institute)\\PSME4\\additionalData\\Immunopeptidomics_KD\\Ctrl4Exc\\HLApsme4KDwithoutCtrl4KD3\\"
Peptides = read.csv(paste0(path,"Analysis\\FilteredData\\HLApsme4KDwithoutCtrl4KD3_peptides.txt"),sep =  "\t",header = 1,stringsAsFactors = FALSE)

path.out = (paste0(path,"output_revision_imprand\\"))
dir.create(path.out)
AAproperties = read.csv(paste0(path,"AA_properties.csv"))
# Peptides= read.csv(paste0(path,"HLApsme4KD\\Analysis\\DataSnapshots\\D_PeptidesAfterFilteringNonValidPeptides.txt"),sep =  "\t",header = 1,stringsAsFactors = FALSE)
# abun_col = colnames(Peptides)[c(81:83,85:87)]
# 
# Peptides[,abun_col][Peptides[,abun_col] == "NA"] = NA
# halfmin = log10((10^min(Peptides[,abun_col],na.rm = TRUE)) / 2)
# Peptides[,abun_col][is.na(Peptides[,abun_col])] = halfmin
# for (i in (1:nrow(Peptides))){
#   Peptides[i,"PVAL"] = t.test(Peptides[i,c(81:83)],Peptides[i,c(85:87)])$p.value
#   Peptides[i,"Ctrl.MED"] = median(unlist(Peptides[i,c(81:83)]))
# 
#   Peptides[i,"KD.MED"] = median(unlist(Peptides[i,c(85:87)]))
# 
# 
# }
abun_col = colnames(Peptides)[c(74:76,77:79)]

for (i in (1:nrow(Peptides))){
  Peptides[i,"PVAL"] = t.test(Peptides[i,c(74:76)],Peptides[i,c(77:79)])$p.value
  Peptides[i,"Ctrl.MED"] = median(unlist(Peptides[i,c(74:76)]))

  Peptides[i,"KD.MED"] = median(unlist(Peptides[i,c(77:79)]))


}
Peptides$FC = log2((10^Peptides$KD.MED)/(10^Peptides$Ctrl.MED))
PEP.PCA = prcomp(t(Peptides[,abun_col]),scale = TRUE,center = TRUE)
PEP.PCA.rot = data.frame(PEP.PCA$x)
PEP.PCA.rot$Treatment = str_split(rownames(PEP.PCA.rot),"_",simplify = TRUE)[,3]
rownames(PEP.PCA.rot) = paste0(str_split(rownames(PEP.PCA.rot),"_",simplify = TRUE)[,3],".",str_split(rownames(PEP.PCA.rot),"_",simplify = TRUE)[,4])
PEP.PCA.rot$Sample = rownames(PEP.PCA.rot)
PEP.PCA.rot %>%
  ggplot(aes(x = PC1,y = PC2,color = Treatment)) +
  geom_point(size = 4) + geom_text_repel(aes(label = Sample)) + themeX +
  scale_color_manual(values = c("darkgrey","darkblue"))
ggsave(filename = paste0(path.out,"PCA_noNorm.tiff"),height =3)
ggsave(filename = paste0(path.out,"PCA_noNorm.pdf"),height =3,useDingbats=FALSE)

#Corr MAT
CORMAT = cor(log10(Peptides[,abun_col]),method = "spearman")

rownames(CORMAT) = colnames(CORMAT) = paste0(str_split(rownames(CORMAT),"_",simplify = TRUE)[,3],".",str_split(rownames(CORMAT),"_",simplify = TRUE)[,4])
pdf(paste0(path.out,"Corplot.pdf"),height = 4)
corrplot(CORMAT,type = 'lower',tl.col = 'black',addCoef.col = 'black')
dev.off()

#Haplotype Comp
haplotype = read.csv(paste0(path,"peptide_haplotype_index.csv"),header = 1,stringsAsFactors = FALSE)
Peptide_hap = merge(Peptides,haplotype,by.x = "Sequence",by.y = "Peptide")

Peptide_hap %>% filter(min.rank <=5) %>%
ggplot(aes(x = FC,y = -log10(PVAL))) +
  geom_point() + geom_hline(yintercept = 1.3) + geom_vline(xintercept = c(-.7,.7))+
 themeX +  xlab("Log2 (PSME4 KD / WT)") + ylab("-Log10(P.value")
ggsave(filename = paste0(path.out,"Volcano_rank5.tiff"),height = 4)
ggsave(filename = paste0(path.out,"Volcano_rank5.pdf"),height = 4)

Peptide_hap %>% filter(min.rank <=2) %>%
  ggplot(aes(x = FC,y = -log10(PVAL))) +
  geom_point() + geom_hline(yintercept = 1.3) + geom_vline(xintercept = c(-.7,.7))+
  themeX +  xlab("Log2 (PSME4 KD / WT)") + ylab("-Log10(P.value")
ggsave(filename = paste0(path.out,"Volcano_rank2.tiff"),height = 4)
ggsave(filename = paste0(path.out,"Volcano_rank2.pdf"),height = 4)

Peptide_hap.long = Peptide_hap %>% select(Ctrl.MED,KD.MED,Sequence) %>% 
  pivot_longer(cols = c(Ctrl.MED,KD.MED), names_to = "Treat",values_to = "Abundance")

Peptide_hap.long %>%
  ggplot(aes(x = Treat,y = Abundance,fill = Treat)) +
  geom_violin() + 
  geom_boxplot(aes(), fill = "white", notch = FALSE,width=0.1,size = 1,
               coef = 0,outlier.size = 0,outlier.stroke = 0,show.legend = FALSE) +
  themeX + scale_fill_manual(values = c("darkgrey","darkblue")) + stat_compare_means()
ggsave(filename = paste0(path.out,"Abun_comp.tiff"),height = 4)
ggsave(filename = paste0(path.out,"Abun_comp.pdf"),height = 4)


Peptide_hap.long %>%
  ggplot(aes(x = Treat,y = Abundance,color = Treat)) +
  geom_beeswarm() + 
  geom_boxplot(aes(), fill = "white", notch = FALSE,width=0.1,size = 1,
               coef = 0,outlier.size = 0,outlier.stroke = 0,show.legend = FALSE) +
  themeX + scale_color_manual(values = c("darkgrey","darkblue")) + stat_compare_means(label = "p.signif")
ggsave(filename = paste0(path.out,"Abun_comp_dots.tiff"),height = 4)
ggsave(filename = paste0(path.out,"Abun_comp_dots.pdf"),height = 4)

#Protein Number
Count =  read.csv(paste0(path,"Analysis\\Statistics_Summary\\C_PeptidesDatasetGeneralStatistics.txt"),sep = "\t", header = 1,stringsAsFactors = FALSE)
Count$Treatment = str_split(Count$X,"_",simplify = TRUE)[,c(3)]
Count$Treatment = factor(Count$Treatment)
Count$Rep = str_split(Count$X,"_",simplify = TRUE)[,c(4)]
Count$Name = paste0(Count$Treatment,".",Count$Rep)
Count %>%
  ggplot(aes(x = Treatment,y = ValidValues)) +
  geom_dotplot(aes(fill = Treatment),binaxis='y', stackdir='center', stackratio=1.5, dotsize=1,lwd = .1) + 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.4,lwd = .1) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom="errorbar", width=0.2,lwd = .1) +
  themeX +
  scale_fill_manual(values = c("darkgrey","darkblue")) + ylim(0,NA)+
  xlab("") + ylab("Peptide Count") + stat_compare_means()
ggsave(file = paste0(path.out,"_peptideCount.pdf"),useDingbats=FALSE,height = 3)
ggsave(file = paste0(path.out,"_peptideCount.tiff"),height = 3)


Peptide_hap.rank = Peptide_hap[order(Peptide_hap$FC),]
Peptide_hap.rank$rank = 1:nrow(Peptide_hap.rank)
Peptide_hap.rank %>%
  ggplot(aes(x = rank,y = FC,fill = -log10(PVAL))) +
  geom_bar(stat = "identity") + themeX +theme(aspect.ratio = .5) +
  xlab("Peptides") +ylab("Log2 KD/WT") 
ggsave(file = paste0(path.out,"barRanked.pdf"),useDingbats=FALSE,height = 2)
ggsave(file = paste0(path.out,"barRanked.tiff"),height = 2)

Peptide_hap %>% filter(min.rank <=5)  %>% group_by(haplotype) %>% summarise(mean = mean(FC))

Peptides_1 = Peptide_hap %>% filter(PVAL <=0.1) %>% filter(min.rank <=5) 

write.csv(Peptides_1,file = paste0(path.out,"Peptide_hap5.csv"))

Peptides_2 = Peptide_hap %>% filter(PVAL <=0.1) %>% filter(min.rank <=2) 

write.csv(Peptides_2,file = paste0(path.out,"Peptide_hap2.csv"))

Peptides_3 = Peptide_hap %>% filter(PVAL <=0.1) %>% filter(min.rank <=2) %>% filter(abs(FC) >.7) 

write.csv(Peptides_3,file = paste0(path.out,"Peptide_hap2FC.csv"))

Peptides_4 = Peptide_hap %>% filter(PVAL <=0.1) %>% filter(min.rank <=5) %>% filter(abs(FC) >.7) 

write.csv(Peptides_4,file = paste0(path.out,"Peptide_hap5FC.csv"))


Peptides_5 = Peptide_hap  %>% filter(min.rank <=5) 

write.csv(Peptides_5,file = paste0(path.out,"Peptide_NoCutoffs.csv"))

# 
# Pep.Test = Peptides_1
# Set.Name  = "Bind5noFC"

# Pep.Test = Peptides_2
# Set.Name  = "Bind2noFC"

# Pep.Test = Peptides_3
# Set.Name  = "Bind2FC50"

# Pep.Test = Peptides_4
# Set.Name  = "Bind5FC50"
# 
Pep.Test = Peptides_5
Set.Name  = "NoCutoff"

Pep.Test$posend = as.factor(str_sub(Pep.Test$Sequence,-1,-1))
KD = Pep.Test %>%filter(FC > 0) %>% group_by(posend)  %>% count() %>% data.frame %>% tibble::column_to_rownames("posend")
WT= Pep.Test %>%filter(FC < 0) %>% group_by(posend)  %>% count() %>% data.frame %>% tibble::column_to_rownames("posend")


KD$percent = (KD$n / sum(KD$n))*100
KD = KD[rev(order(KD$percent)),]
KD[1,"ind"] = 1
KD[1,"cumulative"] = KD[1,"percent"]
for (i in 2:nrow(KD)){
  KD[i,"cumulative"] = KD[i-1,"cumulative"] + KD[i,"percent"]
  KD[i,"ind"] = i
  
}
KD$cond = "KD"
AA_total = KD

WT$percent = (WT$n / sum(WT$n))*100
WT = WT[rev(order(WT$percent)),]
WT[1,"ind"] = 1
WT[1,"cumulative"] = WT[1,"percent"]
for (i in 2:nrow(WT)){
  WT[i,"cumulative"] = WT[i-1,"cumulative"] + WT[i,"percent"]
  WT[i,"ind"] = i
  
}
WT$cond = "WT"
AA_total = rbind(AA_total,WT)

End1 = AA_total %>%
  ggplot(aes(x = ind,y = cumulative,color = cond)) +
  geom_line(size = 1) + geom_point(size = 3) +
  themeX + theme(aspect.ratio = .7) +
  scale_y_continuous(breaks=seq(0,100,20)) + geom_hline(yintercept = 0) +
  scale_color_manual(values = c("darkred","darkgrey")) +xlab("# of Residues") + ylab("cumulative percentage") +ggtitle("End position")

KD_stack = KD %>% arrange(-percent) %>%
  ggplot(aes(y = percent,fill = percent,x = "end")) +
  geom_bar(stat = "identity") + 
  geom_hline(yintercept = 0) + 
  themeX  + theme(aspect.ratio = 3, legend.position = "none") +ggtitle("PSME4 KD") + scale_fill_gradient(low = "white",high = "darkred") + 
  scale_y_continuous(breaks=seq(0,100,20))

WT_stack = WT %>% arrange(-percent) %>%
ggplot(aes(y = percent,fill = percent,x = "end")) +
  geom_bar(stat = "identity") + 
  geom_hline(yintercept = 0) +
  themeX  + theme(aspect.ratio = 3, legend.position = "none") +ggtitle("PSME4 WT") + scale_fill_gradient(low = "white",high = "black") + 
  scale_y_continuous(breaks=seq(0,100,20))
End1+ ggtitle(Set.Name) +KD_stack + WT_stack 
ggsave(filename = paste0(path.out,Set.Name,"_cumulative.tiff"),height = 7)
ggsave(filename = paste0(path.out,Set.Name,"_cumulative.pdf"),height = 7)



#Pie Chart
AAorder = AAproperties$Letter
KD$AA.fact = factor(rownames(KD),levels = AAorder)
cols = (inlmisc::GetColors( 20, scheme = "discrete rainbow",reverse = TRUE))
names(cols) = AAorder
data =  KD
data <- data %>% filter(percent != 0) %>%
  arrange(desc(AA.fact)) %>%
  mutate(prop = percent  / sum(data$percent ) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
KD.PIE = data %>%
  ggplot(aes(x="", y=prop, fill=AA.fact)) +
  geom_bar(stat="identity", width=1, color="white",alpha = .5) +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, x = 1.4, label = AA.fact), color = "black", size=6) +
  scale_fill_manual(values = cols) + ggtitle("KD")

WT$AA.fact = factor(rownames(WT),levels = AAorder)
data =  WT
data <- data %>% filter(percent != 0) %>%
  arrange(desc(AA.fact)) %>%
  mutate(prop = percent  / sum(data$percent ) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
WT.PIE = data %>%
  ggplot(aes(x="", y=prop, fill=AA.fact)) +
  geom_bar(stat="identity", width=1, color="white",alpha = .5) +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, x = 1.4, label = AA.fact), color = "black", size=6) +
  scale_fill_manual(values = cols) + ggtitle("WT")
KD.PIE + ggtitle(paste0("KD         ",Set.Name)) + WT.PIE

ggsave(filename = paste0(path.out,Set.Name,"_PIE.tiff"),height = 4)
ggsave(filename = paste0(path.out,Set.Name,"_PIE.pdf"),height = 4)


Pep.Test %>%
  ggplot(aes(x = haplotype,y = FC,fill = haplotype)) +
  geom_violin() + geom_hline(aes(yintercept = mean(FC)),linetype = "dashed") +
  geom_boxplot(aes(), fill = "white", notch = FALSE,width=0.1,size = 1,
               coef = 0,outlier.size = 0,outlier.stroke = 0,show.legend = FALSE) +
  themeX + scale_fill_manual(values = c("#81B6B6","#2F7979","#909FC3","#3D5387","#9ADB9A")) +
  ylab("Log2 (PSME4 KD / WT)") + stat_compare_means(ref.group = ".all.",label = "p.signif")
ggsave(filename = paste0(path.out,"haplotypePref.tiff"),height = 4)
ggsave(filename = paste0(path.out,"haplotypePref.pdf"),height = 4)

# scatter
Combo = merge(WT,KD,by = "row.names",suffixes = c(".WT",".KD"))
Combo = merge(AAproperties,Combo,by.x = "Letter",by.y = "Row.names")

Combo %>% 
  ggplot(aes(x = percent.WT,y = percent.KD,color = Biochem)) +
  geom_point() + geom_text_repel(aes(label = Letter)) +
  geom_abline() + themeX
ggsave(filename = paste0(path.out,"AAscatter.tiff"),height = 4)
ggsave(filename = paste0(path.out,"AAscatter.pdf"),height = 4)

Pep.Test %>%
  ggplot(aes(x = posend,y = FC,fill = posend)) +
  geom_violin() + geom_hline(yintercept = 0) +
  geom_boxplot(aes(), fill = "white", notch = FALSE,width=0.1,size = 1,
               coef = 0,outlier.size = 0,outlier.stroke = 0,show.legend = FALSE) +
  themeX + theme(aspect.ratio = 3) + coord_flip() +
  ylab("Log2 (PSME4 KD / WT)")
ggsave(filename = paste0(path.out,"AA-end.tiff"),height = 4)
ggsave(filename = paste0(path.out,"AA-end.pdf"),height = 4)
AA.sigframe = data.frame()
for (i in (1:nrow(AAproperties))){
  AAX = AAproperties[i,"Letter"]
  AA.x = Pep.Test %>% filter(posend == AAX)
  AA.sigframe[i,"AA"] = AAX
  AA.sigframe[i,"count"] = nrow(AA.x)
  if (nrow(AA.x) >10){
  AA.sigframe[i,"med"] = mean(AA.x$FC)
  AA.sigframe[i,"sig"] = t.test(AA.x$FC,Pep.Test$FC)$p.value
  }
}
AA.sigframe$neglog10p = -log10(AA.sigframe$sig)

AA.sigframe %>% 
  ggplot(aes(x = med,y = neglog10p)) +
  geom_point() + geom_text_repel(aes(label = AA)) +
    geom_hline(yintercept = 1.3) + geom_vline(xintercept = c(0))+
  themeX +  xlab("Log2 (PSME4 KD / WT)") + ylab("-Log10(P.value)")
ggsave(filename = paste0(path.out,"cleavageVolcano.tiff"),height = 4)
ggsave(filename = paste0(path.out,"cleavageVolcano.pdf"),height = 4)



Pep.Test.ID = merge(Pep.Test,AAproperties,by.x = "posend",by.y = "Letter")

Collapsebycleavage  = Pep.Test.ID %>%
  select(Intensity.4911_HLApsme4KD_CTRL_1_log10:Intensity.4918_HLApsme4KD_KD_4_log10,
         posend,Proteasome.Cleavage.Acitivty) %>%
  pivot_longer(cols = Intensity.4911_HLApsme4KD_CTRL_1_log10:Intensity.4918_HLApsme4KD_KD_4_log10,
               names_to = "treatment",values_to = "Abun") %>%
  group_by(treatment,Proteasome.Cleavage.Acitivty) %>%
  summarise(mAbun = mean((Abun)))
Collapsebycleavage$Group = str_split(Collapsebycleavage$treatment,"_",simplify = TRUE)[,3]

Collapsebycleavage %>%
  ggplot(aes(x= Group,y = mAbun))+
  geom_dotplot(aes(fill = Group),binaxis='y', stackdir='center', stackratio=1.5, dotsize=1,lwd = .1) + 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.4,lwd = .1) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom="errorbar", width=0.2,lwd = .1) +
  themeX + facet_wrap(~Proteasome.Cleavage.Acitivty) +
  stat_compare_means()

Pep.Test.ID.long  = Pep.Test.ID %>%
  select(Intensity.4911_HLApsme4KD_CTRL_1_log10:Intensity.4918_HLApsme4KD_KD_4_log10,
         posend,Proteasome.Cleavage.Acitivty) %>%
  pivot_longer(cols = Intensity.4911_HLApsme4KD_CTRL_1_log10:Intensity.4918_HLApsme4KD_KD_4_log10,
               names_to = "treatment",values_to = "Abun")
Pep.Test.ID.long$Group = str_split(Pep.Test.ID.long$treatment,"_",simplify = TRUE)[,3]

Pep.Test.ID.long %>%
  ggplot(aes(x= Group,y = Abun,fill = Group))+
  geom_violin() + 
  geom_boxplot(aes(), fill = "white", notch = FALSE,width=0.1,size = 1,
               coef = 0,outlier.size = 0,outlier.stroke = 0,show.legend = FALSE) +
  themeX + facet_wrap(~Proteasome.Cleavage.Acitivty) +
  stat_compare_means() +
  scale_fill_manual(values = c("darkgrey","darkblue")) 
ggsave(filename = paste0(path.out,"cleavageViolin.tiff"),height =5)
ggsave(filename = paste0(path.out,"cleavageViolin.pdf"),height = 5)


Pep.Test.ID.median.long = Pep.Test.ID %>%
  select(Ctrl.MED,KD.MED,Proteasome.Cleavage.Acitivty) %>%
  pivot_longer(cols = c(Ctrl.MED,KD.MED),names_to = "Group",values_to = "Abun")
Pep.Test.ID.median.long %>%
  ggplot(aes(x= Group,y = Abun,fill = Group))+
  geom_violin() + 
  geom_boxplot(aes(), fill = "white", notch = FALSE,width=0.1,size = 1,
               coef = 0,outlier.size = 0,outlier.stroke = 0,show.legend = FALSE) +
  themeX + facet_wrap(~Proteasome.Cleavage.Acitivty) +
  stat_compare_means() +
  scale_fill_manual(values = c("darkgrey","darkblue")) 
ggsave(filename = paste0(path.out,"cleavageViolin_medians.tiff"),height = 5)
ggsave(filename = paste0(path.out,"cleavageViolin_medians.pdf"),height = 5)


AA.sigframe.ID = merge(AA.sigframe,AAproperties,by.x = "AA",by.y = "Letter")
color_Activity = c(brewer.pal(3,"Dark2")[c(3,1,2)],"grey")
names(color_Activity) = c("Chymotryptic","Caspase","Tryptic","other")
AA.sigframe.ID %>% 
  ggplot(aes(x = med,y = neglog10p,color = Proteasome.Cleavage.Acitivty )) +
  geom_point() + geom_text_repel(aes(label = AA),size = 10) +
  geom_hline(yintercept = 1.3) + geom_vline(xintercept = c(0))+
  themeX +  xlab("Log2 (PSME4 KD / WT)") + ylab("-Log10(P.value)") +
  scale_color_manual(values = color_Activity)
ggsave(filename = paste0(path.out,"cleavageVolcano.tiff"),height = 4)
ggsave(filename = paste0(path.out,"cleavageVolcano.pdf"),height = 4)

AA.sigframe.ID.forBAR = AA.sigframe.ID %>% filter(!is.na(med))
AA.sigframe.ID.forBAR$Score = sign(AA.sigframe.ID.forBAR$med) * AA.sigframe.ID.forBAR$neglog10p 
AA.sigframe.ID.forBAR = AA.sigframe.ID.forBAR[order(AA.sigframe.ID.forBAR$Score),]
AA.sigframe.ID.forBAR$AA = factor(AA.sigframe.ID.forBAR$AA,levels = AA.sigframe.ID.forBAR$AA)
AA.sigframe.ID.forBAR %>%
  ggplot(aes(y = AA,x = Score,color = Proteasome.Cleavage.Acitivty)) +
  geom_point(aes(size = neglog10p))+
  themeX + theme(aspect.ratio = 2) + geom_vline(xintercept = 0,linetype = "dashed") +
   scale_color_manual(values = color_Activity) +
  ylab("C Terminus") + xlab("PSME4 KD / WT")
ggsave(filename = paste0(path.out,"cleavagebar.tiff"),height = 3)
ggsave(filename = paste0(path.out,"cleavagebar.pdf"),height = 3)

Volcano = Pep.Test.ID %>% filter(min.rank <=5) %>%
  ggplot(aes(x = FC,y = -log10(PVAL))) +
  geom_point(aes(color = Proteasome.Cleavage.Acitivty),pch = 19,size = 1) + geom_hline(yintercept = 1.3) + geom_vline(xintercept = c(-.7,.7))+
  themeX +  xlab("Log2 (PSME4 KD / WT)") + ylab("-Log10(P.value)") +   scale_color_manual(values = color_Activity)


pdf(paste0(path.out,"Volcano_rank2_cleavagecolor.pdf"),height = 4,useDingbats = FALSE)
ggMarginal(Volcano, groupColour = TRUE, groupFill = TRUE)
dev.off()
