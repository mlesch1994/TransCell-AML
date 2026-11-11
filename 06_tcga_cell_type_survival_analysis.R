library(sigatlas)
library(ggplot2)
library(survival)
library(survminer)

# Paths are relative to the project folder: inputs in data/, figures in figures/.
DATA_DIR <- "data"
PLOT_DIR <- "figures/tcga_survival_top50"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

load(file.path(DATA_DIR, "tcga_panc_symbol.RData"))
pheno = read.csv(file.path(DATA_DIR, "phenoDF_withage.csv"))
survival = read.csv(file.path(DATA_DIR, "tcga_pan-cancer_Survival.txt"), header = T, sep = "\t") 

AML_meta = pheno[pheno$sample.type =="primary" & pheno$cancer =="Acute Myeloid Leukemia", ]
AML_meta$sample.type = gsub("-", ".", AML_meta$sample.id)
AML_exp = exp_matrix_symbol[, colnames(exp_matrix_symbol) %in% AML_meta$sample.type]

sample_files = list.files(file.path(DATA_DIR, "markers"), full.names = TRUE)
genelist=list()
features=list()

for (sam in sample_files){
  nom = strsplit(sam, '/')[[1]]
  nom = nom[length(nom)]
  nom = strsplit(nom, '.csv')[[1]]
  cell_type = read.csv(sam, row.names = 1)
  cell_type = cell_type[c(1:50),]
  genelist[[nom]]=cell_type$gene
  features <- append(features, nom)
}
  
gsea = sigatlas(expr = AML_exp, cellmarker = genelist, tissue = "Other", organism = "Other")
gsea = t(gsea)

survival = survival[,c(1:3, 26:29)]
AML_survival = survival[survival$cancer.type.abbreviation == "LAML",]
AML_survival$sample.id = gsub("-", ".", AML_survival$sample)
names(AML_survival)[c(4,5)] = c("status","survival")

####survial analysis
gsea_survival = merge(gsea, AML_survival, by.x=0, by.y="sample.id")
features = unlist(features)

for (feature in features){
  gsea_survival$feature = gsea_survival[, feature]
  
  c <- ggplot(gsea_survival, aes(feature, survival))
  c + stat_smooth(method = "lm") + geom_point() + xlab("expression")
  
  gsea_survival$expressed = NA
  cutoffs = quantile(gsea_survival$feature)
  gsea_survival$expressed = ifelse(gsea_survival$feature > cutoffs[3], "yes", "no")
  
  survdiff(Surv(survival, status) ~ expressed, data=gsea_survival)
  
  my.fit1 <- survfit(Surv(survival, status) ~ expressed, data = gsea_survival) 
  p <- ggsurvplot(my.fit1, xlab = "days", legend.title = feature, pval = T, risk.table = F,
                                   ggtheme = theme_survminer(font.tickslab = c(16, "bold", "black")))
  out_file <- file.path(PLOT_DIR, paste0(feature, "_survival.pdf"))
  ggsave(out_file, plot = p$plot, width = 5.5, height = 4)
  
  ### Hazards ratio calculation
  cox<- coxph(Surv(survival, status) ~ expressed, data = gsea_survival)
                   
  print(feature)
  print(summary(cox))
}

