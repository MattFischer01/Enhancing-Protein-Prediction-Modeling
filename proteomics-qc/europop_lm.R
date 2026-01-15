#This code is to adjust protein expression based on covariates: sidno, age_at_collection, gender, race, genetic PC1-10

#This code reads in the population protein expression data, selects for European population, performs linear model adjustment for covariates, and outputs adjusted protein expression data and averages.

#This code can be adapted for other populations by changing the race1c filter.

library(Matrix)
library(data.table)
library(dplyr)

#Read in all population protein expression data 
proteinexp <- fread("~/2024-11-12_TOPMed_MESA_Proteomics_QC/allcovariate_withproteinexpresion.csv", header= T, sep = ",")

#Select for European population
proteinexp <- filter(proteinexp, race1c == "EUR")

#save the averages for the original protein expression 
average_ogprot <- proteinexp %>% select(sidno, starts_with("OID"))
average_ogprot <- average_ogprot %>% group_by(sidno) %>% summarise(across(starts_with("OID"), mean, na.rm = T))

#write out the average of the original protein expression data
fwrite(average_ogprot, "~/2024-11-12_TOPMed_MESA_Proteomics_QC/europop/data/europop_average_og_prot.csv", col.names = T, sep = ",")

#Read in all genetic pcs 
genetic_pcs <- fread("~/2024-09-10_TOPMed_MESA_WGS_QC/QC_MESAproteomics_allchr/europop/QC_geno0.02_mind0.02/PCAIR/PCAIR.eigenvec", header = T, sep = "\t")

#Keep FID and first 10 PCs of genetic_pcs and add colnames of FID and first 10 PCs to genetic_pcs
genetic_pcs <- genetic_pcs[,1:11]
colnames(genetic_pcs) <- c("FID", "PC1", "PC2", "PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

#Combine proteinexp demographic and ID information with genetic PCS
covariates <- left_join(proteinexp[,1:12], genetic_pcs, by = "FID")

#Make sure the variables in covariates are the correct data types for lm. arrange by sidno
covariates$gender1 <- as.factor(covariates$gender1)
covariates$race1c <- as.factor(covariates$race1c)
covariates <- covariates %>% arrange(sidno)

#Create protein expression matrix with only protein expression data, arrange by sidno:
protein_expression_matrix <- proteinexp %>% select(sidno, starts_with("OID"))
protein_expression_matrix <- protein_expression_matrix %>% arrange(sidno)
protein_expression_matrix <- protein_expression_matrix %>% select(-sidno)
#Make sure its a matrix
protein_expression_matrix <- as.matrix(protein_expression_matrix)

#lm - protein expression adjustments based on covariates
residmat <- matrix(NA,nrow=dim(protein_expression_matrix)[1],ncol=dim(protein_expression_matrix)[2])

#Testing if ranef is adjusting for covariates. Loop with covariates.
for(i in 1:dim(protein_expression_matrix)[2]){
  fit <- lm(protein_expression_matrix[,i] ~ AGE_AT_COLLECTION + gender1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covariates)
  fitresid <- resid(fit)
  residmat[,i] <- fitresid
}

#Save the model 
save(fit, file = "~/2024-11-12_TOPMed_MESA_Proteomics_QC/europop/lm/lm_fit_model.RData")

#Add in sidno and protein colnames to adjusted matrix, then average 
colnames(residmat) <- colnames(protein_expression_matrix)

#write out the adjusted protein expression data with covariates
residmat_covariates <- cbind(covariates, residmat)
fwrite(residmat_covariates, "~/2024-11-12_TOPMed_MESA_Proteomics_QC/europop/lm/europop_adjusted_prot_withcovariates.csv", col.names = T, sep = ",")

#average the adjusted protein expression data by sidno 
residmat_sidno <- cbind(covariates$sidno, residmat)
colnames(residmat_sidno)[1] <- "sidno"
residmat_sidno <- as.data.frame(residmat_sidno)
average_adjustprot <- residmat_sidno %>% group_by(sidno) %>% summarise(across(starts_with("OID"), mean, na.rm = T))

#write out average adjusted protein expression data with sidnos in first column
fwrite(average_adjustprot, "~/2024-11-12_TOPMed_MESA_Proteomics_QC/europop/data/europop_average_adjusted_prot.csv", col.names = T, sep = ",")

covariates_distinct <- covariates %>% distinct(sidno, .keep_all = TRUE)
average_adjusted_covariates <- cbind(covariates_distinct, average_adjustprot)
average_adjusted_covariates <- average_adjusted_covariates[,-23]

#write out average adjusted protein expression data with all covariates 
fwrite(average_adjusted_covariates, "~/2024-11-12_TOPMed_MESA_Proteomics_QC/europop/lm/europop_average_adjusted_prot_withcovariates.csv", col.names = T, sep = ",")

#read in ranefmat from lmer, and only keep sidnos and protein expression data, then write out 
ranefmat_cov <- read_csv("2024-11-12_TOPMed_MESA_Proteomics_QC/europop/lmer/ranefmat_cov.csv")
ranefmat_sidno <- ranefmat_cov %>% select(sidno, starts_with("OID"))
fwrite(ranefmat_sidno, "~/2024-11-12_TOPMed_MESA_Proteomics_QC/europop/data/europop_lmer_ranefmat.csv", col.names = T, sep = ",")