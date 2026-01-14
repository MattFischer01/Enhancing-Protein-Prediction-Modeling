#PC-AiR for European Population in MESA Proteomics WGS Data, can be adapted for other populations as needed
#Date: November 19, 2024

#I will be using PLINK files created post LD Pruning and post removing individuals from the outgroup in the first PCAIR
#PC-AiR - Euro Population
#2963 individuals left after filtering outgroup from 1st pcair 

#Write out new plink files with individuals from the main group, using LD pruned plink files - ran in shell
#plink -bfile ~/2024-09-10_TOPMed_MESA_WGS_QC/QC_MESAproteomics_allchr/europop/QC_geno0.02_mind0.02/PCA/03europop_Ldpruned --keep-fam ~/2024-09-10_TOPMed_MESA_WGS_QC/QC_MESAproteomics_allchr/allpops_merged_geno0.02_mind0.02/PCAIR/ingroup_IDs.txt --make-bed --out ~/2024-09-10_TOPMed_MESA_WGS_QC/QC_MESAproteomics_allchr/europop/QC_geno0.02_mind0.02/PCAIR/europop_filtered

library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(dplyr)
library(tibble)

#Convert PLINK files to GDS
snpgdsBED2GDS(bed.fn = "~/2024-09-10_TOPMed_MESA_WGS_QC/QC_MESAproteomics_allchr/europop/QC_geno0.02_mind0.02/PCAIR/europop_filtered.bed", 
              bim.fn = "~/2024-09-10_TOPMed_MESA_WGS_QC/QC_MESAproteomics_allchr/europop/QC_geno0.02_mind0.02/PCAIR/europop_filtered.bim", 
              fam.fn = "~/2024-09-10_TOPMed_MESA_WGS_QC/QC_MESAproteomics_allchr/europop/QC_geno0.02_mind0.02/PCAIR/europop_filtered.fam", 
              out.gdsfn = "~/2024-09-10_TOPMed_MESA_WGS_QC/QC_MESAproteomics_allchr/europop/QC_geno0.02_mind0.02/PCAIR/genotype.gds")

#Close all previous GDS files
showfile.gds(closeall=TRUE)

#Open GDS file
gdsfile <- "~/2024-09-10_TOPMed_MESA_WGS_QC/QC_MESAproteomics_allchr/europop/QC_geno0.02_mind0.02/PCAIR/genotype.gds"
gds <- snpgdsOpen(gdsfile)

#Calculate the kinship coefficients 
king <- snpgdsIBDKING(gds)

kingMat <- king$kinship 
colnames(kingMat)<-king$sample.id
row.names(kingMat)<-king$sample.id

#Save king RDS file for future use (if needed)
saveRDS(king, file = "~/2024-09-10_TOPMed_MESA_WGS_QC/QC_MESAproteomics_allchr/europop/QC_geno0.02_mind0.02/PCAIR/King_matrix.RDS")

#Close all previous GDS files
showfile.gds(closeall=TRUE)

#Create genotype object needed for pcair function
geno <- GdsGenotypeReader(filename = "~/2024-09-10_TOPMed_MESA_WGS_QC/QC_MESAproteomics_allchr/europop/QC_geno0.02_mind0.02/PCAIR/genotype.gds") 
genoData <- GenotypeData(geno)

#Run PCAIR 
mypcair <- pcair(genoData, kinobj = kingMat, divobj = kingMat)

# plot top 2 PCs
plot(mypcair)
# plot PCs 2 and 3
plot(mypcair, vx = 2, vy = 3)
# plot PCs 3 and 4
plot(mypcair, vx = 3, vy = 4)

#Convert to PCA vec and val files
eigenvec<-mypcair$vectors %>% as.data.frame() %>% rownames_to_column(var="sample_id")
str(eigenvec)
val<-mypcair$values %>% as.data.frame()
fwrite(eigenvec,"~/2024-09-10_TOPMed_MESA_WGS_QC/QC_MESAproteomics_allchr/europop/QC_geno0.02_mind0.02/PCAIR/PCAIR.eigenvec",col.names = T,row.names = F,sep='\t')
fwrite(val,"~/2024-09-10_TOPMed_MESA_WGS_QC/QC_MESAproteomics_allchr/europop/QC_geno0.02_mind0.02/PCAIR/PCAIR.eigenval")