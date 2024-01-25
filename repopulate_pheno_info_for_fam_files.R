#!/usr/bin/env Rscript
# /users/k1642468/commonmind/Rodrigo/TWAS_HERVs_SCZ/repopulate_pheno_info_for_fam_files.R
# This script is gonna extract pheno info from a table and populate unknown Sex/Pheno fields in fam file.
library(data.table)
library(dplyr)
args = commandArgs(TRUE)
# args=list()  # ${CHIP}_genotypes ${CHIP}_genotypes_2 ${CHIP}
# args[1]="h650_genotypes"
# args[2]="h650_genotypes_2"
## args[3]="h650" # not needed anymore
typeof(args)
str(args)
print(args[1])
print(args[2])
print(args[3])

fam <- fread(paste0(args[1],".fam"), stringsAsFactors=F, h=F)
colnames(fam) <- c("FID","IID","Father_ID","Mother_ID","Sex", "Phenotype") # Sex 1 M / 2 F. Pheno 1 CTRL, 2 CASE, -9 (NA) or continuous (for non binary traits)
pd <- fread( "/users/k1642468/commonmind/Rodrigo/TWAS_HERVs_SCZ/pd_910_all", stringsAsFactors=F, h=T)
pd_2 <- dplyr::select(pd, genotype_id,Profile, Gender)
# replace cases with 2, ctrls with 1, females with 2, males with 1 
pd_2 <- pd_2 %>% mutate(Profile=replace(Profile, Profile=="SCZ", 2)) %>% as.data.frame()
pd_2 <- pd_2 %>% mutate(Profile=replace(Profile, Profile=="BD", 2)) %>% as.data.frame()
pd_2 <- pd_2 %>% mutate(Profile=replace(Profile, Profile=="AFF", 2)) %>% as.data.frame()
pd_2 <- pd_2 %>% mutate(Profile=replace(Profile, Profile=="CTRL", 1)) %>% as.data.frame()
pd_2 <- pd_2 %>% mutate(Gender=replace(Gender, Gender=="M", 1)) %>% as.data.frame()
pd_2 <- pd_2 %>% mutate(Gender=replace(Gender, Gender=="F", 2)) %>% as.data.frame()
pd_2$Profile <- as.numeric(as.factor(pd_2$Profile))
pd_2$Gender <- as.numeric(as.factor(pd_2$Gender))

# merge files and filter
head(pd_2)

#if(args[3]=="CMC1") { # if clause removed because sample names in plink files were normalized across chips (for compatibility with script 1, genotype processing prior to imputation)

fam_2 <- fam
fam_2$IID <- sub("^0_", "", fam_2$IID)
order_of_fam <- fam_2$IID
fam_3 <- merge( fam_2, pd_2, by.x="IID", by.y="genotype_id")
fam_3 <- select(fam_3,  FID, IID,Father_ID, Mother_ID, Gender, Profile)
fam_3 <- data.frame(fam_3)
row.names(fam_3)<- fam_3$IID
fam_3 <- fam_3[order_of_fam,]
print(head(fam_3))
print(all.equal(fam_3$IID, fam_2$IID))
write.table(fam_3, paste0(args[2],".fam"), sep="\t", quote=F, col.names=F, row.names=F)
system(paste0("cp ",args[1],".bim ", args[2],".bim"))
system(paste0("cp ",args[1],".bed ", args[2],".bed"))
print("End of script")
# } else { 

# order_of_fam<- fam$IID
# fam_2 <- merge( fam, pd_2, by.x="IID", by.y="genotype_id")
# fam_2 <- select(fam_2,  FID, IID,Father_ID, Mother_ID, Gender, Profile)
# fam_2 <- data.frame(fam_2)
# row.names(fam_2)<- fam_2$IID
# fam_2 <- fam_2[order_of_fam,]
# print("order of fam files match?")
# print(all.equal(fam$IID, fam_2$IID))
# write.table(fam_2, paste0(args[2],".fam"), sep="\t", quote=F, col.names=F, row.names=F)
# system(paste0("cp ",args[1],".bim ", args[2],".bim"))
# system(paste0("cp ",args[1],".bed ", args[2],".bed"))

# } 
