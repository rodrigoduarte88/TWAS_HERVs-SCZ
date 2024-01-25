#!/usr/bin/env Rscript
# /users/k1642468/commonmind/Rodrigo/TWAS_HERVs_SCZ/annotate_rsids_in_bim_file.R 
# This script is gonna extract genomic coordinates to create a temporary SNP name (chr_bp), and it will match this information to a non-redudant database of genomic coordinates containing rsids, based on snp151.
library(data.table)
library(dplyr)
args = commandArgs(TRUE)
# args <- list()
# args[1] <- "h650_genotypes_4" 
# args[2] <- "h650_genotypes_5"
typeof(args)
str(args)
print(args[1])
print(args[2])

bim <- fread(paste0(args[1],".bim"), stringsAsFactors=F, h=F)
print("Loading SNP reference, ~/commonmind/Rodrigo/TWAS_HERVs_SCZ/snp151_common/snp151Common_processed.txt.gz...")
dbsnp <- fread("~/commonmind/Rodrigo/TWAS_HERVs_SCZ/snp151_common/snp151Common_processed.txt.gz", h=T, stringsAsFactors=F)
colnames(bim) <- c("CHR","SNP","DUMMY", "BP", "A1", "A2")
bim$CHR_BP <- paste0(bim$CHR, "_", bim$BP) # to reformat column as CHR_BP. Output of imputation is hg19. dbsnp file is just chr_bp...
print("Removing variants with the same name or genomic coordinates..")
bim <- bim[!(duplicated(bim$SNP) | duplicated(bim$SNP, fromLast=TRUE) ) ,]
bim <- bim[!(duplicated(bim$CHR_BP) | duplicated(bim$CHR_BP, fromLast=TRUE) ) ,]


new_bim <- merge(bim, dbsnp, by.x="CHR_BP", by.y="CHRBP_hg19")

new_bim$allele_match <- ifelse((new_bim$A1==new_bim$A1_dbsnp & new_bim$A2==new_bim$A2_dbsnp), "match",   	
									ifelse((new_bim$A1==new_bim$A2_dbsnp & new_bim$A2==new_bim$A1_dbsnp), "match",   	
										ifelse((new_bim$A1==new_bim$A1_dbsnp_rc & new_bim$A2==new_bim$A2_dbsnp_rc), "match", 
											ifelse((new_bim$A1==new_bim$A2_dbsnp_rc & new_bim$A2==new_bim$A1_dbsnp_rc), "match","err")
											)) )					
table(new_bim$allele_match) 
# new_bim[new_bim$allele_match=="err",] # to see variants that are not matching...
# remove those that had mismatches
new_bim <- new_bim[!new_bim$allele_match=="err",]
print(paste0("there were ",nrow(bim)," variants before. Now there are ", nrow(new_bim),". This is a reduction of ", (nrow(bim)-nrow(new_bim)),", or ",round((nrow(bim)-nrow(new_bim))/nrow(bim)*100, digits=2),"%"))
snp_list <- data.frame(new_bim$SNP)
write.table(snp_list, paste0("snp_list_",args[1]), quote=F, col.names=F, row.names=F)
print("Writing plink files with variants that can be annotated to rs ids..")

system(paste0("plink --bfile ", args[1]," --extract ", paste0("snp_list_",args[1])," --make-bed --out ",args[1],"_temp"))
print("\n\nRelabelling variants in new set of plink files..")
filtered_vars <- fread(paste0(args[1],"_temp.bim"), stringsAsFactors=F, h=F)
colnames(filtered_vars) <- c("CHR","SNP","DUMMY", "BP", "A1", "A2")
filtered_vars$order <- row.names(filtered_vars)
filtered_vars$order <- as.numeric(filtered_vars$order)
filtered_vars_final <- merge( filtered_vars, new_bim, by.x="SNP", by.y="SNP")
filtered_vars_final <- filtered_vars_final[order(order),]
filtered_vars_final <- dplyr::select(filtered_vars_final, CHR.x, RSID_dbsnp, DUMMY.x, BP.x, A1.x, A2.x)
write.table(filtered_vars_final, paste0(args[1],"_temp.bim"), quote=F, sep="\t", col.names=F, row.names=F)

print("\n\nEnsuring only ACTG biallelic variants remain..")
system(paste0("plink --bfile ",args[1],"_temp --snps-only just-acgt --biallelic-only strict --make-bed --out ",args[2]))

print(paste0("The prefix of the output plink files is: ", args[2]))
##########################################################
