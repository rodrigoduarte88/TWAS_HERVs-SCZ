#!/usr/bin/Rscript
# Source: https://github.com/opain/Calculating-FUSION-TWAS-weights-pipeline/
options(verbose=TRUE)
library("optparse")

option_list = list(
  make_option("--PLINK_prefix", action="store", default=NA, type='character',
              help="Prefix of PLINK files for the target sample [required]"),
  make_option("--PLINK_prefix_chr", action="store", default=NA, type='character',
              help="Prefix of per chromosome PLINK files for the target sample [required]"),
  make_option("--phenotype_file", action="store", default=NA, type='character',
              help="File name for normalised and adjusted expression data [required]"),
  make_option("--coordinate_file", action="store", default=NA, type='character',
              help="File name for coordinates of genes/transcripts [required]"),
  make_option("--gene_name", action="store", default=NA, type='character',
              help="Name of gene or transcript to be processed [required]"),
  make_option("--plink", action="store", default=NA, type='character',
              help="Path to PLINK [required]"),
  make_option("--gcta", action="store", default=NA, type='character',
              help="Path to gcta_nr_robust binary [required]"),
  make_option("--gemma", action="store", default=NA, type='character',
              help="Path to gemma binary [required]"),
  make_option("--ld_ref_dir", action="store", default=NA, type='character',
              help="FUSION LD reference directory [required]"),
  make_option("--fusion_software", action="store", default=NA, type='character',
              help="FUSION software directory [required]"),
  make_option("--output_dir", action="store", default=NA, type='character',
              help="Directory name for the output [required]")
)

opt = parse_args(OptionParser(option_list=option_list))
system(paste('mkdir -p ',opt$output_dir,'/Output',sep=''))

root_dir<-getwd()

library(data.table)

# Read in the gene expression data to extract the gene's chromosome and boundary coordinates
Gene_coordinates_file<-fread(opt$coordinate_file, header=T)

CHR<-Gene_coordinates_file$X.Chr[Gene_coordinates_file$ID == opt$gene_name]
Start<-Gene_coordinates_file$start[Gene_coordinates_file$ID == opt$gene_name]-0.5e6
Stop<-Gene_coordinates_file$end[Gene_coordinates_file$ID == opt$gene_name]+0.5e6

############### <RD>
CHR<-CHR[1]
Start<-Start[1]
Stop<-Stop[1]

Start <- ifelse(Start<0, 0, Start)
# if (Start < 0) {
# Start <- 0
# }
################ </RD>

# Extract gene from phenotype file
system(paste('mkdir -p ',opt$output_dir,'/temp',sep=''))  # <RD: added -p to suprress warnings>
system(paste('awk -f ./t.awk c1=FID c2=IID c3=',opt$gene_name,' ',opt$phenotype_file,' > ',opt$output_dir,'/temp/temp_',opt$gene_name,'.pheno',sep=''))

# Using PLINK, extract variants within the specified gene +/- 500kb from start and stop coordinates
if(!is.na(opt$PLINK_prefix)){
	err_1<-system(paste(opt$plink,' --bfile ',opt$PLINK_prefix,' --make-bed --pheno ',opt$output_dir,'/temp/temp_',opt$gene_name,'.pheno', ' --out ', opt$output_dir,'/temp/temp_',opt$gene_name,' --geno 0.02 --chr ',CHR,' --from-bp ',Start,' --to-bp ',Stop,' --extract ',opt$ld_ref_dir,'/1000G.EUR.',CHR,'.bim', sep=''))
} else {
	err_1<-system(paste(opt$plink,' --bfile ',opt$PLINK_prefix_chr,CHR,' --make-bed --pheno ',opt$output_dir,'/temp/temp_',opt$gene_name,'.pheno', ' --out ', opt$output_dir,'/temp/temp_',opt$gene_name,' --geno 0.02 --chr ',CHR,' --from-bp ',Start,' --to-bp ',Stop,' --extract ',opt$ld_ref_dir,'/1000G.EUR.',CHR,'.bim', sep=''))
}

if (err_1 == 13) {
write.table('No SNPs within gene +/- 0.5Mb', paste(opt$output_dir,'/Output/',opt$gene_name,'.err',sep=''), col.names=F, row.names=F, quote=F)
} else {

# Using FUSION, calculate the weights for the genes expression using subset of genotypic data.
setwd(paste(opt$output_dir,'/temp', sep=''))
# The following line is required for gemma. Check FUSION twas manual.  HOWEVER, it stops the script! so cant use it.....?
system(paste('ln -s ./ output', sep=''))
system(paste('Rscript ',opt$fusion_software,'/FUSION.compute_weights.R --bfile ', opt$output_dir,'/temp/temp_',opt$gene_name,' --tmp temp_',opt$gene_name,'.tmp --out ', opt$output_dir,'/Output/',opt$gene_name,' --verbose 2 --save_hsq --PATH_gcta ',opt$gcta,' --PATH_gemma ',opt$gemma,' --PATH_plink ',opt$plink,' --models top1,lasso,enet,blup', sep=''))
}
# You can change which models you want to include, from top1,blup,bslmm,lasso,enet
# Delete the temporary files
system(paste('rm ',opt$output_dir,'/temp/temp_',opt$gene_name,'*', sep=''))

