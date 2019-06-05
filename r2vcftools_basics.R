
########## r2vcftools installation --------------------------------------------------
## Install VCFtools
# GitHub: https://github.com/vcftools/vcftools
# Manual: https://vcftools.github.io/man_latest.html

## Install the r2vcftools library from GitHub
install.packages("devtools")
devtools::install_github("bcm-uga/LEA")
devtools::install_github("nspope/r2vcftools", force=T)

### Load r2vcftools
library(r2vcftools)

########## Help functions ----------------------------------------------------------
?Chrom
?Copy
?Filter
?Fstats
?Geno
?GenotypeMatrix
?Linkage
?Lfmm
?LfmmEnv
?LoadMeta
?Query
?Relatedness
?Save
?Subset

## Load VCFsummary function
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}

########## Load VCF file -----------------------------------------------------------
## Download example VCF file "Imaurandioides.vcf" from figshare: https://ndownloader.figshare.com/files/10990757
snps_raw <-  vcfLink("Imaurandioides.vcf", overwriteID=T)
snps_raw

########## Filter loci (sites) that are not SNPs ---------------------------------------------
snps <- Filter(snps_raw, filterOptions(max.alleles=2, min.alleles=2), indels="remove")

########## Visualize basic information and diversity metrics -----------------------
## Basic info
VCFsummary(snps_raw) 
VCFsummary(snps)
snps@sample_id
snps@site_id
snps@meta

##Chromosome, possitions, and IDs (without a reference genome CHROM = Contig)
Chrom(snps) 

## Genotype Matrix
genotypes <- GenotypeMatrix(snps_unind) 
genotypes[1:10, 1:10] ## -1 is missing; otherwise, gives the number of derived alleles in the genotype -- ie. a 0 and 2 are both homozygotes

## Heterozygosity and inbreeding per individual
HE <- Query(snps_ind_site_miss, type="het")
hist(HE$O.HOM) ## O.HOM, E.HOM, N_SITES, F
hist(HE$E.HOM) 
hist(HE$N_SITES) 
hist(HE$F) ## Inbreeding coefficient

## Nucleotide diversity (PI) per site
PI <- Query(snps_ind_site_miss, type="site-pi")
mean(PI$PI) ## Mean nucleotide divergency per-site
hist(PI$PI)

########## Filter individuals ------------------------------------------------------
## Remove repeated individual
snps@meta
UNIND <- snps@sample_id[snps@sample_id != "IC213_rep_sorted.bam"]
snps_unind <- Subset(snps, samples=UNIND) ## "Subset" is used to keep individuals or loci
snps_unind@meta
VCFsummary(snps_unind)

# Missing data per individual
Missing_ind <- apply(GenotypeMatrix(snps_unind), 1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) 
hist(Missing_ind) # Define threshold for individual missing data tolerance

## Remove individuals with threshold-defined missing data
Missingind.df <- as.data.frame(Missing_ind)
Missingind.df$ID <- row.names(Missingind.df)
Missingind.df[Missingind.df$Missing_ind>30,] 
indtokeep <- Missingind.df[Missingind.df$Missing_ind <= 30,]
snps_ind_miss <- Subset(snps_unind, samples = indtokeep$ID)
snps_ind_miss@meta
VCFsummary(snps_unind)
VCFsummary(snps_ind_miss)

########## Filter loci (sites) by missing data ------------------------------------------------------------
## Missing data per locus (site)
Missing_lc <- apply(GenotypeMatrix(snps_ind_miss), 2, function(x) sum(x < 0)/length(x)*100)
summary(Missing_lc) 
hist(Missing_lc) # Define threshold for site missing data tolerance

## Remove loci (sites) with threshold-defined missing data
snps_ind_site_miss <- Filter(snps_ind_miss, filterOptions(max.missing = 0.7)) ## max.missing = proportion of individuals without missing data to be set as threshold
VCFsummary(snps_ind_miss)
VCFsummary(snps_ind_site_miss) 

########## Filter loci (sites) by quality and depth ---------------------------------------------
## Site quality
quality <- Query(snps_ind_site_miss, type="site-quality")
summary(quality$QUAL)
hist(quality$QUAL) 
hist(quality$QUAL[quality$QUAL<400]) ## Define filtering threshold

## Site depth 
site.depth <- Query(snps_ind_site_miss, type="site-mean-depth")
summary(site.depth$MEAN_DEPTH)
hist(site.depth$MEAN_DEPTH) 
hist(site.depth$MEAN_DEPTH[site.depth$MEAN_DEPTH <100]) ## Define filtering threshold

## Filter
snps_QD <- Filter(snps_ind_site_miss, filterOptions(minQ=30, min.meanDP=20, max.meanDP=800)) 
VCFsummary(snps_ind_site_miss) 
VCFsummary(snps_QD) 

########## Filter loci (sites) by MAF and HWE ---------------------------------------------------
## Allele frequencies per site
freq <- Query(snps_ind_site_miss, type="freq2")
nbelow_maf <- length(freq$X.FREQ.[freq$X.FREQ.<0.05])
hist(freq$X.FREQ., main=paste(nbelow_maf, "loci have a MAF < 0.05")); abline(v=0.05, col="red")

## HWE per site
HWE <- Query(snps_ind_site_miss, type="hardy")
summary(HWE$P_HWE)
hist(HWE$P_HWE)
nbelow_hwe <- length(HWE$P_HWE[HWE$P_HWE<0.0001])
hist(HWE$P_HWE, main=paste(nbelow_hwe, "loci have a HWE_P < 0.05")); abline(v=0.05, col="red")

## Filter
snps_MH <- Filter(snps_ind_site_miss, filterOptions(maf=0.05, hwe=0.0001)) 
VCFsummary(snps_QD) 
VCFsummary(snps_MH) 

########## Filter loci (sites) by Linkage Disequilibrium --------------------------------------------------
### LD within contigs
ld_within <- Linkage(snps_MH, type="geno-r2", linkageOptions(min.r2=0.0))
head(ld_within)
hist(ld_within$R.2) ## Define filtering threshold

ld_within <- Linkage(snps_MH, type="geno-r2", linkageOptions(min.r2=0.5)) ## Sites with high LD
hist(ld_within$R.2)

ld_snps <- ld_within$ID1
nold_snps <- snps_MH@site_id[!(snps_MH@site_id %in% ld_snps)] 
snps_WLD <- Subset(snps_MH, sites=nold_snps) # Keep snps that are not in LD.
VCFsummary(snps_MH) 
VCFsummary(snps_WLD) 

## LD between contigs
ld_between <- Linkage(snps_WLD, type="interchrom-geno-r2", linkageOptions(min.r2=0)) 
hist(ld_between$R.2) ## Define filtering threshold

ld_between <- Linkage(snps_WLD, type="interchrom-geno-r2", linkageOptions(min.r2=0.5)) ## Sites with high LD

ld2_snps <- ld_between$ID1
nold2_snps <- snps_WLD@site_id[!(snps_WLD@site_id %in% ld2_snps)]
snps_BLD <- Subset(snps_WLD, sites=nold2_snps) # Keep snps that are not in LD.
VCFsummary(snps_WLD) 
VCFsummary(snps_BLD)

########## Save filtered vcf -------------------------------------------------------
Save(snps_BLD, "Imaurandioides_filtered.vcf")

########## Genetic diversity metrics -------------------------------------------------------
## Load GenDiv function
GenDiv <- function(vcf){
  HE <- Query(vcf, type="het")
  
  HE$HO <- (HE$N_SITES-HE$O.HOM.)/(HE$N_SITES) ## Observed heterozygosity (HO)
  error_ho <- qt(0.975,df=length(HE$HO)-1)*sd(HE$HO)/sqrt(length(HE$HO))
  ho <- mean(HE$HO)
  left_ho <- ho-error_ho
  right_ho <- ho+error_ho
  HO.df <- data.frame(He_O=ho, low_He_O=left_ho, up_He_O=right_ho)
  
  HE$HE <- (HE$N_SITES-HE$E.HOM.)/(HE$N_SITES) ## Expected heterozygosity (HE)
  error_he <- qt(0.975,df=length(HE$HE)-1)*sd(HE$HE)/sqrt(length(HE$HE))
  he <- mean(HE$HE)
  left_he <- he-error_he
  right_he <- he+error_he
  HE.df <- data.frame(He_E=he, low_He_E=left_he, up_He_E=right_he)
  
  error_f <- qt(0.975,df=length(HE$F)-1)*sd(HE$F)/sqrt(length(HE$F))
  f <- mean(HE$F)
  left_f <- f-error_f
  right_f <- f+error_f
  F.df <- data.frame(F=f, low_F=left_f, up_F=right_f)
  
  PI <- Query(vcf, type="site-pi")
  error_pi <- qt(0.975,df=length(PI$PI)-1)*sd(PI$PI)/sqrt(length(PI$PI))
  pi <- mean(PI$PI)
  left_pi <- pi-error_pi
  right_pi <- pi+error_pi
  PI.df <- data.frame(PI=pi, low_PI=left_pi, up_PI=right_pi)
  
  print(paste0("Genetic diversity metrics")) 
  RES <- cbind(HO.df,  HE.df,  F.df,  PI.df)
  return(RES)
}

GenDiv(snps_BLD)

########## Genetic distances -------------------------------------------------------
## Create dummy population ID and insert into VCF meta data - option 1
snps_BLD@meta$Pop1 <- rep(1:3, each = 10, length.out=nrow(snps_BLD@meta))

## Create dummy population ID and attach it to VCF meta data - option 2
df <- data.frame(ID= snps_BLD@sample_id, Pop2=rep(1:3, each = 10, length.out=nrow(snps_BLD@meta)))
row.names(df) <- snps_BLD@sample_id
snps_BLD_meta <- LoadMeta(snps_BLD, df)
snps_BLD_meta@meta

## Calculate pairwise FST (between populations)
Fstats(snps_BLD, by=snps_BLD@meta$Pop1)

## Calculate pairwise relatedness (between individuals)
Relatedness(snps_BLD)

########## Convert to geno and lfmm formats ----------------------------------------
## Convet VCF file to genotype matrix format
Geno(snps_BLD, output.file = "snps_BLD.geno")

## Converts the VCF file to the format used by LFMM
Lfmm(snps_BLD, output.file = "snps_BLD.lfmm")







