#!/usr/bin/env Rscript
#vcf_dxy.R
#calculates dxy for sites in a vcf between two populations
#intended for unphased, biallelic SNP data
#intended to work similarly to vcftools's --weir-fst-pop syntax
#only works for two populations
#based on Molecular Population Genetics (Hahn 2019) page 87-88
#Equation:
#dxy = Sum( Xi*Yj*Kij )


# more on equation --------------------------------------------------------

#for all alleles i in population X and all alleles j in population Y
#where:
#Xi is the frequency of allele i in population X
#Yj is the frequency of allele j in population Y
#Kij is 0 if i==j and 1 if i!=j

#Example:
#At a locus, the alternative allele frequency is .25 in pop X and 0.75 in pop Y
#the reference allele frequencies are 0.75 in pop X and 0.25 in pop Y
#dxy for this locus is 0.375:

#dxy = 
#   Xalt * Yalt * Kalt.alt +
#   Xalt * Yref * Kalt.ref +
#   Xref * Yalt * Kref.alt +
#   Xref * Yref * Kre.ref

#or 
#   0.25 * 0.75 * 0 +
#   0.25 * 0.25 * 1 +
#   0.75 * 0.75 * 1 +
#   0.75 * 0.25 * 0

#or
#   (Xalt * Yref) + (Xref * Yalt) #this is what gets used here



# PARSE ARGUMENTS ---------------------------------------------------------

suppressMessages(library(optparse))
option_list = list(
  
  make_option(c("--vcf"), type="character", default=NULL, 
              help="The vcf to use as input"),
  
  make_option(c("--popX"), type="character", default=NULL, 
              help="Single column table of individuals from VCF in population X"),
  
  make_option(c("--popY"), type="character", default=NULL, 
              help="Single column table of individuals from VCF in population Y"),

  make_option(c("--out"), type="character", default='out.dxy', 
              help="Output name")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
vcfInput = opt$vcf
inPopX = opt$popX
inPopY = opt$popY
outName = opt$out


# IMPORT AND FORMAT -------------------------------------------------------------

#pops
popx = read.table(inPopX, stringsAsFactors=F)$V1
popy = read.table(inPopY, stringsAsFactors=F)$V1


#upload and format vcf data
suppressMessages(library(tidyverse))
vdat = read_tsv(file=vcfInput, comment='##')
pdat = vdat[,1:9]
sdat = vdat[,10:ncol(vdat)]
geno = data.frame(lapply(sdat, function(x) substr(x, 1,3)))

#split into alternative allele counts
a1 = data.frame(lapply(sdat, function(x) substr(x, 1,1)), stringsAsFactors=F)
a2 = data.frame(lapply(sdat, function(x) substr(x, 3,3)), stringsAsFactors=F)
colnames(a1) = paste(colnames(a1), 'allele1', sep='_')
colnames(a2) = paste(colnames(a2), 'allele2', sep='_')
counts = cbind(a1,a2)


#split genotype data by pop
keep1 = append(paste(popx, 'allele1', sep='_'), paste(popx, 'allele2', sep='_'))
dat1 = counts %>% 
  select(keep1) %>% 
  as_tibble()

keep2 = append(paste(popy, 'allele1', sep='_'), paste(popy, 'allele2', sep='_'))
dat2 = counts %>% 
  select(keep2) %>% 
  as_tibble()


# GET DXY -----------------------------------------------------------------
#get allele frequencies
get_p = function(a){
  a[a=='.']<-NA
  x=na.omit(as.numeric(a))
  p = sum(x) / length(x)
  return(p)
}

#set up variables for dxy
Xalt = apply(dat1, 1, function(x) get_p(x))
Yalt = apply(dat2, 1, function(x) get_p(x))
Xref = 1-Xalt
Yref = 1-Yalt


#get dxy
dxy = Xalt*Yref + Yalt*Xref

#write out
colnames(pdat)[1]='CHROM'
pdat %>%
    select(CHROM, POS) %>%
    cbind(dxy) %>% 
    write_tsv(path=outName)
  




