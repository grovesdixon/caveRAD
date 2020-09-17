#angsd_ibs_pca.R
#This is the exploratory version of the script that looks at lots of plot types
#use plot_angsd_ibs_pca.R to make main figure

library(tidyverse)
library(vegan)
library(adegenet) 
library(cowplot)
library(ggdendro)
theme_set(theme_cowplot())



# SELECT SPECIES TO RUN ON ------------------------------------------------

#NESTICUS
rm(list=ls())
bamsFile="ibs_pca_admix/nesticus/bams"
ibsMatFile="ibs_pca_admix/nesticus/ibsResults.ibsMat"
pdatOut='ibs_pca_admix/nesticus/pdat.Rdata'
pcaCovFile="ibs_pca_admix/nesticus/PCAResults.covMat"

#or

#PTOMAPHAGUS
rm(list=ls())
bamsFile="ibs_pca_admix/ptomaphagus/bams"
ibsMatFile="ibs_pca_admix/ptomaphagus/ibsResults.ibsMat"
pdatOut='ibs_pca_admix/ptomaphagus/pdat.Rdata'
pcaCovFile="ibs_pca_admix/ptomaphagus/PCAResults.covMat"




# FIRST RUN FOR NESTICUS WITH EVERYTHING ----------------------------------

bams0=read.table(bamsFile)[,1] # list of bam files that matches the angsd runs
bams=sub('.bam', '', as.character(bams0))
goods=c(1:length(bams))
length(bams)

#--------------------
# build individual to population correspondences
site = substr(bams, start=3, stop=4) #'site' here refers to population
i2p = data.frame('sample'=bams, 'pop'=site)
rownames(i2p) = bams
head(i2p)

# settign up colors for plotting
palette(rainbow(length(unique(site))))
colors=as.numeric(as.factor(site))
colpops=as.numeric(as.factor(sort(unique(site))))

#-------------
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)
ma = as.matrix(read.table(ibsMatFile))
ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.7)
dd.row=ggdendrogram(hc, rotate = FALSE, size = 2)
dd.row=as.dendrogram(hc)
ddata_x <- dendro_data(dd.row) 
p2 <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())
labs <- label(ddata_x) %>% 
  mutate(site=substr(label, start=3, stop=4))
p2 + geom_text(data=labs,
               aes(label=label, x=x, y=0, colour=labs$site), angle=90, hjust = 1) +
  lims(y=c(-.1, 0.4)) +
  labs(color='')


# optimal number of clusters ----------------------------------------------
library(mclust)
d_clust2 <- Mclust(as.dist(ma), G=1:5)
d_clust2$BIC
plot(d_clust2$BIC)
m.best2 <- dim(d_clust2$z)[2]

cat("model-based optimal number of clusters:", m.best2, "\n")


# performing PCoA and CAP
conds=data.frame(cbind(site))
pp0=capscale(ma~1)
pp=capscale(ma~site,conds)

# significance of by-site divergence
adonis(ma~site,conds)

# eigenvectors
plot(pp0$CA$eig) 

axes2plot=c(1,2)  
cmd=pp0
plot(cmd,choices=axes2plot,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=transp(colors,alpha=0.7))
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cmd,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cmd,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops,label=T)


# unscaled, to identify outliers
plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=colpops,label=T)

#-------------------
#MDS AS SHOWN ON ANGSD WEBSITE (http://www.popgen.dk/angsd/index.php/PCA_MDS#Model)
#This indeed gives exact same thing as capscale
m <- as.matrix(read.table(ibsMatFile))
mds <- cmdscale(as.dist(m))
dim(mds)
rownames(mds)=bams
mds %>% 
  data.frame() %>% 
  mutate(sample=bams,
         cave=site) %>% 
  ggplot(aes(x=X1, y=X2, fill=site)) +
  geom_point(size=4, pch=21, color='black') + 
  labs(x='MDS1', y='MDS2') +
  theme(legend.position='none')

# PLOT FIGURE -------------------------------------------------------------


point.size=4
names(cmd)
pdat = data.frame(cmd$Ybar)
pdat$Cave=site
mds12 = pdat %>% 
  ggplot(aes(x=Dim1, y=Dim2, fill=Cave)) +
  geom_point(size=4, pch=21, color='black') + 
  labs(x='MDS1', y='MDS2') +
  theme(legend.position='none')

mds13 = pdat %>% 
  ggplot(aes(x=Dim1, y=Dim3, fill=Cave)) +
  geom_point(size=3, pch=21, color='black') + 
  labs(x='MDS1', y='MDS3')

l=cowplot::get_legend(mds13)


plot_grid(mds12, mds13+theme(legend.position = 'none'))
save(pdat, hc, site, file=pdatOut)

#--------------------
# covariance / PCA (not really needed, I prefer IBS)

# choose either of the following two covarince matrices:
co = as.matrix(read.table(pcaCovFile)) # covariance based on single-read resampling
#co = as.matrix(read.table("ok.covar")) # covariance by ngsCovar
co =co[goods,goods]
dimnames(co)=list(bams[goods],bams[goods])

# PCoA and CAP (constranied analysis of proximities)  
conds=data.frame(cbind(site))
pp0=capscale(as.dist(1-cov2cor(co))~1) # PCoA
pp=capscale(as.dist(1-cov2cor(co))~site,conds) # CAP

# significance of by-site divergence, based on 1-correlation as distance
adonis(as.dist(1-cov2cor(co))~site,conds)

# eigenvectors: how many are interesting?
plot(pp0$CA$eig) 

axes2plot=c(1,3)  
cc=pp0
plot(cc,choices=axes2plot,type="n") # choices - axes to display
points(cc,choices=axes2plot,pch=19,col=colors)
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cc,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cc,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops,label=T)


#-------------
# t-SNE:  machine learning to identify groups of samples 
# based on genotypes' correlations
# (only makes sense if you have hundreds of samples)

library(Rtsne)
library(vegan)
library(adegenet)
quartz()

# perplexity:  expected number fo neighbors. Set to 0.5x N(samples per pop)
perp=8
rt = Rtsne(as.dist(1-cov2cor(co)), perplexity=perp,max_iter=2,is_distance=T)
for (i in 1:100){
  rt = Rtsne(as.dist(1-cov2cor(co)), perplexity=perp,max_iter=10,Y_init=rt$Y,is_distance=T)
  plot(rt$Y,col=colors,pch=16,cex=0.8,main=i*10)
}
# ordispider(rt$Y,groups=site,col="grey80",alpha=0.01)
ordiellipse(rt$Y,groups= site,draw="polygon",col=colpops,label=T)







