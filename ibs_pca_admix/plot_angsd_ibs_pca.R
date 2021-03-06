#plot_angsd_ibs_pca.R
library(ggdendro)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
rm(list=ls())

# LOAD RESULTS FROM angsd_ibs_pca.R ---------------------------------------

#nesticus
ll=load('ibs_pca_admix/nesticus/pdat.Rdata')
ll
n.pdat = pdat
n.hc = hc
n.site=site

#ptomaphagus
ll=load('ibs_pca_admix/ptomaphagus/pdat.Rdata')
ll
p.pdat = pdat
p.hc = hc
p.site = site


# PLOT HIERARCHICAL CLUSTERING -------------------------------------------
watershedList = list('BT'='Guntersville Lake',
                     'GV'='Guntersville Lake',
                     'SB'='Upper Elk River',
                     'ST'='Upper Elk River')
PCA.SHAPE.CHOICE = c(25, 21)
TREE.SHAPE.CHOICE = c(21, 25)

my_plot_hclust = function(hc, groupLabel, watershedList,tree.point.size=2){
  dd.row=ggdendrogram(hc, rotate = FALSE, size = 2)
  dd.row=as.dendrogram(hc)
  ddata_x <- dendro_data(dd.row)
  ys=append(ddata_x$segments$y, ddata_x$segments$yend)
  yshift = -max(ys)/20
  p2 <- ggplot(segment(ddata_x)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank())
  labs <- label(ddata_x) %>% 
    mutate(site=substr(label, start=3, stop=4),
           label2=substr(label, start=5, stop=8),
           watershed = unlist(watershedList[site]))
  p2 + 
    # geom_text(data=labs, aes(label=label2, x=x, y=yshift, colour=labs$site), angle=0, hjust = 0.5) +
    geom_point(data=labs,aes(x=x,y=0, fill=labs$site, shape=watershed), size=tree.point.size,color='black') +
    scale_shape_manual(values=TREE.SHAPE.CHOICE) +
    labs(color='',
         y='height') +
    theme(legend.position = 'none')
}

NYLIM = c(-0.07, 0.36)
PYLIM = c(-0.07, 0.20)
REL_HEIGHTS = c(1,8)

#nesticus
n.tree0 = my_plot_hclust(n.hc, n.site, watershedList)
n.tree1 = n.tree0 +
  lims(y=NYLIM)
n.main = ggdraw() + draw_label('N. barri', fontface='italic')
n.tree= plot_grid(n.main, n.tree1, nrow=2, rel_heights = REL_HEIGHTS)
n.tree

#ptomaphagus
p.tree0 = my_plot_hclust(p.hc, p.site, watershedList)
p.tree1 = p.tree0 +
  lims(y=PYLIM)
p.main = ggdraw() + draw_label('P. hatchi', fontface='italic')
p.tree= plot_grid(p.main, p.tree1, nrow=2, rel_heights = REL_HEIGHTS)
p.tree



# MULTIDIMENTIONAL SCALING -----------------------------------------------

get_mod_lims=function(vector){
  lims = c(min(vector), max(vector))
  range=max(vector)-min(vector)
  buffer=range/10
  mod.lims = c( (min(vector)-buffer), (max(vector)+buffer) )
}

plot_mds = function(pdat, point.size=4,
                    letters,
                    shapeChoice = c(21,23),
                    invert1 = 1){
  pdat = pdat %>% 
    mutate(watershed = unlist(watershedList[Cave]),
           Dim1=Dim1*invert1,
           watershed = factor(watershed, levels=c('Upper Elk River', 'Guntersville Lake')))
  mod.lims1 = get_mod_lims(pdat$Dim1)
  mod.lims2 = get_mod_lims(pdat$Dim2)
  mod.lims3 = get_mod_lims(pdat$Dim3)
  
      
  mds12 = pdat %>% 
    ggplot(aes(x=Dim1, y=Dim2, fill=Cave, shape=watershed)) +
    geom_point(size=point.size, color='black') + 
    scale_shape_manual(values=shapeChoice) +
    labs(x='MDS1', y='MDS2') +
    theme(legend.position='none',
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    lims(x=mod.lims1,
         y=mod.lims2)
  
  mds13 = pdat %>% 
    ggplot(aes(x=Dim1, y=Dim3, fill=Cave, shape=watershed)) +
    geom_point(size=point.size, color='black') + 
    scale_shape_manual(values=shapeChoice) +
    labs(x='MDS1', y='MDS3', color='cave') +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position='bottom',
          legend.justification = 'center') +
    lims(x=mod.lims1,
         y=mod.lims3)
  plt=plot_grid(mds12, mds13+theme(legend.position = 'none'), labels=letters)
  
  #set up the legends
  fillPlt = pdat %>% 
    ggplot(aes(x=Dim1, y=Dim3, fill=Cave)) +
    geom_point(size=point.size, color='black', pch=21) +
    labs(fill='cave') +
    theme(legend.position='bottom',
          legend.justification = 'center')
  fillLegend = cowplot::get_legend(fillPlt)
  
  shapePlt = pdat %>% 
    ggplot(aes(x=Dim1, y=Dim3, shape=watershed)) +
    geom_point(size=4, color='black') +
    scale_shape_manual(values=PCA.SHAPE.CHOICE) +
    theme(legend.position='bottom',
          legend.justification = 'center')
  shapeLegend = cowplot::get_legend(shapePlt)
  l = plot_grid(fillLegend,
                shapeLegend)

  return(list('legend'=l,
              'plot'=plt))
}

n.mds.list = plot_mds(n.pdat, letters=c('b', 'c'), shapeChoice = PCA.SHAPE.CHOICE)
p.mds.list = plot_mds(p.pdat, letters=c('e', 'f'), shapeChoice = PCA.SHAPE.CHOICE,  invert1 = -1)
n.mds=n.mds.list[['plot']]
p.mds = p.mds.list[['plot']]


# ASSEMBLE ----------------------------------------------------------------

n.pan = plot_grid(n.tree, n.mds, nrow=2, labels=c('a'))
p.pan = plot_grid(p.tree, p.mds, nrow=2, labels=c('d'))
plts=plot_grid(n.pan, p.pan, nrow=1)
plot_grid(plts, n.mds.list[['legend']], nrow=2, rel_heights=c(15,1))

