#plot_admixture.R
library(tidyverse)
library(rlang)
library(cowplot)
theme_set(theme_cowplot())
rm(list=ls())
source('ibs_pca_admix/admixture_functions.R')



# GLOBAL VARS -------------------------------------------------------------

KPOPS=5
PREFIX='mydata_'
XLAB = 'cave'
YLAB = 'ancestry'
sppLabRelWidth=c(1,40)
mainColorList = list('BT' = '#7CAE00',
                  'GV' = '#00BFC4',
                  'SB' = '#F8766D',
                  'ST' = '#C77CFF')
mcols = unlist(mainColorList)


# BUILD FOR NESTICUS ------------------------------------------------------

#read in the data
dfList = read_admixture(samplesFile = 'ibs_pca_admix/nesticus/samples.txt',
                          qoptDir = 'ibs_pca_admix/nesticus/',
                          kpops = KPOPS,
                          prefix = PREFIX)

#apply the plotting function to the list of dataframes
n.splots=lapply(names(dfList), function(x) plot_admix_samples(x))
n.splots[[2]] = n.splots[[2]] + theme(axis.text.y = element_blank())
n.splots[[3]] = n.splots[[3]] + theme(axis.text.y = element_blank())
n.splots[[4]] = n.splots[[4]] + theme(axis.text.y = element_blank())
names(n.splots)=c('k2', 'k3', 'k4', 'k5')

#fix colors for k=4 to match rest of figures
#(I could not figure out a proper way to do this)
mcols
n.splots[[3]]
modCols = as.character(c(mcols[1],#bt
                         mcols[3],#gv
                         mcols[4],#st
                         mcols[2]))#sb
n.splots[[3]] + scale_fill_manual(values=modCols)
n.splots[[3]] = n.splots[[3]] + scale_fill_manual(values=modCols)

#plot the list of plots
nYLAB = 'N. barri'
n.plts=plot_grid(plotlist=n.splots, nrow=1)
n.ylab = ggdraw() + draw_label(nYLAB, angle=90, fontface='italic')
n.top0=plot_grid(n.ylab, n.plts, nrow=1, rel_widths=sppLabRelWidth)
# title=ggdraw() + draw_label('N. barri', fontface='italic')
# n.plt=plot_grid(title, top0, nrow=2, rel_heights=c(1,15))
# n.plt




# BUILD FOR PTOMAPHAGUS ---------------------------------------------------

#read in the data
dfList = read_admixture(samplesFile = 'ibs_pca_admix/ptomaphagus/samples.txt',
                        qoptDir = 'ibs_pca_admix/ptomaphagus/',
                        kpops = KPOPS,
                        prefix = PREFIX)

#apply the plotting function to the list of dataframes
p.splots=lapply(names(dfList), function(x) plot_admix_samples(x))
p.splots[[2]] = p.splots[[2]] + theme(axis.text.y = element_blank())
p.splots[[3]] = p.splots[[3]] + theme(axis.text.y = element_blank())
p.splots[[4]] = p.splots[[4]] + theme(axis.text.y = element_blank()) 

#fix colors for k=4 to match rest of figures
mcols
p.splots[[3]]
modCols = as.character(c(mcols[2],#sb
                         mcols[1],#gv
                         mcols[4],#st
                         mcols[3]))#bt
p.splots[[3]] + scale_fill_manual(values=modCols)
p.splots[[3]] = p.splots[[3]] + scale_fill_manual(values=modCols)


#plot the list of plots
pYLAB='P. hatchi'
p.plts=plot_grid(plotlist=p.splots, nrow=1)
p.ylab = ggdraw() + draw_label(pYLAB, angle=90, fontface='italic')
p.top0=plot_grid(p.ylab, p.plts, nrow=1, rel_widths=sppLabRelWidth)


#ASSEMBLE TOGETHER
tops = plot_grid(n.top0, p.top0, nrow=2)
mainYlab = ggdraw() + draw_label(YLAB, angle=90)
final = plot_grid(mainYlab, tops, nrow=1, rel_widths = c(1,40))
final


