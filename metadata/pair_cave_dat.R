#pair_cave_dat.R

#upload cave data
cdat = read_csv('metadata/cave_data.csv') %>% 
  mutate(length_m = length_ft*0.3048) %>% 
  select(-length_ft) %>% 
  data.frame()
rownames(cdat)=cdat$cave

#format into pairs
pcdat = data.frame()
cnames = c('cave1', 'cave2', 'ws1', 'ws2', 'dist', 'l1', 'l2')
caves = sort(unique(cdat$cave))
for (i in 1:(length(caves)-1) ){
  c1 = caves[i]
  others = caves[(i+1):length(caves)]
  for (c2 in others){
    w1=cdat[c1,'watershed']
    w2=cdat[c2,'watershed']
    dist=cdat[c1,paste('dist', c2, sep='_')]
    l1=cdat[c1,'length_m']
    l2=cdat[c2,'length_m']
    res.row = data.frame(c1,c2,w1,w2,dist,l1,l2)
    colnames(res.row)=cnames
    pcdat = rbind(pcdat, res.row)
  }
}

save(cdat,pcdat, file='metadata/cave_data.Rdata')
