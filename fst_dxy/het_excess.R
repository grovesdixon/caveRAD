#het_excess.R

CUT=0.2


# for nesticus ------------------------------------------------------------

#read in
hdat = read.table('mpileup_analyses/nesticus/out.hwe', header = T)
head(hdat)

#check density and CUT
hdat %>% 
  ggplot(aes(x=P_HET_EXCESS)) +
  geom_density() +
  geom_vline(xintercept = CUT, lty=2)

#filter
keep = hdat %>% 
  filter(P_HET_EXCESS > CUT)

#output
keep %>% 
  select(CHR, POS) %>% 
  write_tsv(path='mpileup_analyses/nesticus/sites_with_acceptable_heterozygosity.tsv', col_names = F)



# for ptomaphagus ---------------------------------------------------------


#read in
hdat = read.table('mpileup_analyses/ptomaphagus/out.hwe', header = T)
head(hdat)

#check density and CUT
hdat %>% 
  ggplot(aes(x=P_HET_EXCESS)) +
  geom_density() +
  geom_vline(xintercept = CUT, lty=2)

#filter
keep = hdat %>% 
  filter(P_HET_EXCESS > CUT)

#output
keep %>% 
  select(CHR, POS) %>% 
  write_tsv(path='mpileup_analyses/ptomaphagus/sites_with_acceptable_heterozygosity.tsv', col_names = F)
