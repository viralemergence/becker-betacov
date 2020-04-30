## phylofactor of betacov host status

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(ape)
library(phylofactor)
library(tidyverse)
library(data.table)
library(ggtree)

## load phylogeny
setwd("~/Desktop/batcov_gpf")
tree=readRDS('STFull.rds')

## load greg matching host taxa
setwd("~/Desktop/batcov_gpf")
data=readRDS('GregPredictions.rds')

## just bats
bats=data[data$hOrder=='Chiroptera',]

## get taxa info
btaxa=bats[c('hOrder','hFamily','hGenus')]
btaxa=btaxa[!duplicated(btaxa$hGenus),]

## just bat genera
gkeep=as.character(btaxa$hGenus)

## get tree tips
tdata=data.frame(tree$tip.label)
names(tdata)='tips'

## get genus
tdata$genus=sapply(strsplit(as.character(tdata$tips),'_'),function(x) x[1])

## trim to bats
tdata2=tdata[tdata$genus%in%gkeep,]

## trim
btree=keep.tip(tree,as.character(tdata2$tips))

## fix names
tdata=tdata2
rm(tree,tdata2,data,bats,gkeep)
tdata$unique_names=tdata$tips
tdata$tips=NULL

## fix taxa
names(tdata)=c('hGenus','unique_name')

## merge with taxa
tdata=merge(tdata,btaxa,by='hGenus',all.x=T)

## taxa
tdata$taxonomy=with(tdata,paste(hOrder,hFamily,hGenus,unique_name,sep='; '))

## load in CoV data
setwd("~/Desktop/batcov_gpf/Data")
library(tidyverse)
read_csv('BatCoV-assoc.csv') %>% filter(origin == 'Anthony') -> batcov
read_csv('BatCoV-assoc.csv') -> batcov

## just beta
batcov=batcov[which(batcov$virus_genus=='Betacoronavirus'),]

## make beta and subgenus
batcov %>% mutate(betacov = as.numeric(virus_genus == 'Betacoronavirus'),
                  sarbecov = as.numeric(virus_subgenus == 'Sarbecovirus')) -> batcov

## Replace NA's for subgenus with 0 only if they're not betacoronaviruses
batcov$sarbecov[is.na(batcov$sarbecov) & !(batcov$virus_genus == 'Betacoronavirus')] <- 0
batcov$sarbecov[is.na(batcov$sarbecov) & (batcov$virus_genus == 'Betacoronavirus')] <- -1

## collapse
shorthand <- function(x) {
  if(sum(x==1)>0) {1} else {
    if(sum(x==-1)>0) {-1} else {
      0
    }
  }
}
batcov %>% group_by(host_species) %>% 
  dplyr::summarize(betacov = max(betacov),
                   sarbecov = shorthand(sarbecov)) -> batcov

## fix host species
library(plyr)
batcov$host_species=revalue(batcov$host_species,
                            c('Dermanura phaeotis'='Artibeus phaeotis',
                              'Hipposideros commersonii'='Hipposideros commersoni'))

## unique_name
batcov$unique_name=gsub(" ","_",batcov$host_species)

## merge
batdf=merge(tdata,batcov,by='unique_name',all.x = T)
table(batdf$betacov)==nrow(batcov)

## replace na
batdf$betacov=replace_na(batdf$betacov,0)
batdf$sarbecov=na_if(batdf$sarbecov,-1)
batdf$sarbecov=replace_na(batdf$sarbecov,0)

## merge into phylogeny order
bdata=batdf[match(btree$tip.label,batdf$unique_name),]

## merge with caper
library(caper)
cdata=comparative.data(phy=btree,data=bdata,names.col=unique_name,vcv=T,na.omit=F,warn.dropped=T)

## add unique
cdata$data$unique_name=rownames(cdata$data)
cdata$data$Species=rownames(cdata$data)

## Holm rejection procedure
HolmProcedure <- function(pf,FWER=0.05){
  if (pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients['phyloS','Pr(>|t|)'])
  } else {
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients['phyloS','Pr(>|z|)'])
  }
  D <- length(pf$tree$tip.label)
  
  ## this is the line for Holm's sequentially rejective cutoff
  keepers <- pvals<=(FWER/(2*D-3 - 2*(0:(pf$nfactors-1))))
  
  
  if (!all(keepers)){
    nfactors <- min(which(!keepers))-1
  } else {
    nfactors <- pf$nfactors
  }
  return(nfactors)
}

## get species in a clade
cladeget=function(pf,factor){
  spp=pf$tree$tip.label[pf$groups[[factor]][[1]]]
  return(spp)
}

## set taxonomy
taxonomy=data.frame(cdata$data$taxonomy)
names(taxonomy)="taxonomy"
taxonomy$Species=rownames(cdata$data)
taxonomy=taxonomy[c("Species","taxonomy")]
taxonomy$taxonomy=as.character(taxonomy$taxonomy)

## gpf
set.seed(1)
pf=gpf(Data=cdata$data,tree=cdata$phy,
       frmla.phylo=betacov~phylo,
       family=binomial,algorithm='phylo',nfactors=5)
keep=HolmProcedure(pf,FWER=0.05)
pf.tree(pf,factors=1:keep,size=0.1,layout="circular")$ggplot

## refit with correct number
pf2=gpf(Data=cdata$data,tree=cdata$phy,
       frmla.phylo=betacov~phylo,
       family=binomial,algorithm='phylo',nfactors=keep)

## set key
setkey(pf2$Data,'Species')

## summarize
pf.taxa(pf2,taxonomy,factor=1)$group1

## get clade
cdata$data$pf_beta=ifelse(cdata$data$unique_name%in%cladeget(pf2,1),'clade','other')

## table
table(cdata$data$pf_beta,cdata$data$betacov)
prop.table(table(cdata$data$pf_beta,cdata$data$betacov))*100

## make data
set=cdata$data

## hist of pred
set$preds_allbetacov=predict(pf2,type='response')
hist(set$preds_allbetacov)

## refit to just the sabre
set.seed(1)
pf=gpf(Data=cdata$data,tree=cdata$phy,
       frmla.phylo=sarbecov~phylo,
       family=binomial,algorithm='phylo',nfactors=5)
keep=HolmProcedure(pf,FWER=0.05)
pf.tree(pf,factors=1:keep,size=0.1,layout="circular")$ggplot

## refit with correct number
pf2=gpf(Data=cdata$data,tree=cdata$phy,
        frmla.phylo=sarbecov~phylo,
        family=binomial,algorithm='phylo',nfactors=keep)

## set key
setkey(pf2$Data,'Species')

## summarize
pf.taxa(pf2,taxonomy,factor=1)$group1

## get clade
cdata$data$pf_sab=ifelse(cdata$data$unique_name%in%cladeget(pf2,1),'clade','other')

## table
table(cdata$data$pf_sab,cdata$data$sarbecov)
prop.table(table(cdata$data$pf_sab,cdata$data$sarbecov))*100

## predict
set$preds_sarbecovirus=predict(pf2,type='response')
hist(set$preds_sarbecovirus)

## trim set
set2=set[c('unique_name','betacov','sarbecov','preds_allbetacov','preds_sarbecovirus')]

## export
setwd("~/Desktop/batcov_gpf")
write.csv(set2,'bat betacov_phylofactor preds.csv')