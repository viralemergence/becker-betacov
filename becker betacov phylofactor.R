## phylofactor of betacoronavirus host status
## danbeck@iu.edu
## last updated 05/13/2020

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(ape)
library(data.table)
library(ggtree)
library(plyr)
library(phylofactor)
library(tidyr)

## load master data
setwd("~/Desktop/virionette/03_interaction_data")
data=read.csv('virionette.csv',header=T)
data$tree=gsub(' ','_',data$host_species)

## betacovs only
data=data[which(data$virus_genus=='Betacoronavirus'),]
data=data.frame(data)

## betacov
data$betacov=1

## aggregate
data=data[!duplicated(data$host_species),]

## load mammal supertree
setwd("~/Desktop/virionette/04_predictors")
tree=readRDS('Full Supertree.rds')

## data
tdata=data.frame(tree=tree$tip.label)

## merge
data=merge(tdata,data,by='tree',all=T)
rm(tdata)

## fix betacov
data$betacov=replace_na(data$betacov,0)

## load bat supertree
btree=readRDS('bat-supertree_clean.rds')

## get bats
data$bats=ifelse(data$tree%in%btree$tip.label,'bats','other')

## clean
rm(btree)

## load in cites
cites=read.csv('Citations.csv',header=T)
cites$X=NULL
cites$tree=gsub(' ','_',cites$name)
cites$name=NULL

## merge 
data=merge(data,cites,by='tree')
rm(cites)

## merge into phylogeny order
data=data[match(tree$tip.label,data$tree),]

## fix names
data$treenames=data$tree
data$tree=NULL

## merge with caper
library(caper)
cdata=comparative.data(phy=tree,data=data,names.col=treenames,vcv=T,na.omit=F,warn.dropped=T)

## add correct names
cdata$data$treenames=rownames(cdata$data)
cdata$data$Species=rownames(cdata$data)

## sqrt cites
cdata$data$scites=sqrt(cdata$data$cites)

## subset to bats
bdata=cdata[which(cdata$data$bats=='bats'),]

## Holm rejection procedure
HolmProcedure <- function(pf,FWER=0.5){
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

## bat-only phylofactor, uncorrected
set.seed(1)
bat_pf=gpf(Data=bdata$data,tree=bdata$phy,
       frmla.phylo=betacov~phylo,
       family=binomial,algorithm='phylo',nfactors=3)
bat_keep=HolmProcedure(bat_pf)
pf.tree(bat_pf,factors=1:bat_keep,size=0.1,layout="circular")$ggplot

## refit to correct clades
set.seed(1)
bat_pf=gpf(Data=bdata$data,tree=bdata$phy,
           frmla.phylo=betacov~phylo,
           family=binomial,algorithm='phylo',nfactors=bat_keep)

## bat-only phylofactor, corrected for cites
set.seed(1)
cbat_pf=gpf(Data=bdata$data,tree=bdata$phy,
           frmla.phylo=betacov~phylo,
           family=binomial,algorithm='phylo',nfactors=3,
           weights=sqrt(bdata$data$cites))
cbat_keep=HolmProcedure(cbat_pf)

## all mammal phylofactor, uncorrected
set.seed(1)
mam_pf=gpf(Data=cdata$data,tree=cdata$phy,
           frmla.phylo=betacov~phylo,
           family=binomial,algorithm='phylo',nfactors=3)
mam_keep=HolmProcedure(mam_pf)
pf.tree(mam_pf,factors=1:mam_keep,size=0.1,layout="circular")$ggplot

## refit to correct clades
set.seed(1)
mam_pf=gpf(Data=cdata$data,tree=cdata$phy,
           frmla.phylo=betacov~phylo,
           family=binomial,algorithm='phylo',nfactors=mam_keep)

## all mammal phylofactor, corrected for cites
set.seed(1)
cmam_pf=gpf(Data=cdata$data,tree=cdata$phy,
            frmla.phylo=betacov~phylo,
            family=binomial,algorithm='phylo',nfactors=3,
            weights=sqrt(cdata$data$cites))
cmam_keep=HolmProcedure(cmam_pf)

## BeckerBatsUncorrected.csv
set=data.frame(host_species=bdata$data$treenames,
               preds=predict(bat_pf,type='response'))

## write
setwd("~/Desktop/becker-betacov")
write.csv(set,'BeckerBatsUncorrected.csv')

## BeckerBatsCitations.csv
set=data.frame(host_species=bdata$data$treenames,
               preds=ifelse(HolmProcedure(cbat_pf)==0,NA,predict(cbat_pf,type='response')))

## write
setwd("~/Desktop/becker-betacov")
write.csv(set,'BeckerBatsCitations.csv')

## BeckerMammalUncorrected.csv
set=data.frame(host_species=cdata$data$treenames,
               preds=predict(mam_pf,type='response'))

## write
setwd("~/Desktop/becker-betacov")
write.csv(set,'BeckerMammalUncorrected.csv')

## BeckerMammalCitations.csv
set=data.frame(host_species=cdata$data$treenames,
               preds=ifelse(HolmProcedure(cmam_pf)==0,NA,predict(cmam_pf,type='response')))

## write
setwd("~/Desktop/becker-betacov")
write.csv(set,'BeckerMammalCitations.csv')