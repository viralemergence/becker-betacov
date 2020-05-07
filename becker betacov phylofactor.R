## phylofactor of betacoronavirus host status
## danbeck@iu.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(ape)
library(phylofactor)
library(tidyverse)
library(data.table)
library(ggtree)
library(plyr)

## load clean data
setwd("~/Desktop/cleanbats_betacov/clean data")
data=read.csv('bat-phylo-traits_clean.csv',header=T)

## load bat supertree
btree=readRDS('bat-supertree_clean.rds')

## reduce data down to virus, species, taxonomy
bdata=data[c('clean_hostnames','betacov','sarbecov')]

## get genus
bdata$genus=sapply(strsplit(as.character(bdata$clean_hostnames),'_'),function(x) x[1])

## load taxonomy
setwd("~/Desktop/becker-betacov")
taxa=read.csv('bat taxonomy.csv',header=T)
taxa$X=NULL
taxa$hOrder=NULL
names(taxa)=c('family','genus')

## merge into bdata
bdata=merge(bdata,taxa,by='genus',all.x=T)
rm(taxa)

## fix missing families
bdata$family=as.character(bdata$family)
bdata$family2=revalue(bdata$genus,
                      c('Aproteles'='Pteropodidae',
                        'Paracoelops'='Hipposideridae'))
bdata$family=ifelse(is.na(bdata$family),bdata$family2,bdata$family)
bdata$family2=NULL

## taxa
bdata$taxonomy=with(bdata,paste(family,genus,clean_hostnames,sep='; '))

## merge into phylogeny order
bdata=bdata[match(btree$tip.label,bdata$clean_hostnames),]

## merge with caper
library(caper)
cdata=comparative.data(phy=btree,data=bdata,names.col=clean_hostnames,vcv=T,na.omit=F,warn.dropped=T)

## add correct names
cdata$data$clean_hostnames=rownames(cdata$data)
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
keep=HolmProcedure(pf)
pf.tree(pf,factors=1:keep,size=0.1,layout="circular")$ggplot

## refit with retained number of significant clades
pf2=gpf(Data=cdata$data,tree=cdata$phy,
       frmla.phylo=betacov~phylo,
       family=binomial,algorithm='phylo',nfactors=keep)

## set key
setkey(pf2$Data,'Species')

## summarize
pf.taxa(pf2,taxonomy,factor=1)$group1

## get clade
cdata$data$pf_beta=ifelse(cdata$data$clean_hostnames%in%cladeget(pf2,1),'clade','other')

## table
table(cdata$data$pf_beta,cdata$data$betacov)
prop.table(table(cdata$data$pf_beta,cdata$data$betacov))*100

## make data
set=cdata$data

## get predictions
set$Prediction=predict(pf2,type='response')

## export
setwd("~/Desktop/becker-betacov")
write.csv(set,'PhylofactorPredictions.csv')