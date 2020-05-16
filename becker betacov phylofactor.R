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
library(plotrix)

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

## trim
bb=data[c('tree','bats')]
setwd("~/Desktop")
write.csv(bb,'bats vs other.csv')

## clean
rm(btree)

## load in cites
cites=read.csv('Citations.csv',header=T)
cites$X=NULL
cites$tree=gsub(' ','_',cites$name)
cites$name=NULL

## merge 
data=merge(data,cites,by='tree',all.x=T)
rm(cites)
table(is.na(data$cites))

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

## label
cdata$data$label=cdata$data$treenames
cdata$data$betafac=factor(cdata$data$betacov)

## subset to bats
bdata=cdata[which(cdata$data$bats=='bats'),]

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
pf.tree(cbat_pf,factors=1:cbat_keep,size=0.1,layout="circular")$ggplot

## refit to correct clades
set.seed(1)
cbat_pf=gpf(Data=bdata$data,tree=bdata$phy,
           frmla.phylo=betacov~phylo,
           family=binomial,algorithm='phylo',nfactors=cbat_keep,
           weights=sqrt(bdata$data$cites))

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

## load in taxonomy
setwd("~/Desktop/becker-betacov")
taxonomy=read.csv('mammal taxonomy.csv',header=T)
taxonomy$X=NULL
taxonomy$Sp=NULL
taxonomy=taxonomy[!duplicated(taxonomy$hGenus),]

## get tree data
tdata=cdata$data
tdata=tdata[c('treenames')]

## get genus
tdata$hGenus=sapply(strsplit(tdata$treenames,'_'),function(x) x[1])

## merge
test=merge(tdata,taxonomy,all.x=T,by='hGenus')
a=test[is.na(test$hFamily),]
length(unique(a$hGenus))

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

## combine tree and data
library(treeio)
all_tree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")
bat_tree=treeio::full_join(as.treedata(bdata$phy),bdata$data,by="label")

## set x max
plus=5

## make bat base
bat_base=ggtree(bat_tree,size=0.1,layout="fan")
bat_base=bat_base$data

## tips only
bat_base=bat_base[which(bat_base$isTip==T),]

## make data frame
bat_preds=data.frame(x=bat_base$x,
                 y=bat_base$y,
                 yend=bat_base$y,
                 xend=ifelse(bat_base$betacov==0,unique(bat_base$x)[1],
                             unique(bat_base$x)[1]+plus),
                 betafac=factor(bat_base$betacov))

## repeat for mammals
mam_base=ggtree(all_tree,size=0.1,layout="fan")
mam_base=mam_base$data
mam_base=mam_base[which(mam_base$isTip==T),]

## make data frame
mam_preds=data.frame(x=mam_base$x,
                     y=mam_base$y,
                     yend=mam_base$y,
                     xend=ifelse(mam_base$betacov==0,unique(mam_base$x)[1],
                                 unique(mam_base$x)[1]+plus*2),
                     betafac=factor(mam_base$betacov))

## make phylofactor figure
library(patchwork)

## specify lines
lwd=0.1

## fix palette
AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
afun=function(x){
  a=AlberColours[1:x]
  return(a)
}

## bat raw
p1=pf.tree(bat_pf,factors=1:bat_keep,size=lwd,color.fcn=afun,
           layout="circular",alphas=0.5)$ggplot+
  ggtitle('bat phylofactor')+
  geom_segment(data=bat_preds,
               aes(x=x,y=y,xend=xend,yend=yend,
                   colour=betafac),size=0.25)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c('white','black'))

## bat citations
p2=pf.tree(cbat_pf,factors=1:cbat_keep,size=lwd,color.fcn=afun,
           layout="circular",alphas=0.5)$ggplot+
  ggtitle('bat phylofactor (with citations)')+
  geom_segment(data=bat_preds,
               aes(x=x,y=y,xend=xend,yend=yend,
                   colour=betafac),size=0.25)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c('white','black'))

## mammal raw
p3=pf.tree(mam_pf,factors=1:mam_keep,size=lwd/4,color.fcn=afun,
           layout="circular",alphas=rep(0.5,mam_keep))$ggplot+
  ggtitle('mammal phylofactor')+
  geom_segment(data=mam_preds,
               aes(x=x,y=y,xend=xend,yend=yend,
                   colour=betafac),size=0.25)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c('white','black'))

## mammal citations
p4=ggtree(cdata$phy,size=lwd/2,layout="circular")+
  ggtitle('mammal phylofactor (with citations)')+
  geom_segment(data=mam_preds,
               aes(x=x,y=y,xend=xend,yend=yend,
                   colour=betafac),size=0.25)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c('white','black'))

## combine
setwd("~/Dropbox (Personal)/GBatNet/prediction phylo")
png("phylo betacov.png",width=7,height=7,units="in",res=300)
p1+p2+p3+p4
dev.off()

