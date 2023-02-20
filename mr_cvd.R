rm(list=ls())

library(chron)
library(colorspace)
library(mime)
library(munsell)
library(labeling)
library(rlang)
library(stringi)
library(evaluate)
library(highr)
library(markdown)
library(yaml)
library(backports)
library(jsonlite)
library(digest)
library(plyr)
library(reshape2)
library(scales)
library(tibble)
library(lazyeval)
library(RColorBrewer)
library(stringr)
library(knitr)
library(magrittr)
library(checkmate)
library(htmlwidgets)
library(viridisLite)
library(Rcpp)
library(Formula)
library(ggplot2)
library(latticeExtra)
library(acepack)
library(gtable)
library(data.table)
library(htmlTable)
library(viridis)
library(htmltools)
library(base64enc)
library(minqa)
library(RcppEigen)
library(lme4)
library(SparseM)
library(MatrixModels)
library(pbkrtest)
library(quantreg)
library(car)
library(htmlTable)
library(Hmisc)
library(survival)
library(foreign)
library(bitops)
library(caTools)
library(gplots)
library(ROCR)
library(mice)
library(writexl)
library(HardyWeinberg)
library(officer)
library(uuid)
library(compareGroups)
library(nlme)
library(vcd)
library(boot)
library(tibble)
library(haven)
library(icenReg)
library(MASS)
library(sandwich)   
library(lmtest)
library(gam)
library(smoothHR)
library(metafor)
library(DBI)
library(mitools)
library(RcppArmadillo)
library(miceadds)
library(dplyr)
library(estimatr)
library(lubridate)
library(snakecase)
library(janitor)
library(miceadds)
library(fmsb)
library(usethis)

library(devtools)
library(googleAuthR)
library(MendelianRandomization)
library(mr.raps)
library(meta)
library(MRPRESSO)
library(MRInstruments)
library(MRMix)
library(RadialMR)
library(ieugwasr)
library(TwoSampleMR)


RutinesLocals<- "N:/durable/Syntax/ah/routines"
source(file.path(RutinesLocals,"carrega.llibreria.r"))
source(file.path(RutinesLocals,"merge2.r"))
source(file.path(RutinesLocals,"fix2.r"))
source(file.path(RutinesLocals,"table2.r"))
source(file.path(RutinesLocals,"subset2.r"))
source(file.path(RutinesLocals,"format2.r"))
source(file.path(RutinesLocals,"order2.r"))
source(file.path(RutinesLocals,"intervals.r"))


### GUAPAS ###
##############

guapa<-function(x)
{
  redondeo<-ifelse(abs(x)<0.00001,signif(x,1),
                   ifelse(abs(x)<0.0001,signif(x,1),
                          ifelse(abs(x)<0.001,signif(x,1),
                                 ifelse(abs(x)<0.1,sprintf("%.3f",round(x,3)),
                                        ifelse(abs(x)<1,sprintf("%.2f",round(x,2)),
                                               ifelse(abs(x)<10,sprintf("%.2f",round(x,2)),
                                                      ifelse(abs(x)<100,sprintf("%.1f",round(x,1)),
                                                             ifelse(abs(x)>=100,round(x,0)))))))))
  return(redondeo)
}

ic_guapa<-function(x,y,z)
{
  ic<-paste(x," [",y,"; ",z,"]",sep="")
  return(ic)
}

pval_guapa<-function(x)
{
  pval<-ifelse(x<0.00001,"<0.00001",
               ifelse(x<0.001,"<0.001",round(x,3)))
  return(pval)
}

mean_ic_guapa <- function(x, na.rm=FALSE) 
{
  if (na.rm) x <- na.omit(x)
  se<-sqrt(var(x)/length(x))
  z<-qnorm(1-0.05/2)
  media<-mean(x)
  ic95a<-guapa(media-(z*se))
  ic95b<-guapa(media+(z*se))
  media<-guapa(media)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

mean_sd_guapa <- function(x) 
{
  media<-guapa(mean(x, na.rm=TRUE))
  sd<-guapa(sd(x, na.rm=TRUE))
  end<-paste(media," (",sd,")",sep="")
  return(end)
}

beta_se_ic_guapa <- function(x, y) 
{
  z<-qnorm(1-0.05/2)
  ic95a<-guapa(x-(z*y))
  ic95b<-guapa(x+(z*y))
  media<-guapa(x)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa2 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-round(exp(x),5)
  ic95a<-round(exp(x-(z*y)),5)
  ic95b<-round(exp(x+(z*y)),5)
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }

# Non-linear MR function
# source("N:/data/durable/Projects/Magnus_MR_BMI/R/Old/nlme_summ_aes.R")

setwd("N:/durable/Projects/Hernaez_MR_CVD")


#########################################
### CALCULATION OF GRS - DATA SOURCES ###
#########################################

# CAD: van der Harst P et al., Circ Res, 2018 (n=547261, European ancestry (UKBB + CARDIoGRAMC4D), 148 SNPs)


#############################
### GENERATION OF CAD GRS ###
#############################

# MoBa SELECTED SNPs ON CAD #

dat<-as.data.frame(fread("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/cad/cad_grs.raw",header=TRUE,sep="\t",sep2="\t"))
dat<-dat[,c(2,7:dim(dat)[2])]


# SNPs RELATED TO CAD #

cad_ma<-read.csv2("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/cad/gwas_cad_varderharst2018.csv",header=TRUE,
                  sep=";",dec=".")
names(cad_ma)<-tolower(names(cad_ma))
cad_ma<-rename.vars(cad_ma,
                    from=c("effect_allele","eaf"),
                    to=c("tested_allele_ma","maf_ma"))
cad_ma$tested_allele_ma<-tolower(cad_ma$tested_allele_ma)


### CHECK THAT THE 141 SNPs IN DAT BELONG TO THE META-ANALYSIS RESULTS ###

bbb<-as.character(sort(cad_ma$rsid,decreasing=FALSE))
dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))
write.table(ccc,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/descriptive/snps_cad_1smr.csv",
            sep=";", col.names=FALSE ,row.names=FALSE)
length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND cad_ma MATCH


# TABLE WITH NUMBER OF MISSING VALUES OF ALL VARIABLES IN dat #

vars01<-names(dat)
tab_na<-NULL

for(i in 1:length(vars01))
  
{
  x<-length(which(is.na(dat[,vars01[i]])))
  tab_na<-as.data.frame(rbind(tab_na,cbind(x)))
}

colnames(tab_na)<-c("NAs")
rownames(tab_na)<-vars01
table(tab_na$NAs)


# EXTRACT THE INFORMATION ON THE ALLELE TESTED IN MOBA #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
cad_ma<-merge2(cad_ma,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
cad_ma<-na.omit(cad_ma)


# SNP FLIPPING IF THE ALLELE ORDER IN MoBa AND THE META-ANALYSIS DIFFER #

cad_ma$same_coding<-with(cad_ma,ifelse(tested_allele_moba==tested_allele_ma,1,
                                       ifelse(tested_allele_moba!=tested_allele_ma,0,NA)))

cad_ma$flip<-with(cad_ma,ifelse(beta>=0 & same_coding==1,0,
                                ifelse(beta<0 & same_coding==0,0,
                                       ifelse(beta>=0 & same_coding==0,1,
                                              ifelse(beta<0 & same_coding==1,1,NA)))))

cad_ma$maf<-with(cad_ma,ifelse(same_coding==1,maf_ma,
                               ifelse(same_coding==0,1-maf_ma,NA)))
cad_ma$coef<-with(cad_ma,ifelse(beta>0,beta,
                                ifelse(beta<0,-beta,NA)))


### CALCULATION OF DM GRS ###

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("IID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-cad_ma[cad_ma$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-cad_ma[cad_ma$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$cad_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","cad_grs")]

save(dat,file="N:/data/durable/Projects/Hernaez_MR_CVD/source_files/cad_grs.RData")


######################################
### GENERATION OF WORKING DATABASE ###
######################################

### MERGING GRSs ###

setwd("N:/durable/Projects/Hernaez_MR_CVD")

load("./source_files/cad_grs.RData")
cad_grs<-dat
dat$genotype<-1
genotype<-dat[,c("id","genotype")]

load("N:/durable/Projects/Hernaez_MR_BMI/R/MoBa_raw.RData")
dat$sentrixid_mom<-NULL
dat$sentrixid_dad<-NULL

moms<-spss.get("N:/durable/RAW/3.MoBa_genetics/key/2022_03_31_MoBaGeneticsTot_Mother_PDB2374sav.sav",
               use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(moms)<-tolower(names(moms))
moms$m_id_2374<-gsub(" ", "", moms$m_id_2374)
moms$sentrixid_mom<-gsub(" ", "", moms$sentrix_id)
moms$batch_mom<-gsub(" ", "", moms$batch)
moms<-moms[,c("m_id_2374","sentrixid_mom","batch_mom")]
names(genotype)<-c("sentrixid_mom","genotype")
moms<-merge2(moms,genotype,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
moms$genotype<-with(moms,ifelse(is.na(genotype),0,genotype))
moms<-moms[order(moms$m_id_2374,-abs(moms$genotype)),]
moms<-moms[!duplicated(moms$m_id_2374),]
moms$genotype<-NULL
dat<-merge2(dat,moms,by.id=c("m_id_2374"),all.x=TRUE,sort=FALSE)

dads<-spss.get("N:/durable/RAW/3.MoBa_genetics/key/2022_03_31_MoBaGeneticsTot_Father_PDB2374sav.sav",
               use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(dads)<-tolower(names(dads))
dads$f_id_2374<-gsub(" ", "", dads$f_id_2374)
dads$sentrixid_dad<-gsub(" ", "", dads$sentrix_id)
dads$batch_dad<-gsub(" ", "", dads$batch)
dads<-dads[,c("f_id_2374","sentrixid_dad","batch_dad")]
names(genotype)<-c("sentrixid_dad","genotype")
dads<-merge2(dads,genotype,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)
dads$genotype<-with(dads,ifelse(is.na(genotype),0,genotype))
dads<-dads[order(dads$f_id_2374,-abs(dads$genotype)),]
dads<-dads[!duplicated(dads$f_id_2374),]
dads$genotype<-NULL
dat<-merge2(dat,dads,by.id=c("f_id_2374"),all.x=TRUE,sort=FALSE)


cad_grs<-rename.vars(cad_grs,from=c("id","cad_grs"),to=c("sentrixid_mom","cad_grs_mom"))
dat<-merge2(dat,cad_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
cad_grs<-rename.vars(cad_grs,from=c("sentrixid_mom","cad_grs_mom"),to=c("sentrixid_dad","cad_grs_dad"))
dat<-merge2(dat,cad_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

cad_grs<-NULL
length(which(!is.na(dat$cad_grs_mom)))
length(which(!is.na(dat$cad_grs_dad)))


### MERGING QC INDICATORS ###

qc_indic<-as.data.frame(read.delim("N:/durable/Projects/Hernaez_MR_BMI/Stata/Marker_QC_Ind_Alex_group.txt",
                                   header=TRUE,sep=","))
qc_indic$V2<-NULL
qc_indic<-rename.vars(qc_indic,from=c("V1"),to=c("sentrixid_mom"))
qc_indic$qc_gen_mom<-1
dat<-merge2(dat,qc_indic,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
qc_indic<-rename.vars(qc_indic,from=c("sentrixid_mom","qc_gen_mom"),to=c("sentrixid_dad","qc_gen_dad"))
dat<-merge2(dat,qc_indic,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)
dat$qc_gen_mom<-with(dat,ifelse(is.na(qc_gen_mom),0,qc_gen_mom))
dat$qc_gen_dad<-with(dat,ifelse(is.na(qc_gen_dad),0,qc_gen_dad))

length(which(!is.na(dat$batch_mom)))
length(which(!is.na(dat$batch_dad)))


### BMI AND HEIGHT MEASURED VALUES ###

dat$bmi_mom<-dat$aa85/((dat$aa87/100)^2)
dat$bmi_dad<-dat$aa89/((dat$aa88/100)^2)

dat$bmi_mom<-with(dat,ifelse(bmi_mom<13 | bmi_mom>60,NA,bmi_mom))
dat$bmi_dad<-with(dat,ifelse(bmi_dad<13 | bmi_dad>60,NA,bmi_dad))


### HIGHEST EDUCATIONAL LEVEL, COMPLETED OR ONGOING ###

# Mothers #

dat$aa1125<-with(dat,ifelse(!is.na(aa1124) & !is.na(aa1125) & (aa1125>aa1124),aa1125,
                            ifelse(!is.na(aa1124) & !is.na(aa1125) & (aa1125<aa1124),NA,aa1125)))
dat$aa1128x<-with(dat,ifelse(aa1128==1 | aa1129==1,1,0))
dat$aa1128x<-with(dat,ifelse(is.na(aa1128x),0,aa1128x))

dat$edu_mom<-with(dat,ifelse(aa1124==1,1,
                             ifelse(aa1124==2,2,
                                    ifelse(aa1124==3,3,
                                           ifelse(aa1124==4,4,
                                                  ifelse(aa1124==5,5,
                                                         ifelse(aa1124==6,6,NA)))))))
dat$edu_mom<-with(dat,ifelse(is.na(aa1125),edu_mom,
                             ifelse(aa1125==1,1,
                                    ifelse(aa1125==2,2,
                                           ifelse(aa1125==3,3,
                                                  ifelse(aa1125==4,4,
                                                         ifelse(aa1125==5,5,
                                                                ifelse(aa1125==6,6,NA))))))))
dat$edu_mom<-with(dat,ifelse(is.na(edu_mom),0,edu_mom))

dat$edu_mom<-with(dat,ifelse(edu_mom==0 & aa1128x==1,3,
                             ifelse((edu_mom>0 & edu_mom<3) & aa1128x==1,3,
                                    ifelse(edu_mom>=3 & aa1128x==1,edu_mom,
                                           ifelse(edu_mom==0 & aa1128x==0,NA,
                                                  ifelse((edu_mom>0 & edu_mom<3) & aa1128x==0,edu_mom,
                                                         ifelse(edu_mom>=3 & aa1128x==0,edu_mom,NA)))))))

dat$eduyears_mom<-with(dat,ifelse(edu_mom==1,10,
                                  ifelse(edu_mom==2,10,
                                         ifelse(edu_mom==3,13,
                                                ifelse(edu_mom==4,13,
                                                       ifelse(edu_mom==5,19,
                                                              ifelse(edu_mom==6,20,
                                                                     ifelse(is.na(edu_mom),NA,NA))))))))

dat$edu_mom<-with(dat,ifelse(edu_mom==1,1,
                             ifelse(edu_mom==2,1,
                                    ifelse(edu_mom==3,2,
                                           ifelse(edu_mom==4,2,
                                                  ifelse(edu_mom==5,3,
                                                         ifelse(edu_mom==6,4,
                                                                ifelse(is.na(edu_mom),NA,NA))))))))

# Fathers #

dat$aa1127<-with(dat,ifelse(!is.na(aa1126) & !is.na(aa1127) & (aa1127>aa1126),aa1127,
                            ifelse(!is.na(aa1126) & !is.na(aa1127) & (aa1127<aa1126),NA,aa1127)))
dat$aa1130x<-with(dat,ifelse(aa1130==1 | aa1131==1,1,0))
dat$aa1130x<-with(dat,ifelse(is.na(aa1130x),0,aa1130x))

dat$edu_dad<-with(dat,ifelse(aa1126==1,1,
                             ifelse(aa1126==2,2,
                                    ifelse(aa1126==3,3,
                                           ifelse(aa1126==4,4,
                                                  ifelse(aa1126==5,5,
                                                         ifelse(aa1126==6,6,NA)))))))
dat$edu_dad<-with(dat,ifelse(is.na(aa1127),edu_dad,
                             ifelse(aa1127==1,1,
                                    ifelse(aa1127==2,2,
                                           ifelse(aa1127==3,3,
                                                  ifelse(aa1127==4,4,
                                                         ifelse(aa1127==5,5,
                                                                ifelse(aa1127==6,6,NA))))))))
dat$edu_dad<-with(dat,ifelse(is.na(edu_dad),0,edu_dad))

dat$edu_dad<-with(dat,ifelse(edu_dad==0 & aa1130x==1,3,
                             ifelse((edu_dad>0 & edu_dad<3) & aa1130x==1,3,
                                    ifelse(edu_dad>=3 & aa1130x==1,edu_dad,
                                           ifelse(edu_dad==0 & aa1130x==0,NA,
                                                  ifelse((edu_dad>0 & edu_dad<3) & aa1130x==0,edu_dad,
                                                         ifelse(edu_dad>=3 & aa1130x==0,edu_dad,NA)))))))

dat$eduyears_dad<-with(dat,ifelse(edu_dad==1,10,
                                  ifelse(edu_dad==2,10,
                                         ifelse(edu_dad==3,13,
                                                ifelse(edu_dad==4,13,
                                                       ifelse(edu_dad==5,19,
                                                              ifelse(edu_dad==6,20,
                                                                     ifelse(is.na(edu_dad),NA,NA))))))))
dat$edu_dad<-with(dat,ifelse(edu_dad==1,1,
                             ifelse(edu_dad==2,1,
                                    ifelse(edu_dad==3,2,
                                           ifelse(edu_dad==4,2,
                                                  ifelse(edu_dad==5,3,
                                                         ifelse(edu_dad==6,4,
                                                                ifelse(is.na(edu_dad),NA,NA))))))))


### DEFINITION OF PARITY (number of previous deliveries) ###

dat$parity<-with(dat,ifelse(paritet_5==0,0,
                            ifelse(paritet_5==1,1,
                                   ifelse(paritet_5==2,2,
                                          ifelse(paritet_5==3,3,
                                                 ifelse(paritet_5==4,3,NA))))))

length(which(is.na(dat$parity)))


### DEFINITION OF EVER SMOKERS (smkinit) ###

# Mothers #
dat$smkinit_mom<-with(dat,ifelse(aa1355==2 | aa1356>1 | aa1357!=0 | aa1358!=0 | aa1359>1 | aa1360!=0 | aa1361!=0 | 
                                   !is.na(aa1362) | aa1363==2 | !is.na(aa1364) | !is.na(aa1365),1,0))
dat$smkinit_mom<-with(dat,ifelse(is.na(smkinit_mom),0,smkinit_mom))
dat$smkinit_mom<-with(dat,ifelse(is.na(aa1355) & is.na(aa1356) & is.na(aa1357) & is.na(aa1358) & is.na(aa1359) & is.na(aa1360) & is.na(aa1361) & 
                                   is.na(aa1362) & is.na(aa1363) & is.na(aa1364) & is.na(aa1365),NA,smkinit_mom))

# Fathers #
dat$smkinit_dad<-with(dat,ifelse(aa1353==2 | aa1354==2 | ff214==2 | ff215>1 | ff216!=0 | ff217!=0 | ff218>1 | ff219!=0 | ff220!=0,1,0))
dat$smkinit_dad<-with(dat,ifelse(is.na(smkinit_dad),0,smkinit_dad))
dat$smkinit_dad<-with(dat,ifelse(is.na(aa1353) & is.na(aa1354) & is.na(ff214) & is.na(ff215) & is.na(ff216) & is.na(ff217) & is.na(ff218) & is.na(ff219) & is.na(ff220),NA,smkinit_dad))


### MERGING PRINCIPAL COMPONENTS ###

pcs<-as.data.frame(read.delim("N:/durable/RAW/3.MoBa_genetics/p471/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt",
                              header=TRUE,sep="\t"))[,c("SENTRIXID","genotyping_batch_num","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                                                        "PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")]
names(pcs)<-c("sentrixid_mom","genotype_batch_mom","pc01_mom","pc02_mom","pc03_mom","pc04_mom","pc05_mom","pc06_mom","pc07_mom","pc08_mom","pc09_mom","pc10_mom",
              "pc11_mom","pc12_mom","pc13_mom","pc14_mom","pc15_mom","pc16_mom","pc17_mom","pc18_mom","pc19_mom","pc20_mom")
pcs$sentrixid_mom<-as.character(pcs$sentrixid_mom)
dat<-merge2(dat,pcs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)

pcs<-rename.vars(pcs,
                 from=c("sentrixid_mom","genotype_batch_mom","pc01_mom","pc02_mom","pc03_mom","pc04_mom","pc05_mom","pc06_mom","pc07_mom","pc08_mom","pc09_mom","pc10_mom",
                        "pc11_mom","pc12_mom","pc13_mom","pc14_mom","pc15_mom","pc16_mom","pc17_mom","pc18_mom","pc19_mom","pc20_mom"),
                 to=c("sentrixid_dad","genotype_batch_dad","pc01_dad","pc02_dad","pc03_dad","pc04_dad","pc05_dad","pc06_dad","pc07_dad","pc08_dad","pc09_dad","pc10_dad",
                      "pc11_dad","pc12_dad","pc13_dad","pc14_dad","pc15_dad","pc16_dad","pc17_dad","pc18_dad","pc19_dad","pc20_dad"))
pcs$sentrixid_dad<-as.character(pcs$sentrixid_dad)
dat<-merge2(dat,pcs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

length(which(!is.na(dat$pc01_mom)))
length(which(!is.na(dat$pc02_mom)))

length(which(!is.na(dat$pc01_dad)))
length(which(!is.na(dat$pc02_dad)))


### DEFINITION OF SUBFERTILITY ###

dat$exclude<-with(dat,ifelse(versjon_skjema1_tbl1=="",1,0))
dat$art<-with(dat,ifelse(!is.na(art) & art>0,1,0))

dat$subf<-with(dat,ifelse(aa48<=12,0,
                          ifelse(aa48>12,1,NA)))
dat$subf<-with(dat,ifelse(art==0,subf,
                          ifelse(art==1,1,NA)))
dat$subf<-with(dat,ifelse(exclude==0 & is.na(subf),0,
                          ifelse(exclude==0 & subf==0,0,
                                 ifelse(exclude==0 & subf==1,1,
                                        ifelse(exclude==1 & is.na(subf),NA,
                                               ifelse(exclude==1 & subf==0,NA,
                                                      ifelse(exclude==1 & subf==1,NA,NA)))))))

dat$subf_12plus<-with(dat,ifelse(aa48<12,0,
                                 ifelse(aa48>=12,1,NA)))
dat$subf_12plus<-with(dat,ifelse(art==0,subf_12plus,
                                 ifelse(art==1,1,NA)))
dat$subf_12plus<-with(dat,ifelse(exclude==0 & is.na(subf_12plus),0,
                                 ifelse(exclude==0 & subf_12plus==0,0,
                                        ifelse(exclude==0 & subf_12plus==1,1,
                                               ifelse(exclude==1 & is.na(subf_12plus),NA,
                                                      ifelse(exclude==1 & subf_12plus==0,NA,
                                                             ifelse(exclude==1 & subf_12plus==1,NA,NA)))))))


### DEFINITION OF HISTORY OF MISCARRIAGE ###

vars<-c("aa95","aa96","aa101","aa102","aa107","aa108","aa113","aa114",
        "aa119","aa120","aa125","aa126","aa131","aa132","aa137","aa138",
        "aa143","aa144","aa149","aa150")

for(i in 1:length(vars))
  
{
  dat[,vars[i]]<-with(dat,ifelse(is.na(dat[,vars[i]]),0,dat[,vars[i]]))
}

dat$miscarriage<-with(dat,ifelse(aa95==2 & aa96<20,1,
                                 ifelse(aa101==2 & aa102<20,1,
                                        ifelse(aa107==2 & aa108<20,1,
                                               ifelse(aa113==2 & aa114<20,1,
                                                      ifelse(aa119==2 & aa120<20,1,
                                                             ifelse(aa125==2 & aa126<20,1,
                                                                    ifelse(aa131==2 & aa132<20,1,
                                                                           ifelse(aa137==2 & aa138<20,1,
                                                                                  ifelse(aa143==2 & aa144<20,1,
                                                                                         ifelse(aa149==2 & aa150<20,1,0)))))))))))


### DEFINITION OF HISTORY OF PREECLAMPSIA/ECLAMPSIA ###

# All pregnancies of the MoBa mothers #

mbrn<-spss.get("N:/durable/RAW/1.MoBa_questionnaires_MBRN/2.MBRN_Birth_Registry/If_you_want_records_MoBa_parents_siblings_children/MFRdata_202257_2020Q4.sav",use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(mbrn)<-tolower(names(mbrn))
mbrn$birthyear_mom<-with(mbrn,faar-mors_alder)
mbrn$birthyear_dad<-with(mbrn,faar-fars_alder)


### PREECLAMPSIA / ECLAMPSIA ###

vars<-c("preekl","preekltidl","eklampsi","hellp","hypertensjon_alene","hypertensjon_kronisk","svlen","diabetes_mellitus",
        "spabort_12_5","spabort_23_5","dodfodte_5","dodkat","fstart")

for(i in 1:length(vars))
{
  mbrn[,vars[i]]<-with(mbrn,ifelse(is.na(mbrn[,vars[i]]),0,mbrn[,vars[i]]))
}

mbrn$eclampsia<-with(mbrn,ifelse(preekl>0,1,
                                         ifelse(preekltidl>0,1,
                                                ifelse(eklampsi>0,1,
                                                       ifelse(hellp>0,1,0)))))
mbrn$hta_preg<-with(mbrn,ifelse(preekl>0,1,
                                         ifelse(preekltidl>0,1,
                                                ifelse(eklampsi>0,1,
                                                       ifelse(hellp>0,1,
                                                              ifelse(hypertensjon_alene>0,1,0))))))
mbrn$hta_chronic<-mbrn$hypertensjon_kronisk


### PRETERM BIRTH ###                                               

mbrn$preterm<-with(mbrn,ifelse(svlen<37,1,0))
mbrn$preterm_sp<-with(mbrn,ifelse(preterm==1 & fstart==1,1,
                                  ifelse(preterm==1 & fstart==2,NA,
                                         ifelse(preterm==1 & fstart==3,NA,
                                                ifelse(preterm==0,0,NA)))))

### GESTATIONAL DIABETES ###                                               

mbrn$gdm<-with(mbrn,ifelse(diabetes_mellitus==4,1,0))


### STILLBIRTHS (ACCORDING TO DODKAT = 7/8/9) ###

mbrn$stillbirth<-with(mbrn,ifelse(dodkat==7,1,
                                  ifelse(dodkat==8,1,
                                         ifelse(dodkat==9,1,0))))


### MISCARRIAGES / STILLBIRTHS ACCORDING TO MBRN ###

mbrn$miscarriage12_mbrn<-with(mbrn,ifelse(mbrn$spabort_12_5>0,1,0))
mbrn$miscarriage24_mbrn<-with(mbrn,ifelse(mbrn$spabort_23_5>0,1,0))
mbrn$miscarriage_mbrn<-with(mbrn,ifelse(mbrn$miscarriage12_mbrn>0 | mbrn$miscarriage24_mbrn>0,1,0))
mbrn$stillbirth24_mbrn<-with(mbrn,ifelse(mbrn$dodfodte_5>0,1,0))


### SMALL AND LARGE FOR GESTATIONAL AGE (SGA, LGA) ### 

# Create unique IDs for children without id ("lopenr_barn" is NA) #

mbrn$series<-1:dim(mbrn)[1]
mbrn$lopenr_barn2<-with(mbrn,ifelse(is.na(lopenr_barn),paste(lopenr_mor,series,sep="_"),lopenr_barn))
mbrn$lopenr_barn<-with(mbrn,ifelse(is.na(lopenr_barn) & is.na(lopenr_mor),NA,lopenr_barn2))
mbrn<-subset2(mbrn,"!is.na(mbrn$lopenr_barn)")

# Erase children with sex not specified (0), uncertain (3) or missing (9) #

mbrn$kjonn<-with(mbrn,ifelse(kjonn==1,1,
                                     ifelse(kjonn==2,2,NA)))
mbrn<-subset2(mbrn,"!is.na(mbrn$kjonn)")

mbrn$svlen2<-with(mbrn,ifelse(svlen<25,24,
                              ifelse(svlen>42,43,svlen)))
mbrn$vekt<-with(mbrn,ifelse(vekt<500,NA,
                              ifelse(vekt>6500,NA,vekt)))

mbrn_boys<-subset2(mbrn,"mbrn$kjonn==1")
mbrn_boys<-mbrn_boys[,c("lopenr_barn","svlen2","vekt")]
mbrn_girls<-subset2(mbrn,"mbrn$kjonn==2")
mbrn_girls<-mbrn_girls[,c("lopenr_barn","svlen2","vekt")]

# SGA and LGA in boys #

mbrn_boys$s24<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==24,mbrn_boys$vekt,NA))
mbrn_boys$s25<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==25,mbrn_boys$vekt,NA))
mbrn_boys$s26<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==26,mbrn_boys$vekt,NA))
mbrn_boys$s27<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==27,mbrn_boys$vekt,NA))
mbrn_boys$s28<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==28,mbrn_boys$vekt,NA))
mbrn_boys$s29<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==29,mbrn_boys$vekt,NA))
mbrn_boys$s30<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==30,mbrn_boys$vekt,NA))
mbrn_boys$s31<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==31,mbrn_boys$vekt,NA))
mbrn_boys$s32<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==32,mbrn_boys$vekt,NA))
mbrn_boys$s33<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==33,mbrn_boys$vekt,NA))
mbrn_boys$s34<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==34,mbrn_boys$vekt,NA))
mbrn_boys$s35<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==35,mbrn_boys$vekt,NA))
mbrn_boys$s36<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==36,mbrn_boys$vekt,NA))
mbrn_boys$s37<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==37,mbrn_boys$vekt,NA))
mbrn_boys$s38<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==38,mbrn_boys$vekt,NA))
mbrn_boys$s39<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==39,mbrn_boys$vekt,NA))
mbrn_boys$s40<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==40,mbrn_boys$vekt,NA))
mbrn_boys$s41<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==41,mbrn_boys$vekt,NA))
mbrn_boys$s42<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==42,mbrn_boys$vekt,NA))
mbrn_boys$s43<-with(mbrn_boys,ifelse(mbrn_boys$svlen2==43,mbrn_boys$vekt,NA))

mbrn_boys$sga24<-as.numeric(ntile(mbrn_boys$s24,10))
mbrn_boys$sga25<-as.numeric(ntile(mbrn_boys$s25,10))
mbrn_boys$sga26<-as.numeric(ntile(mbrn_boys$s26,10))
mbrn_boys$sga27<-as.numeric(ntile(mbrn_boys$s27,10))
mbrn_boys$sga28<-as.numeric(ntile(mbrn_boys$s28,10))
mbrn_boys$sga29<-as.numeric(ntile(mbrn_boys$s29,10))
mbrn_boys$sga30<-as.numeric(ntile(mbrn_boys$s30,10))
mbrn_boys$sga31<-as.numeric(ntile(mbrn_boys$s31,10))
mbrn_boys$sga32<-as.numeric(ntile(mbrn_boys$s32,10))
mbrn_boys$sga33<-as.numeric(ntile(mbrn_boys$s33,10))
mbrn_boys$sga34<-as.numeric(ntile(mbrn_boys$s34,10))
mbrn_boys$sga35<-as.numeric(ntile(mbrn_boys$s35,10))
mbrn_boys$sga36<-as.numeric(ntile(mbrn_boys$s36,10))
mbrn_boys$sga37<-as.numeric(ntile(mbrn_boys$s37,10))
mbrn_boys$sga38<-as.numeric(ntile(mbrn_boys$s38,10))
mbrn_boys$sga39<-as.numeric(ntile(mbrn_boys$s39,10))
mbrn_boys$sga40<-as.numeric(ntile(mbrn_boys$s40,10))
mbrn_boys$sga41<-as.numeric(ntile(mbrn_boys$s41,10))
mbrn_boys$sga42<-as.numeric(ntile(mbrn_boys$s42,10))
mbrn_boys$sga43<-as.numeric(ntile(mbrn_boys$s43,10))

vars<-c("sga24","sga25","sga26","sga27","sga28","sga29","sga30","sga31","sga32","sga33",
        "sga34","sga35","sga36","sga37","sga38","sga39","sga40","sga41","sga42","sga43")

for(i in 1:length(vars))
  
{
  mbrn_boys[,vars[i]]<-with(mbrn_boys,ifelse(is.na(mbrn_boys[,vars[i]]),0,mbrn_boys[,vars[i]]))
}


mbrn_boys$sga_boys<-with(mbrn_boys,ifelse(sga24==1,1,
                                          ifelse(sga25==1,1,
                                                 ifelse(sga26==1,1,
                                                        ifelse(sga27==1,1,
                                                               ifelse(sga28==1,1,
                                                                      ifelse(sga29==1,1,
                                                                             ifelse(sga30==1,1,
                                                                                    ifelse(sga31==1,1,
                                                                                           ifelse(sga32==1,1,
                                                                                                  ifelse(sga33==1,1,
                                                                                                         ifelse(sga34==1,1,
                                                                                                                ifelse(sga35==1,1,
                                                                                                                       ifelse(sga36==1,1,
                                                                                                                              ifelse(sga37==1,1,
                                                                                                                                     ifelse(sga38==1,1,
                                                                                                                                            ifelse(sga39==1,1,
                                                                                                                                                   ifelse(sga40==1,1,
                                                                                                                                                          ifelse(sga41==1,1,
                                                                                                                                                                 ifelse(sga42==1,1,
                                                                                                                                                                        ifelse(sga43==1,1,0)))))))))))))))))))))
mbrn_boys$lga_boys<-with(mbrn_boys,ifelse(sga24==10,1,
                                          ifelse(sga25==10,1,
                                                 ifelse(sga26==10,1,
                                                        ifelse(sga27==10,1,
                                                               ifelse(sga28==10,1,
                                                                      ifelse(sga29==10,1,
                                                                             ifelse(sga30==10,1,
                                                                                    ifelse(sga31==10,1,
                                                                                           ifelse(sga32==10,1,
                                                                                                  ifelse(sga33==10,1,
                                                                                                         ifelse(sga34==10,1,
                                                                                                                ifelse(sga35==10,1,
                                                                                                                       ifelse(sga36==10,1,
                                                                                                                              ifelse(sga37==10,1,
                                                                                                                                     ifelse(sga38==10,1,
                                                                                                                                            ifelse(sga39==10,1,
                                                                                                                                                   ifelse(sga40==10,1,
                                                                                                                                                          ifelse(sga41==10,1,
                                                                                                                                                                 ifelse(sga42==10,1,
                                                                                                                                                                        ifelse(sga43==10,1,0)))))))))))))))))))))
mbrn_boys<-mbrn_boys[,c("lopenr_barn","sga_boys","lga_boys")]


# SGA and LGA in girls #

mbrn_girls$s24<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==24,mbrn_girls$vekt,NA))
mbrn_girls$s25<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==25,mbrn_girls$vekt,NA))
mbrn_girls$s26<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==26,mbrn_girls$vekt,NA))
mbrn_girls$s27<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==27,mbrn_girls$vekt,NA))
mbrn_girls$s28<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==28,mbrn_girls$vekt,NA))
mbrn_girls$s29<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==29,mbrn_girls$vekt,NA))
mbrn_girls$s30<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==30,mbrn_girls$vekt,NA))
mbrn_girls$s31<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==31,mbrn_girls$vekt,NA))
mbrn_girls$s32<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==32,mbrn_girls$vekt,NA))
mbrn_girls$s33<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==33,mbrn_girls$vekt,NA))
mbrn_girls$s34<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==34,mbrn_girls$vekt,NA))
mbrn_girls$s35<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==35,mbrn_girls$vekt,NA))
mbrn_girls$s36<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==36,mbrn_girls$vekt,NA))
mbrn_girls$s37<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==37,mbrn_girls$vekt,NA))
mbrn_girls$s38<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==38,mbrn_girls$vekt,NA))
mbrn_girls$s39<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==39,mbrn_girls$vekt,NA))
mbrn_girls$s40<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==40,mbrn_girls$vekt,NA))
mbrn_girls$s41<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==41,mbrn_girls$vekt,NA))
mbrn_girls$s42<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==42,mbrn_girls$vekt,NA))
mbrn_girls$s43<-with(mbrn_girls,ifelse(mbrn_girls$svlen2==43,mbrn_girls$vekt,NA))

mbrn_girls$sga24<-as.numeric(ntile(mbrn_girls$s24,10))
mbrn_girls$sga25<-as.numeric(ntile(mbrn_girls$s25,10))
mbrn_girls$sga26<-as.numeric(ntile(mbrn_girls$s26,10))
mbrn_girls$sga27<-as.numeric(ntile(mbrn_girls$s27,10))
mbrn_girls$sga28<-as.numeric(ntile(mbrn_girls$s28,10))
mbrn_girls$sga29<-as.numeric(ntile(mbrn_girls$s29,10))
mbrn_girls$sga30<-as.numeric(ntile(mbrn_girls$s30,10))
mbrn_girls$sga31<-as.numeric(ntile(mbrn_girls$s31,10))
mbrn_girls$sga32<-as.numeric(ntile(mbrn_girls$s32,10))
mbrn_girls$sga33<-as.numeric(ntile(mbrn_girls$s33,10))
mbrn_girls$sga34<-as.numeric(ntile(mbrn_girls$s34,10))
mbrn_girls$sga35<-as.numeric(ntile(mbrn_girls$s35,10))
mbrn_girls$sga36<-as.numeric(ntile(mbrn_girls$s36,10))
mbrn_girls$sga37<-as.numeric(ntile(mbrn_girls$s37,10))
mbrn_girls$sga38<-as.numeric(ntile(mbrn_girls$s38,10))
mbrn_girls$sga39<-as.numeric(ntile(mbrn_girls$s39,10))
mbrn_girls$sga40<-as.numeric(ntile(mbrn_girls$s40,10))
mbrn_girls$sga41<-as.numeric(ntile(mbrn_girls$s41,10))
mbrn_girls$sga42<-as.numeric(ntile(mbrn_girls$s42,10))
mbrn_girls$sga43<-as.numeric(ntile(mbrn_girls$s43,10))

vars<-c("sga24","sga25","sga26","sga27","sga28","sga29","sga30","sga31","sga32","sga33",
        "sga34","sga35","sga36","sga37","sga38","sga39","sga40","sga41","sga42","sga43")

for(i in 1:length(vars))
  
{
  mbrn_girls[,vars[i]]<-with(mbrn_girls,ifelse(is.na(mbrn_girls[,vars[i]]),0,mbrn_girls[,vars[i]]))
}


mbrn_girls$sga_girls<-with(mbrn_girls,ifelse(sga24==1,1,
                                             ifelse(sga25==1,1,
                                                    ifelse(sga26==1,1,
                                                           ifelse(sga27==1,1,
                                                                  ifelse(sga28==1,1,
                                                                         ifelse(sga29==1,1,
                                                                                ifelse(sga30==1,1,
                                                                                       ifelse(sga31==1,1,
                                                                                              ifelse(sga32==1,1,
                                                                                                     ifelse(sga33==1,1,
                                                                                                            ifelse(sga34==1,1,
                                                                                                                   ifelse(sga35==1,1,
                                                                                                                          ifelse(sga36==1,1,
                                                                                                                                 ifelse(sga37==1,1,
                                                                                                                                        ifelse(sga38==1,1,
                                                                                                                                               ifelse(sga39==1,1,
                                                                                                                                                      ifelse(sga40==1,1,
                                                                                                                                                             ifelse(sga41==1,1,
                                                                                                                                                                    ifelse(sga42==1,1,
                                                                                                                                                                           ifelse(sga43==1,1,0)))))))))))))))))))))

mbrn_girls$lga_girls<-with(mbrn_girls,ifelse(sga24==10,1,
                                             ifelse(sga25==10,1,
                                                    ifelse(sga26==10,1,
                                                           ifelse(sga27==10,1,
                                                                  ifelse(sga28==10,1,
                                                                         ifelse(sga29==10,1,
                                                                                ifelse(sga30==10,1,
                                                                                       ifelse(sga31==10,1,
                                                                                              ifelse(sga32==10,1,
                                                                                                     ifelse(sga33==10,1,
                                                                                                            ifelse(sga34==10,1,
                                                                                                                   ifelse(sga35==10,1,
                                                                                                                          ifelse(sga36==10,1,
                                                                                                                                 ifelse(sga37==10,1,
                                                                                                                                        ifelse(sga38==10,1,
                                                                                                                                               ifelse(sga39==10,1,
                                                                                                                                                      ifelse(sga40==10,1,
                                                                                                                                                             ifelse(sga41==10,1,
                                                                                                                                                                    ifelse(sga42==10,1,
                                                                                                                                                                           ifelse(sga43==10,1,0)))))))))))))))))))))

mbrn_girls<-mbrn_girls[,c("lopenr_barn","sga_girls","lga_girls")]

mbrn<-merge2(mbrn,mbrn_boys,by.id=c("lopenr_barn"),all.x=TRUE,sort=FALSE)
mbrn<-merge2(mbrn,mbrn_girls,by.id=c("lopenr_barn"),all.x=TRUE,sort=FALSE)

vars<-c("sga_boys","sga_girls","lga_boys","lga_girls")
for(i in 1:length(vars))
{
  mbrn[,vars[i]]<-with(mbrn,ifelse(is.na(mbrn[,vars[i]]),0,mbrn[,vars[i]]))
}

mbrn$sga<-with(mbrn,ifelse(sga_boys==1,1,
                           ifelse(sga_girls==1,1,0)))
mbrn$lga<-with(mbrn,ifelse(lga_boys==1,1,
                           ifelse(lga_girls==1,1,0)))


### Create variables that inform of any event in any pregnancy (included or not in MoBa) ###

mbrn$mobaid_mfr_2374<-with(mbrn,ifelse(mobaid_mfr_2374=="        ",NA,mobaid_mfr_2374))
mbrn$flerfodsel<-with(mbrn,ifelse(is.na(flerfodsel),0,flerfodsel))

# Mothers #

mbrn_mom<-subset2(mbrn,"!is.na(mbrn$lopenr_mor)")

mom_agemax<-mbrn_mom[,c("lopenr_mor","mors_alder")]
mom_agemax<-mom_agemax[order(mom_agemax$lopenr_mor,-abs(mom_agemax$mors_alder)),]
mom_agemax<-mom_agemax[!duplicated(mom_agemax$lopenr_mor),]

mom_yearpregmax<-mbrn_mom[,c("lopenr_mor","faar")]
mom_yearpregmax<-mom_yearpregmax[order(mom_yearpregmax$lopenr_mor,abs(mom_yearpregmax$faar)),]
mom_yearpregmax<-mom_yearpregmax[!duplicated(mom_yearpregmax$lopenr_mor),]

mom_birthyear<-mbrn_mom[,c("lopenr_mor","birthyear_mom")]
mom_birthyear<-mom_birthyear[order(mom_birthyear$lopenr_mor,abs(mom_birthyear$birthyear_mom)),]
mom_birthyear<-mom_birthyear[!duplicated(mom_birthyear$lopenr_mor),]

mom_paritymax<-mbrn_mom[,c("lopenr_mor","paritet_5")]
mom_paritymax<-mom_paritymax[order(mom_paritymax$lopenr_mor,-abs(mom_paritymax$paritet_5)),]
mom_paritymax<-mom_paritymax[!duplicated(mom_paritymax$lopenr_mor),]

mom_eclampsia<-mbrn_mom[,c("lopenr_mor","eclampsia")]
mom_eclampsia<-mom_eclampsia[order(mom_eclampsia$lopenr_mor,-abs(mom_eclampsia$eclampsia)),]
mom_eclampsia<-mom_eclampsia[!duplicated(mom_eclampsia$lopenr_mor),]

mom_hta_preg<-mbrn_mom[,c("lopenr_mor","hta_preg")]
mom_hta_preg<-mom_hta_preg[order(mom_hta_preg$lopenr_mor,-abs(mom_hta_preg$hta_preg)),]
mom_hta_preg<-mom_hta_preg[!duplicated(mom_hta_preg$lopenr_mor),]

mom_hta_chronic<-mbrn_mom[,c("lopenr_mor","hta_chronic")]
mom_hta_chronic<-mom_hta_chronic[order(mom_hta_chronic$lopenr_mor,-abs(mom_hta_chronic$hta_chronic)),]
mom_hta_chronic<-mom_hta_chronic[!duplicated(mom_hta_chronic$lopenr_mor),]

mom_preterm<-mbrn_mom[,c("lopenr_mor","preterm")]
mom_preterm<-mom_preterm[order(mom_preterm$lopenr_mor,-abs(mom_preterm$preterm)),]
mom_preterm<-mom_preterm[!duplicated(mom_preterm$lopenr_mor),]

mom_preterm_sp<-mbrn_mom[,c("lopenr_mor","preterm_sp")]
mom_preterm_sp<-mom_preterm_sp[order(mom_preterm_sp$lopenr_mor,-abs(mom_preterm_sp$preterm_sp)),]
mom_preterm_sp<-mom_preterm_sp[!duplicated(mom_preterm_sp$lopenr_mor),]

mom_gdm<-mbrn_mom[,c("lopenr_mor","gdm")]
mom_gdm<-mom_gdm[order(mom_gdm$lopenr_mor,-abs(mom_gdm$gdm)),]
mom_gdm<-mom_gdm[!duplicated(mom_gdm$lopenr_mor),]

mom_sga<-mbrn_mom[,c("lopenr_mor","sga")]
mom_sga<-mom_sga[order(mom_sga$lopenr_mor,-abs(mom_sga$sga)),]
mom_sga<-mom_sga[!duplicated(mom_sga$lopenr_mor),]

mom_lga<-mbrn_mom[,c("lopenr_mor","lga")]
mom_lga<-mom_lga[order(mom_lga$lopenr_mor,-abs(mom_lga$lga)),]
mom_lga<-mom_lga[!duplicated(mom_lga$lopenr_mor),]

mom_miscarriage12_mbrn<-mbrn_mom[,c("lopenr_mor","miscarriage12_mbrn")]
mom_miscarriage12_mbrn<-mom_miscarriage12_mbrn[order(mom_miscarriage12_mbrn$lopenr_mor,-abs(mom_miscarriage12_mbrn$miscarriage12_mbrn)),]
mom_miscarriage12_mbrn<-mom_miscarriage12_mbrn[!duplicated(mom_miscarriage12_mbrn$lopenr_mor),]

mom_miscarriage24_mbrn<-mbrn_mom[,c("lopenr_mor","miscarriage24_mbrn")]
mom_miscarriage24_mbrn<-mom_miscarriage24_mbrn[order(mom_miscarriage24_mbrn$lopenr_mor,-abs(mom_miscarriage24_mbrn$miscarriage24_mbrn)),]
mom_miscarriage24_mbrn<-mom_miscarriage24_mbrn[!duplicated(mom_miscarriage24_mbrn$lopenr_mor),]

mom_miscarriage_mbrn<-mbrn_mom[,c("lopenr_mor","miscarriage_mbrn")]
mom_miscarriage_mbrn<-mom_miscarriage_mbrn[order(mom_miscarriage_mbrn$lopenr_mor,-abs(mom_miscarriage_mbrn$miscarriage_mbrn)),]
mom_miscarriage_mbrn<-mom_miscarriage_mbrn[!duplicated(mom_miscarriage_mbrn$lopenr_mor),]

mom_stillbirth24_mbrn<-mbrn_mom[,c("lopenr_mor","stillbirth24_mbrn")]
mom_stillbirth24_mbrn<-mom_stillbirth24_mbrn[order(mom_stillbirth24_mbrn$lopenr_mor,-abs(mom_stillbirth24_mbrn$stillbirth24_mbrn)),]
mom_stillbirth24_mbrn<-mom_stillbirth24_mbrn[!duplicated(mom_stillbirth24_mbrn$lopenr_mor),]

mom_stillbirth<-mbrn_mom[,c("lopenr_mor","stillbirth")]
mom_stillbirth<-mom_stillbirth[order(mom_stillbirth$lopenr_mor,-abs(mom_stillbirth$stillbirth)),]
mom_stillbirth<-mom_stillbirth[!duplicated(mom_stillbirth$lopenr_mor),]

mbrn_mom<-mbrn_mom[,c("mobaid_mfr_2374","lopenr_mor","flerfodsel")]
mbrn_mom<-merge2(mbrn_mom,mom_eclampsia,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_hta_preg,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_hta_chronic,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_preterm,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_preterm_sp,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_gdm,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_sga,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_lga,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_miscarriage12_mbrn,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_miscarriage24_mbrn,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_miscarriage_mbrn,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_stillbirth24_mbrn,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_stillbirth,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_agemax,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_yearpregmax,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_birthyear,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_paritymax,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-subset2(mbrn_mom,"!is.na(mbrn_mom$mobaid_mfr_2374)")
names(mbrn_mom)<-c("lopenr_mor","mobaid_mfr_2374","flerfodsel","eclampsia_mom","hta_preg_mom","hta_chronic_mom",
                   "preterm_mom","preterm_sp_mom","gdm_mom","sga_mom","lga_mom",
                   "miscarriage12_mbrn_mom","miscarriage24_mbrn_mom","miscarriage_mbrn_mom",
                   "stillbirth24_mbrn_mom","stillbirth_mom","agemax_mom","yearpregmax_mom","birthyear_mom","paritymax_mom")
mbrn_mom$agemax_mom<-mbrn_mom$yearpregmax_mom-mbrn_mom$birthyear_mom


# Fathers #

mbrn_dad<-subset2(mbrn,"!is.na(mbrn$lopenr_far)")

dad_agemax<-mbrn_dad[,c("lopenr_far","fars_alder")]
dad_agemax<-dad_agemax[order(dad_agemax$lopenr_far,-abs(dad_agemax$fars_alder)),]
dad_agemax<-dad_agemax[!duplicated(dad_agemax$lopenr_far),]

dad_yearpregmax<-mbrn_dad[,c("lopenr_far","faar")]
dad_yearpregmax<-dad_yearpregmax[order(dad_yearpregmax$lopenr_far,abs(dad_yearpregmax$faar)),]
dad_yearpregmax<-dad_yearpregmax[!duplicated(dad_yearpregmax$lopenr_far),]

dad_birthyear<-mbrn_dad[,c("lopenr_far","birthyear_dad")]
dad_birthyear<-dad_birthyear[order(dad_birthyear$lopenr_far,abs(dad_birthyear$birthyear_dad)),]
dad_birthyear<-dad_birthyear[!duplicated(dad_birthyear$lopenr_far),]

dad_paritymax<-mbrn_dad[,c("lopenr_far","paritet_5")]
dad_paritymax<-dad_paritymax[order(dad_paritymax$lopenr_far,-abs(dad_paritymax$paritet_5)),]
dad_paritymax<-dad_paritymax[!duplicated(dad_paritymax$lopenr_far),]

dad_eclampsia<-mbrn_dad[,c("lopenr_far","eclampsia")]
dad_eclampsia<-dad_eclampsia[order(dad_eclampsia$lopenr_far,-abs(dad_eclampsia$eclampsia)),]
dad_eclampsia<-dad_eclampsia[!duplicated(dad_eclampsia$lopenr_far),]

dad_hta_preg<-mbrn_dad[,c("lopenr_far","hta_preg")]
dad_hta_preg<-dad_hta_preg[order(dad_hta_preg$lopenr_far,-abs(dad_hta_preg$hta_preg)),]
dad_hta_preg<-dad_hta_preg[!duplicated(dad_hta_preg$lopenr_far),]

dad_preterm<-mbrn_dad[,c("lopenr_far","preterm")]
dad_preterm<-dad_preterm[order(dad_preterm$lopenr_far,-abs(dad_preterm$preterm)),]
dad_preterm<-dad_preterm[!duplicated(dad_preterm$lopenr_far),]

dad_preterm_sp<-mbrn_dad[,c("lopenr_far","preterm_sp")]
dad_preterm_sp<-dad_preterm_sp[order(dad_preterm_sp$lopenr_far,-abs(dad_preterm_sp$preterm_sp)),]
dad_preterm_sp<-dad_preterm_sp[!duplicated(dad_preterm_sp$lopenr_far),]

dad_gdm<-mbrn_dad[,c("lopenr_far","gdm")]
dad_gdm<-dad_gdm[order(dad_gdm$lopenr_far,-abs(dad_gdm$gdm)),]
dad_gdm<-dad_gdm[!duplicated(dad_gdm$lopenr_far),]

dad_sga<-mbrn_dad[,c("lopenr_far","sga")]
dad_sga<-dad_sga[order(dad_sga$lopenr_far,-abs(dad_sga$sga)),]
dad_sga<-dad_sga[!duplicated(dad_sga$lopenr_far),]

dad_lga<-mbrn_dad[,c("lopenr_far","lga")]
dad_lga<-dad_lga[order(dad_lga$lopenr_far,-abs(dad_lga$lga)),]
dad_lga<-dad_lga[!duplicated(dad_lga$lopenr_far),]

dad_miscarriage12_mbrn<-mbrn_dad[,c("lopenr_far","miscarriage12_mbrn")]
dad_miscarriage12_mbrn<-dad_miscarriage12_mbrn[order(dad_miscarriage12_mbrn$lopenr_far,-abs(dad_miscarriage12_mbrn$miscarriage12_mbrn)),]
dad_miscarriage12_mbrn<-dad_miscarriage12_mbrn[!duplicated(dad_miscarriage12_mbrn$lopenr_far),]

dad_miscarriage24_mbrn<-mbrn_dad[,c("lopenr_far","miscarriage24_mbrn")]
dad_miscarriage24_mbrn<-dad_miscarriage24_mbrn[order(dad_miscarriage24_mbrn$lopenr_far,-abs(dad_miscarriage24_mbrn$miscarriage24_mbrn)),]
dad_miscarriage24_mbrn<-dad_miscarriage24_mbrn[!duplicated(dad_miscarriage24_mbrn$lopenr_far),]

dad_miscarriage_mbrn<-mbrn_dad[,c("lopenr_far","miscarriage_mbrn")]
dad_miscarriage_mbrn<-dad_miscarriage_mbrn[order(dad_miscarriage_mbrn$lopenr_far,-abs(dad_miscarriage_mbrn$miscarriage_mbrn)),]
dad_miscarriage_mbrn<-dad_miscarriage_mbrn[!duplicated(dad_miscarriage_mbrn$lopenr_far),]

dad_stillbirth24_mbrn<-mbrn_dad[,c("lopenr_far","stillbirth24_mbrn")]
dad_stillbirth24_mbrn<-dad_stillbirth24_mbrn[order(dad_stillbirth24_mbrn$lopenr_far,-abs(dad_stillbirth24_mbrn$stillbirth24_mbrn)),]
dad_stillbirth24_mbrn<-dad_stillbirth24_mbrn[!duplicated(dad_stillbirth24_mbrn$lopenr_far),]

dad_stillbirth<-mbrn_dad[,c("lopenr_far","stillbirth")]
dad_stillbirth<-dad_stillbirth[order(dad_stillbirth$lopenr_far,-abs(dad_stillbirth$stillbirth)),]
dad_stillbirth<-dad_stillbirth[!duplicated(dad_stillbirth$lopenr_far),]

mbrn_dad<-mbrn_dad[,c("mobaid_mfr_2374","lopenr_far","flerfodsel")]
mbrn_dad<-merge2(mbrn_dad,dad_eclampsia,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_hta_preg,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_preterm,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_preterm_sp,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_gdm,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_sga,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_lga,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_miscarriage12_mbrn,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_miscarriage24_mbrn,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_miscarriage_mbrn,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_stillbirth24_mbrn,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_stillbirth,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_agemax,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_yearpregmax,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_birthyear,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_paritymax,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-subset2(mbrn_dad,"!is.na(mbrn_dad$mobaid_mfr_2374)")
names(mbrn_dad)<-c("lopenr_far","mobaid_mfr_2374","flerfodsel","eclampsia_dad","hta_preg_dad",
                   "preterm_dad","preterm_sp_dad","gdm_dad","sga_dad","lga_dad",
                   "miscarriage12_mbrn_dad","miscarriage24_mbrn_dad","miscarriage_mbrn_dad",
                   "stillbirth24_mbrn_dad","stillbirth_dad","agemax_dad","yearpregmax_dad","birthyear_dad","paritymax_dad")
mbrn_dad$agemax_dad<-mbrn_dad$yearpregmax_dad-mbrn_dad$birthyear_dad


mbrn_key<-spss.get("N:/durable/RAW/1.MoBa_questionnaires_MBRN/2.MBRN_Birth_Registry/If_you_want_records_MoBa_parents_siblings_children/Koblingsbro_MFR_MoBa.sav",use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(mbrn_key)<-tolower(names(mbrn_key))
mbrn_key<-mbrn_key[,c("preg_id_2374","mobaid_mfr_2374")]

mbrn_mom<-merge2(mbrn_mom,mbrn_key,by.id=c("mobaid_mfr_2374"),all.x=TRUE,sort=FALSE)
mbrn_mom<-subset2(mbrn_mom,"mbrn_mom$flerfodsel==0")
mbrn_mom$flerfodsel<-NULL
mbrn_mom$mobaid_mfr_2374<-NULL

mbrn_dad<-merge2(mbrn_dad,mbrn_key,by.id=c("mobaid_mfr_2374"),all.x=TRUE,sort=FALSE)
mbrn_dad<-subset2(mbrn_dad,"mbrn_dad$flerfodsel==0")
mbrn_dad$flerfodsel<-NULL
mbrn_dad$mobaid_mfr_2374<-NULL


### FINAL FORMAT OF DATA ###

dat<-rename.vars(dat,
                 from=c("mors_alder","fars_alder"),
                 to=c("agedelivery_mom","agedelivery_dad"))

# dim(dat)[1] = 114629
dat$exclude<-with(dat,ifelse(flerfodsel==1,1,0))
dat$exclude<-with(dat,ifelse(is.na(exclude),0,exclude))
table(dat$exclude)
dat<-subset2(dat,"dat$exclude==0")

# dim(dat)[1] = 110663
dat$exclude<-with(dat,ifelse(versjon_skjema1_tbl1=="",1,0))
dat$exclude<-with(dat,ifelse(is.na(exclude),0,exclude))
table(dat$exclude)
dat<-subset2(dat,"dat$exclude==0")

# dim(dat)[1] = 100315

dat<-dat[,c("sentrixid_mom","sentrixid_dad","m_id_2374","f_id_2374","preg_id_2374","batch_mom","batch_dad",
            "agedelivery_mom","agedelivery_dad","bmi_mom","bmi_dad","eduyears_mom","eduyears_dad","smkinit_mom","smkinit_dad",
            "subf","subf_12plus","miscarriage","cad_grs_mom","cad_grs_dad",
            "genotype_batch_mom","pc01_mom","pc02_mom","pc03_mom","pc04_mom","pc05_mom","pc06_mom","pc07_mom","pc08_mom","pc09_mom","pc10_mom",
            "pc11_mom","pc12_mom","pc13_mom","pc14_mom","pc15_mom","pc16_mom","pc17_mom","pc18_mom","pc19_mom","pc20_mom",
            "genotype_batch_dad","pc01_dad","pc02_dad","pc03_dad","pc04_dad","pc05_dad","pc06_dad","pc07_dad","pc08_dad","pc09_dad","pc10_dad",
            "pc11_dad","pc12_dad","pc13_dad","pc14_dad","pc15_dad","pc16_dad","pc17_dad","pc18_dad","pc19_dad","pc20_dad")]

dat<-merge2(dat,mbrn_mom,by.id=c("preg_id_2374"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,mbrn_dad,by.id=c("preg_id_2374"),all.x=TRUE,sort=FALSE)


### EXCLUSIONS - REMOVAL OF CONSENT ###

excl<-spss.get("N:/durable/RAW/1.MoBa_questionnaires_MBRN/INFO_files_start_your_analysis_here/PDB2374_SV_INFO_V12_20220824.sav",
               use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(excl)<-tolower(names(excl))
excl$consent<-1
excl<-excl[,c("preg_id_2374","consent")]
dat<-merge2(dat,excl,by.id=c("preg_id_2374"),all.x=TRUE,sort=FALSE)
dat$consent<-with(dat,ifelse(is.na(consent),0,consent))
dat<-dat[dat$consent==1,]

# dim(dat)[1] = 99973


# Unique variables in MoBa #

mom_subf<-dat[,c("m_id_2374","subf_12plus")]
mom_subf<-mom_subf[order(mom_subf$m_id_2374,-abs(mom_subf$subf_12plus)),]
mom_subf<-mom_subf[!duplicated(mom_subf$m_id_2374),]
mom_subf<-rename.vars(mom_subf,from=c("subf_12plus"),to=c("subf_mom"))

mom_misc<-dat[,c("m_id_2374","miscarriage")]
mom_misc<-mom_misc[order(mom_misc$m_id_2374,-abs(mom_misc$miscarriage)),]
mom_misc<-mom_misc[!duplicated(mom_misc$m_id_2374),]
mom_misc<-rename.vars(mom_misc,from=c("miscarriage"),to=c("miscarriage_moba_mom"))

mom_bmimax<-dat[,c("m_id_2374","bmi_mom")]
mom_bmimax<-mom_bmimax[order(mom_bmimax$m_id_2374,-abs(mom_bmimax$bmi_mom)),]
mom_bmimax<-mom_bmimax[!duplicated(mom_bmimax$m_id_2374),]
mom_bmimax<-rename.vars(mom_bmimax,from=c("bmi_mom"),to=c("bmimax_mom"))
summary(mom_bmimax$bmimax_mom)

mom_eduyearsmax<-dat[,c("m_id_2374","eduyears_mom")]
mom_eduyearsmax<-mom_eduyearsmax[order(mom_eduyearsmax$m_id_2374,-abs(mom_eduyearsmax$eduyears_mom)),]
mom_eduyearsmax<-mom_eduyearsmax[!duplicated(mom_eduyearsmax$m_id_2374),]
mom_eduyearsmax<-rename.vars(mom_eduyearsmax,from=c("eduyears_mom"),to=c("eduyearsmax_mom"))

mom_smkinitmax<-dat[,c("m_id_2374","smkinit_mom")]
mom_smkinitmax<-mom_smkinitmax[order(mom_smkinitmax$m_id_2374,-abs(mom_smkinitmax$smkinit_mom)),]
mom_smkinitmax<-mom_smkinitmax[!duplicated(mom_smkinitmax$m_id_2374),]
mom_smkinitmax<-rename.vars(mom_smkinitmax,from=c("smkinit_mom"),to=c("smkinitmax_mom"))

moba_mom<-merge2(mom_subf,mom_misc,by.id=c("m_id_2374"),all.x=TRUE,sort=FALSE)
moba_mom<-merge2(moba_mom,mom_bmimax,by.id=c("m_id_2374"),all.x=TRUE,sort=FALSE)
moba_mom<-merge2(moba_mom,mom_eduyearsmax,by.id=c("m_id_2374"),all.x=TRUE,sort=FALSE)
moba_mom<-merge2(moba_mom,mom_smkinitmax,by.id=c("m_id_2374"),all.x=TRUE,sort=FALSE)


dad_subf<-dat[,c("f_id_2374","subf_12plus")]
dad_subf<-dad_subf[order(dad_subf$f_id_2374,-abs(dad_subf$subf_12plus)),]
dad_subf<-dad_subf[!duplicated(dad_subf$f_id_2374),]
dad_subf<-rename.vars(dad_subf,from=c("subf_12plus"),to=c("subf_dad"))
dad_subf<-dad_subf[dad_subf$f_id_2374!="", ]

dad_misc<-dat[,c("f_id_2374","miscarriage")]
dad_misc<-dad_misc[order(dad_misc$f_id_2374,-abs(dad_misc$miscarriage)),]
dad_misc<-dad_misc[!duplicated(dad_misc$f_id_2374),]
dad_misc<-rename.vars(dad_misc,from=c("miscarriage"),to=c("miscarriage_moba_dad"))
dad_misc<-dad_misc[dad_misc$f_id_2374!="", ]

dad_bmimax<-dat[,c("f_id_2374","bmi_dad")]
dad_bmimax<-dad_bmimax[order(dad_bmimax$f_id_2374,-abs(dad_bmimax$bmi_dad)),]
dad_bmimax<-dad_bmimax[!duplicated(dad_bmimax$f_id_2374),]
dad_bmimax<-rename.vars(dad_bmimax,from=c("bmi_dad"),to=c("bmimax_dad"))
dad_bmimax<-dad_bmimax[dad_bmimax$f_id_2374!="", ]

dad_eduyearsmax<-dat[,c("f_id_2374","eduyears_dad")]
dad_eduyearsmax<-dad_eduyearsmax[order(dad_eduyearsmax$f_id_2374,-abs(dad_eduyearsmax$eduyears_dad)),]
dad_eduyearsmax<-dad_eduyearsmax[!duplicated(dad_eduyearsmax$f_id_2374),]
dad_eduyearsmax<-rename.vars(dad_eduyearsmax,from=c("eduyears_dad"),to=c("eduyearsmax_dad"))
dad_eduyearsmax<-dad_eduyearsmax[dad_eduyearsmax$f_id_2374!="", ]

dad_smkinitmax<-dat[,c("f_id_2374","smkinit_dad")]
dad_smkinitmax<-dad_smkinitmax[order(dad_smkinitmax$f_id_2374,-abs(dad_smkinitmax$smkinit_dad)),]
dad_smkinitmax<-dad_smkinitmax[!duplicated(dad_smkinitmax$f_id_2374),]
dad_smkinitmax<-rename.vars(dad_smkinitmax,from=c("smkinit_dad"),to=c("smkinitmax_dad"))
dad_smkinitmax<-dad_smkinitmax[dad_smkinitmax$f_id_2374!="", ]

moba_dad<-merge2(dad_subf,dad_misc,by.id=c("f_id_2374"),all.x=TRUE,sort=FALSE)
moba_dad<-merge2(moba_dad,dad_bmimax,by.id=c("f_id_2374"),all.x=TRUE,sort=FALSE)
moba_dad<-merge2(moba_dad,dad_eduyearsmax,by.id=c("f_id_2374"),all.x=TRUE,sort=FALSE)
moba_dad<-merge2(moba_dad,dad_smkinitmax,by.id=c("f_id_2374"),all.x=TRUE,sort=FALSE)


dat<-merge2(dat,moba_mom,by.id=c("m_id_2374"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,moba_dad,by.id=c("f_id_2374"),all.x=TRUE,sort=FALSE)


### DEFINITION OF STUDY FLOW CHART ###

# Total singleton pregnancies
length(which(!is.na(dat$preg_id_2374)))

# Unique IDs of mothers/fathers
length(unique(dat[which(!is.na(dat$subf_12plus==0)),c("m_id_2374")]))
length(unique(dat[which(!is.na(dat$subf_12plus==0)),c("f_id_2374")]))

# IDs of mothers/fathers with genotype data
length(dat[which(!is.na(dat$cad_grs_mom)),c("sentrixid_mom")])
length(dat[which(!is.na(dat$cad_grs_dad)),c("sentrixid_dad")])

# IDs of unique mothers/fathers with genotype data
length(unique(dat[which(!is.na(dat$cad_grs_mom)),c("sentrixid_mom")]))
length(unique(dat[which(!is.na(dat$cad_grs_dad)),c("sentrixid_dad")]))


datx<-subset2(dat,"!is.na(dat$cad_grs_mom)")
jpeg(filename="./Outputs/descriptive/cad_grs_mom.jpg",
     width=9000,height=9000,res=1200,pointsize=11.5)
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~cad_grs_mom,data=datx))
dev.off()

datx<-subset2(dat,"!is.na(dat$cad_grs_dad)")
jpeg(filename="./Outputs/descriptive/cad_grs_dad.jpg",
     width=9000,height=9000,res=1200,pointsize=11.5)
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~cad_grs_dad,data=datx))
dev.off()

attributes(dat$bmi_mom)$label<-c("bmi_mom")
attributes(dat$bmi_dad)$label<-c("bmi_dad")

save(dat,file="N:/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")


#####################
### MAIN ANALYSES ###
#####################

setwd("N:/durable/Projects/Hernaez_MR_CVD")

dir.create("./Outputs")
dir.create("./Outputs/descriptive")
dir.create("./Outputs/results")

setwd("N:/durable/Projects/Hernaez_MR_CVD/Outputs")


### POPULATION DESCRIPTION ###
##############################

load("N:/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")

datx<-subset2(dat,"!is.na(dat$cad_grs_mom)")
datx<-datx[!duplicated(datx$m_id_2374),]

xxx<-datx[,c("birthyear_mom","agemax_mom","eduyearsmax_mom","bmimax_mom","smkinitmax_mom","paritymax_mom",
             "hta_preg_mom","eclampsia_mom","gdm_mom","preterm_mom","preterm_sp_mom","sga_mom","miscarriage_mbrn_mom","stillbirth_mom",
             "subf_mom","lga_mom")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx, method=c("birthyear_mom"=2,"agemax_mom"=2,"bmimax_mom"=2,"smkinitmax_mom"=3,"paritymax_mom"=3,
                                             "subf_mom"=3,"miscarriage_mbrn_mom"=3,"stillbirth_mom"=3,"eclampsia_mom"=3,
                                             "hta_preg_mom"=3,"preterm_mom"=3,"preterm_sp_mom"=3,"sga_mom"=3,"lga_mom"=3,"gdm_mom"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(cbind(all$descr[,1]))
colnames(tab1)<-c("Mothers-All")
write.table(tab1,file="./descriptive/descriptive_mothers.csv",sep=";",col.names=NA)


datx<-subset2(dat,"!is.na(dat$cad_grs_dad)")
datx<-datx[!duplicated(datx$f_id_2374),]

xxx<-datx[,c("birthyear_dad","agemax_dad","eduyearsmax_dad","bmimax_dad","smkinitmax_dad","paritymax_dad",
             "hta_preg_dad","eclampsia_dad","gdm_dad","preterm_dad","preterm_sp_dad","sga_dad","miscarriage_mbrn_dad","stillbirth_dad",
             "subf_dad","lga_dad")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx, method=c("birthyear_dad"=2,"agemax_dad"=2,"bmimax_dad"=2,"smkinitmax_dad"=3,"paritymax_dad"=3,
                                             "subf_dad"=3,"miscarriage_mbrn_dad"=3,"stillbirth_dad"=3,"eclampsia_dad"=3,
                                             "hta_preg_dad"=3,"preterm_dad"=3,"preterm_sp_dad"=3,"sga_dad"=3,"lga_dad"=3,"gdm_dad"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(cbind(all$descr[,1]))
colnames(tab1)<-c("Fathers-All")
write.table(tab1,file="./descriptive/descriptive_fathers.csv",sep=";",col.names=NA)


### SELECTION BIAS - INCLUDED vs NON-INCLUDED ###

datx<-dat
datx$exclusion<-with(datx,ifelse(!is.na(datx$cad_grs_mom),0,1))
datx<-datx[!duplicated(datx$m_id_2374),]
datx$paritymax_mom<-datx$paritymax_mom+1

xxx<-datx[,c("exclusion","birthyear_mom","agemax_mom","eduyearsmax_mom","bmimax_mom","smkinitmax_mom","paritymax_mom",
             "hta_preg_mom","eclampsia_mom","gdm_mom","preterm_mom","preterm_sp_mom","sga_mom","miscarriage_mbrn_mom","stillbirth_mom",
             "subf_mom","lga_mom")]

all<-NULL
all<-createTable(compareGroups(exclusion~.,
                               xxx, method=c("birthyear_mom"=2,"agemax_mom"=2,"bmimax_mom"=2,"smkinitmax_mom"=3,"paritymax_mom"=3,
                                             "subf_mom"=3,"miscarriage_mbrn_mom"=3,"stillbirth_mom"=3,"eclampsia_mom"=3,
                                             "hta_preg_mom"=3,"preterm_mom"=3,"preterm_sp_mom"=3,"sga_mom"=3,"lga_mom"=3,"gdm_mom"=3)),
                 show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Included","Non-included","P-value","N")
write.table(tab1,file="./descriptive/selectionbias_mothers.csv",sep=";",col.names=NA)


datx<-dat
datx$exclusion<-with(datx,ifelse(!is.na(datx$cad_grs_dad),0,1))
datx<-datx[!duplicated(datx$f_id_2374),]
datx$paritymax_dad<-datx$paritymax_dad+1

xxx<-datx[,c("exclusion","birthyear_dad","agemax_dad","eduyearsmax_dad","bmimax_dad","smkinitmax_dad","paritymax_dad",
             "hta_preg_dad","eclampsia_dad","gdm_dad","preterm_dad","preterm_sp_dad","sga_dad","miscarriage_mbrn_dad","stillbirth_dad",
             "subf_dad","lga_dad")]

all<-NULL
all<-createTable(compareGroups(exclusion~.,
                               xxx, method=c("birthyear_dad"=2,"agemax_dad"=2,"bmimax_dad"=2,"smkinitmax_dad"=3,"paritymax_dad"=3,
                                             "subf_dad"=3,"miscarriage_mbrn_dad"=3,"stillbirth_dad"=3,"eclampsia_dad"=3,
                                             "hta_preg_dad"=3,"preterm_dad"=3,"preterm_sp_dad"=3,"sga_dad"=3,"lga_dad"=3,"gdm_dad"=3)),
                 show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Included","Non-included","P-value","N")
write.table(tab1,file="./descriptive/selectionbias_fathers.csv",sep=";",col.names=NA)


### LOGISTIC REGRESSION / MENDELIAN RANDOMIZATION: LINEAR ASSOCIATIONS ###
##########################################################################

load("N:/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")
dat$adj<-0

vars00<-c("cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad",
          "cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad",
          "cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad",
          "cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad",
          "cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad")
vars01<-c("cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z",
          "cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z",
          "cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z",
          "cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z",
          "cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z")
vars02<-c("hta_preg_mom","hta_preg_dad","eclampsia_mom","eclampsia_dad","gdm_mom","gdm_dad",
          "preterm_mom","preterm_dad","preterm_sp_mom","preterm_sp_dad","sga_mom","sga_dad",
          "miscarriage_mbrn_mom","miscarriage_mbrn_dad","stillbirth_mom","stillbirth_dad",
          "subf_mom","subf_dad","lga_mom","lga_dad")
vars03<-c("hta_chronic_mom","hta_chronic_mom","hta_chronic_mom","hta_chronic_mom",
          "adj","adj","adj","adj","adj","adj","adj","adj","adj","adj","adj","adj","adj","adj","adj","adj")
vars04<-c("cad_grs_dad","cad_grs_mom","cad_grs_dad","cad_grs_mom",
          "cad_grs_dad","cad_grs_mom","cad_grs_dad","cad_grs_mom",
          "cad_grs_dad","cad_grs_mom","cad_grs_dad","cad_grs_mom",
          "cad_grs_dad","cad_grs_mom","cad_grs_dad","cad_grs_mom",
          "cad_grs_dad","cad_grs_mom","cad_grs_dad","cad_grs_mom")
vars05<-c("cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z",
          "cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z",
          "cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z",
          "cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z",
          "cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z")
vars08<-c("pc01_mom","pc01_dad","pc01_mom","pc01_dad","pc01_mom","pc01_dad","pc01_mom","pc01_dad","pc01_mom","pc01_dad",
          "pc01_mom","pc01_dad","pc01_mom","pc01_dad","pc01_mom","pc01_dad","pc01_mom","pc01_dad","pc01_mom","pc01_dad")
vars09<-c("pc02_mom","pc02_dad","pc02_mom","pc02_dad","pc02_mom","pc02_dad","pc02_mom","pc02_dad","pc02_mom","pc02_dad",
          "pc02_mom","pc02_dad","pc02_mom","pc02_dad","pc02_mom","pc02_dad","pc02_mom","pc02_dad","pc02_mom","pc02_dad")
vars10<-c("pc03_mom","pc03_dad","pc03_mom","pc03_dad","pc03_mom","pc03_dad","pc03_mom","pc03_dad","pc03_mom","pc03_dad",
          "pc03_mom","pc03_dad","pc03_mom","pc03_dad","pc03_mom","pc03_dad","pc03_mom","pc03_dad","pc03_mom","pc03_dad")
vars11<-c("pc04_mom","pc04_dad","pc04_mom","pc04_dad","pc04_mom","pc04_dad","pc04_mom","pc04_dad","pc04_mom","pc04_dad",
          "pc04_mom","pc04_dad","pc04_mom","pc04_dad","pc04_mom","pc04_dad","pc04_mom","pc04_dad","pc04_mom","pc04_dad")
vars12<-c("pc05_mom","pc05_dad","pc05_mom","pc05_dad","pc05_mom","pc05_dad","pc05_mom","pc05_dad","pc05_mom","pc05_dad",
          "pc05_mom","pc05_dad","pc05_mom","pc05_dad","pc05_mom","pc05_dad","pc05_mom","pc05_dad","pc05_mom","pc05_dad")
vars13<-c("pc06_mom","pc06_dad","pc06_mom","pc06_dad","pc06_mom","pc06_dad","pc06_mom","pc06_dad","pc06_mom","pc06_dad",
          "pc06_mom","pc06_dad","pc06_mom","pc06_dad","pc06_mom","pc06_dad","pc06_mom","pc06_dad","pc06_mom","pc06_dad")
vars14<-c("pc07_mom","pc07_dad","pc07_mom","pc07_dad","pc07_mom","pc07_dad","pc07_mom","pc07_dad","pc07_mom","pc07_dad",
          "pc07_mom","pc07_dad","pc07_mom","pc07_dad","pc07_mom","pc07_dad","pc07_mom","pc07_dad","pc07_mom","pc07_dad")
vars15<-c("pc08_mom","pc08_dad","pc08_mom","pc08_dad","pc08_mom","pc08_dad","pc08_mom","pc08_dad","pc08_mom","pc08_dad",
          "pc08_mom","pc08_dad","pc08_mom","pc08_dad","pc08_mom","pc08_dad","pc08_mom","pc08_dad","pc08_mom","pc08_dad")
vars16<-c("pc09_mom","pc09_dad","pc09_mom","pc09_dad","pc09_mom","pc09_dad","pc09_mom","pc09_dad","pc09_mom","pc09_dad",
          "pc09_mom","pc09_dad","pc09_mom","pc09_dad","pc09_mom","pc09_dad","pc09_mom","pc09_dad","pc09_mom","pc09_dad")
vars17<-c("pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad",
          "pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad")
vars18<-c("genotype_batch_mom","genotype_batch_dad","genotype_batch_mom","genotype_batch_dad","genotype_batch_mom","genotype_batch_dad",
          "genotype_batch_mom","genotype_batch_dad","genotype_batch_mom","genotype_batch_dad","genotype_batch_mom","genotype_batch_dad",
          "genotype_batch_mom","genotype_batch_dad","genotype_batch_mom","genotype_batch_dad","genotype_batch_mom","genotype_batch_dad",
          "genotype_batch_mom","genotype_batch_dad")
vars19<-c("m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374",
          "m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374")
vars20<-c("CAD-GRS - HT in pregnancy, mothers","CAD-GRS - HT in pregnancy, fathers",
          "CAD-GRS - Pre-/Eclampsia, mothers","CAD-GRS - Pre-/Eclampsia, fathers",
          "CAD-GRS - Gestat. diabetes, mothers","CAD-GRS - Gestat. diabetes, fathers",
          "CAD-GRS - Preterm, mothers","CAD-GRS - Preterm, fathers",
          "CAD-GRS - Spontaneous preterm, mothers","CAD-GRS - Spontaneous preterm, fathers",
          "CAD-GRS - SGA, mothers","CAD-GRS - SGA, fathers",
          "CAD-GRS - Miscarriage, mothers","CAD-GRS - Miscarriage, fathers",
          "CAD-GRS - Stillbirth, mothers","CAD-GRS - Stillbirth, fathers",
          "CAD-GRS - Subfertility, mothers","CAD-GRS - Subfertility, fathers",
          "CAD-GRS - LGA, mothers","CAD-GRS - LGA, fathers")

tab<-NULL
z<-qnorm(1-0.05/2)
for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"!is.na(dat[,vars00[i]]) & !is.na(dat[,vars02[i]]) & dat[,vars03[i]]==0 & !is.na(dat[,vars18[i]])")
  datx<-datx[!duplicated(datx[,vars19[i]]),]
  datx[,vars01[i]]<-as.numeric(with(datx,scale(datx[,vars00[i]])))
  sample<-dim(datx)[1]
  
  mod01<-glm(as.factor(datx[,vars02[i]])~datx[,vars01[i]],
             data=datx, family="binomial")
  estimate01<-as.numeric(summary(mod01)$coefficients[2,1])
  se01<-as.numeric(summary(mod01)$coefficients[2,2])
  coef01<-risk_se_ic_guapa(estimate01,se01)
  pval01<-pval_guapa(as.numeric(summary(mod01)$coefficients[2,4]))
  
  mod02<-glm(as.factor(datx[,vars02[i]])~datx[,vars01[i]]
             +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
             +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
             +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]]+as.factor(datx[,vars18[i]]),
             data=datx, family="binomial")
  estimate02<-as.numeric(summary(mod02)$coefficients[2,1])
  se02<-as.numeric(summary(mod02)$coefficients[2,2])
  coef02<-risk_se_ic_guapa(estimate02,se02)
  pval02<-pval_guapa(as.numeric(summary(mod02)$coefficients[2,4]))
  or02<-round(exp(estimate02),5)
  ic95lo02<-round(exp(estimate02-(z*se02)),5)
  ic95hi02<-round(exp(estimate02+(z*se02)),5)

  dat2<-subset2(datx,"!is.na(datx[,vars04[i]])")
  dat2[,vars05[i]]<-as.numeric(with(dat2,scale(dat2[,vars04[i]])))
  sample2<-dim(dat2)[1]
  
  mod03<-glm(as.factor(dat2[,vars02[i]])~dat2[,vars01[i]]
             +dat2[,vars08[i]]+dat2[,vars09[i]]+dat2[,vars10[i]]+dat2[,vars11[i]]
             +dat2[,vars12[i]]+dat2[,vars13[i]]+dat2[,vars14[i]]
             +dat2[,vars15[i]]+dat2[,vars16[i]]+dat2[,vars17[i]]+as.factor(dat2[,vars18[i]])+dat2[,vars05[i]],
             data=dat2, family="binomial")
  estimate03<-as.numeric(summary(mod03)$coefficients[2,1])
  se03<-as.numeric(summary(mod03)$coefficients[2,2])
  coef03<-risk_se_ic_guapa(estimate03,se03)
  pval03<-pval_guapa(as.numeric(summary(mod03)$coefficients[2,4]))
  or03<-round(exp(estimate03),5)
  ic95lo03<-round(exp(estimate03-(z*se03)),5)
  ic95hi03<-round(exp(estimate03+(z*se03)),5)
  
  aaa<-datx[,vars01[i]]
  mod_base<-glm(formula=as.factor(datx[,vars02[i]])~datx[,vars08[i]]
                +datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
                +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
                +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]]+as.factor(datx[,vars18[i]]),
                data=datx, family="binomial")
  mod_lin<-glm(formula=as.factor(datx[,vars02[i]])~aaa
               +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
               +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
               +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]]+as.factor(datx[,vars18[i]]),
               data=datx, family="binomial")
  mod_nlin<-glm(formula=as.factor(datx[,vars02[i]])~bs(aaa,df=4)
                +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
                +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
                +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]]+as.factor(datx[,vars18[i]]),
                data=datx, family="binomial")
  p_lrtest02<-pval_guapa(lrtest(mod_lin,mod_nlin)[2,5])
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02,p_lrtest02,estimate02,se02,or02,ic95lo02,ic95hi02,sample,
                       coef03,pval03,estimate03,se03,or03,ic95lo03,ic95hi03,sample2))
}

colnames(tab)<-c("coef (raw)","pval (raw)","coef (adj.)","pval (adj.)","p_lrtest",
                 "beta","se","or","ic95lo","ic95hi","sample",
                 "coef (partner adj.)","pval (partner adj.)","beta (partner adj.)","se (partner adj.)",
                 "or (partner adj.)","ic95lo (partner adj.)","ic95hi (partner adj.)","sample (partner adj.)")
rownames(tab)<-vars20
write.table(tab,file="./results/linear.csv",sep=";",col.names=NA)


### HORIZONTAL PLEIOTROPY: ASSOCIATION OF GRS WITH COVARIATES ###
#################################################################

load("N:/data/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")

vars00<-c("cad_grs_mom","cad_grs_dad")
vars01<-c("cad_grs_mom_z","cad_grs_dad_z")
vars04<-c("birthyear_mom","birthyear_dad")
vars05<-c("eduyearsmax_mom","eduyearsmax_dad")
vars06<-c("bmimax_mom","bmimax_dad")
vars07<-c("smkinitmax_mom","smkinitmax_dad")
vars08<-c("paritymax_mom","paritymax_dad")
vars09<-c("m_id_2374","f_id_2374")

tab<-NULL
tab2<-NULL

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"!is.na(dat[,vars00[i]])")
  datx<-datx[!duplicated(datx[,vars09[i]]),]
  
  datx[,vars01[i]]<-as.numeric(with(datx,scale(datx[,vars00[i]])))
  mod01<-lm_robust(datx[,vars04[i]]~datx[,vars01[i]], 
                   data=datx, clusters=datx[,vars09[i]], se_type="stata")
  mod02<-lm_robust(datx[,vars05[i]]~datx[,vars01[i]], 
                   data=datx, clusters=datx[,vars09[i]], se_type="stata")
  mod03<-lm_robust(datx[,vars06[i]]~datx[,vars01[i]], 
                   data=datx, clusters=datx[,vars09[i]], se_type="stata")
  mod04<-lm_robust(datx[,vars08[i]]~datx[,vars01[i]], 
                   data=datx, clusters=datx[,vars09[i]], se_type="stata")
  mod05<-miceadds::glm.cluster(formula=as.factor(datx[,vars07[i]])~datx[,vars01[i]],
                               data=datx, cluster=datx[,vars09[i]], family="binomial")
  
  coef01<-paste(guapa(summary(mod01)$coefficients[2,1])," [",
                guapa(summary(mod01)$coefficients[2,5]),"; ",
                guapa(summary(mod01)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod01)$coefficients[2,4]),")",sep="")
  coef02<-paste(guapa(summary(mod02)$coefficients[2,1])," [",
                guapa(summary(mod02)$coefficients[2,5]),"; ",
                guapa(summary(mod02)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod02)$coefficients[2,4]),")",sep="")
  coef03<-paste(guapa(summary(mod03)$coefficients[2,1])," [",
                guapa(summary(mod03)$coefficients[2,5]),"; ",
                guapa(summary(mod03)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod03)$coefficients[2,4]),")",sep="")
  coef04<-paste(guapa(summary(mod04)$coefficients[2,1])," [",
                guapa(summary(mod04)$coefficients[2,5]),"; ",
                guapa(summary(mod04)$coefficients[2,6]),"] (P=",
                pval_guapa(summary(mod04)$coefficients[2,4]),")",sep="")
  coef05<-paste(risk_se_ic_guapa(as.numeric(summary(mod05)[2,1]),as.numeric(summary(mod05)[2,2]))," (P=",
                pval_guapa(as.numeric(summary(mod05)[2,4])),")",sep="")
  
  beta01<-round(summary(mod01)$coefficients[2,1],6)
  se01<-round(summary(mod01)$coefficients[2,2],6)
  beta02<-round(summary(mod02)$coefficients[2,1],6)
  se02<-round(summary(mod02)$coefficients[2,2],6)
  beta03<-round(summary(mod03)$coefficients[2,1],6)
  se03<-round(summary(mod03)$coefficients[2,2],6)
  beta04<-round(summary(mod04)$coefficients[2,1],6)
  se04<-round(summary(mod04)$coefficients[2,2],6)
  beta05<-round(summary(mod05)[2,1],6)
  se05<-round(summary(mod05)[2,2],6)
  
  tab<-rbind(tab,cbind(coef01,coef02,coef03,coef04,coef05))
  tab2<-rbind(tab2,cbind(beta01,se01,beta02,se02,beta03,se03,beta04,se04,beta05,se05))
}

rownames(tab)<-vars01
colnames(tab)<-c("age","eduyears","bmi","parity","smk_init")
rownames(tab2)<-vars01
colnames(tab2)<-c("age_beta","age_se","eduyears_beta","eduyears_se","bmi_beta","bmi_se","parity_beta","parity_se","smk_beta","smk_se")
write.table(tab,file="./descriptive/hp_covars.csv",sep=";",col.names=NA)
write.table(tab2,file="./descriptive/hp_covars_metaanalysis.csv",sep=";",col.names=NA)


###################################
### GWAS FOR PREGNANCY OUTCOMES ###
###################################

dir.create("N:/data/durable/Projects/Hernaez_MR_CVD/gwas")
dir.create("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/preg_outc")
dir.create("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/preg_outc/results")

setwd("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/preg_outc")


###########################################
### GENERATION OF DATASETS FOR THE GWAS ###
###########################################

### Maternal and paternal IDs ###

dat<-as.data.frame(fread("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/cad/cad_grs.raw",header=TRUE,sep="\t",sep2="\t"))
xxx<-dat[,c("FID","IID")]
load("N:/data/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")
dat<-subset2(dat,"!is.na(dat$cad_grs_mom)")
dat<-as.data.frame(dat[,c("sentrixid_mom")])
names(dat)<-c("IID")
dat<-as.data.frame(dat[!duplicated(dat$IID),])
names(dat)<-c("IID")
dat<-merge2(dat,xxx,by.id=c("IID"),all.x=TRUE,sort=FALSE)
mom<-dat[,c("FID","IID")]
write.table(mom,"N:/data/durable/Projects/Hernaez_MR_CVD/gwas/preg_outc/maternal_ids_batch2022.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat<-as.data.frame(fread("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/cad/cad_grs.raw",header=TRUE,sep="\t",sep2="\t"))
xxx<-dat[,c("FID","IID")]
load("N:/data/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")
dat<-subset2(dat,"!is.na(dat$cad_grs_dad)")
dat<-as.data.frame(dat[,c("sentrixid_dad")])
names(dat)<-c("IID")
dat<-as.data.frame(dat[!duplicated(dat$IID),])
names(dat)<-c("IID")
dat<-merge2(dat,xxx,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dad<-dat[,c("FID","IID")]
write.table(dad,"N:/data/durable/Projects/Hernaez_MR_CVD/gwas/preg_outc/paternal_ids_batch2022.txt",sep="\t",row.names=FALSE, quote=FALSE)


### Maternal and paternal covariates ###

dat<-as.data.frame(fread("N:/data/durable/RAW/3.MoBa_genetics/p471/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt",header=TRUE,sep="\t",sep2="\t"))
momcv<-mom
dat<-dat[,c("IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
momcv<-merge2(momcv,dat,by.id=c("IID"),all.x=TRUE,sort=FALSE)
momcv<-momcv[,c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
write.table(momcv,"N:/data/durable/Projects/Hernaez_MR_CVD/gwas/preg_outc/maternal_cov.txt",sep="\t",row.names=FALSE, quote=FALSE)

dadcv<-dad
dadcv<-merge2(dadcv,dat,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dadcv<-dadcv[,c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
write.table(dadcv,"N:/data/durable/Projects/Hernaez_MR_CVD/gwas/preg_outc/paternal_cov.txt",sep="\t",row.names=FALSE, quote=FALSE)


### Phenotypes ###

load("N:/data/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")
dat<-subset2(dat,"!is.na(dat$cad_grs_mom)")

dat2<-dat[,c("sentrixid_mom","miscarriage_mbrn_mom")]
dat2<-rename.vars(dat2,from=c("sentrixid_mom","miscarriage_mbrn_mom"),to=c("IID","miscarriage_mbrn"))
dat2<-dat2[order(dat2$IID,-abs(dat2$miscarriage_mbrn)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,mom,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","miscarriage_mbrn")]
write.table(dat2,"./maternal_miscarriage_mbrn.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_mom","eclampsia_mom")]
dat2<-rename.vars(dat2,from=c("sentrixid_mom","eclampsia_mom"),to=c("IID","eclampsia"))
dat2<-dat2[order(dat2$IID,-abs(dat2$eclampsia)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,mom,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","eclampsia")]
write.table(dat2,"./maternal_eclampsia.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_mom","hta_preg_mom")]
dat2<-rename.vars(dat2,from=c("sentrixid_mom","hta_preg_mom"),to=c("IID","htapreg"))
dat2<-dat2[order(dat2$IID,-abs(dat2$htapreg)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,mom,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","htapreg")]
write.table(dat2,"./maternal_htapreg.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_mom","preterm_mom")]
dat2<-rename.vars(dat2,from=c("sentrixid_mom","preterm_mom"),to=c("IID","preterm"))
dat2<-dat2[order(dat2$IID,-abs(dat2$preterm)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,mom,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","preterm")]
write.table(dat2,"./maternal_preterm.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_mom","preterm_sp_mom")]
dat2<-rename.vars(dat2,from=c("sentrixid_mom","preterm_sp_mom"),to=c("IID","preterm_sp"))
dat2<-dat2[order(dat2$IID,-abs(dat2$preterm_sp)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,mom,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","preterm_sp")]
write.table(dat2,"./maternal_preterm_sp.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_mom","sga_mom")]
dat2<-rename.vars(dat2,from=c("sentrixid_mom","sga_mom"),to=c("IID","sga"))
dat2<-dat2[order(dat2$IID,-abs(dat2$sga)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,mom,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","sga")]
write.table(dat2,"./maternal_sga.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_mom","lga_mom")]
dat2<-rename.vars(dat2,from=c("sentrixid_mom","lga_mom"),to=c("IID","lga"))
dat2<-dat2[order(dat2$IID,-abs(dat2$lga)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,mom,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","lga")]
write.table(dat2,"./maternal_lga.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_mom","gdm_mom")]
dat2<-rename.vars(dat2,from=c("sentrixid_mom","gdm_mom"),to=c("IID","gdm"))
dat2<-dat2[order(dat2$IID,-abs(dat2$gdm)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,mom,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","gdm")]
write.table(dat2,"./maternal_gdm.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_mom","stillbirth_mom")]
dat2<-rename.vars(dat2,from=c("sentrixid_mom","stillbirth_mom"),to=c("IID","stillbirth"))
dat2<-dat2[order(dat2$IID,-abs(dat2$stillbirth)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,mom,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","stillbirth")]
write.table(dat2,"./maternal_stillbirth.txt",sep="\t",row.names=FALSE, quote=FALSE)


load("N:/data/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")
dat<-subset2(dat,"!is.na(dat$cad_grs_dad)")

dat2<-dat[,c("sentrixid_dad","miscarriage_mbrn_dad")]
dat2<-rename.vars(dat2,from=c("sentrixid_dad","miscarriage_mbrn_dad"),to=c("IID","miscarriage_mbrn"))
dat2<-dat2[order(dat2$IID,-abs(dat2$miscarriage_mbrn)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,dad,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","miscarriage_mbrn")]
write.table(dat2,"./paternal_miscarriage_mbrn.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_dad","eclampsia_dad")]
dat2<-rename.vars(dat2,from=c("sentrixid_dad","eclampsia_dad"),to=c("IID","eclampsia"))
dat2<-dat2[order(dat2$IID,-abs(dat2$eclampsia)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,dad,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","eclampsia")]
write.table(dat2,"./paternal_eclampsia.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_dad","hta_preg_dad")]
dat2<-rename.vars(dat2,from=c("sentrixid_dad","hta_preg_dad"),to=c("IID","htapreg"))
dat2<-dat2[order(dat2$IID,-abs(dat2$htapreg)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,dad,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","htapreg")]
write.table(dat2,"./paternal_htapreg.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_dad","preterm_dad")]
dat2<-rename.vars(dat2,from=c("sentrixid_dad","preterm_dad"),to=c("IID","preterm"))
dat2<-dat2[order(dat2$IID,-abs(dat2$preterm)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,dad,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","preterm")]
write.table(dat2,"./paternal_preterm.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_dad","preterm_sp_dad")]
dat2<-rename.vars(dat2,from=c("sentrixid_dad","preterm_sp_dad"),to=c("IID","preterm_sp"))
dat2<-dat2[order(dat2$IID,-abs(dat2$preterm_sp)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,dad,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","preterm_sp")]
write.table(dat2,"./paternal_preterm_sp.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_dad","sga_dad")]
dat2<-rename.vars(dat2,from=c("sentrixid_dad","sga_dad"),to=c("IID","sga"))
dat2<-dat2[order(dat2$IID,-abs(dat2$sga)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,dad,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","sga")]
write.table(dat2,"./paternal_sga.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_dad","lga_dad")]
dat2<-rename.vars(dat2,from=c("sentrixid_dad","lga_dad"),to=c("IID","lga"))
dat2<-dat2[order(dat2$IID,-abs(dat2$lga)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,dad,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","lga")]
write.table(dat2,"./paternal_lga.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_dad","gdm_dad")]
dat2<-rename.vars(dat2,from=c("sentrixid_dad","gdm_dad"),to=c("IID","gdm"))
dat2<-dat2[order(dat2$IID,-abs(dat2$gdm)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,dad,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","gdm")]
write.table(dat2,"./paternal_gdm.txt",sep="\t",row.names=FALSE, quote=FALSE)

dat2<-dat[,c("sentrixid_dad","stillbirth_dad")]
dat2<-rename.vars(dat2,from=c("sentrixid_dad","stillbirth_dad"),to=c("IID","stillbirth"))
dat2<-dat2[order(dat2$IID,-abs(dat2$stillbirth)),]
dat2<-dat2[!duplicated(dat2$IID),]
dat2<-merge2(dat2,dad,by.id=c("IID"),all.x=TRUE,sort=FALSE)
dat2<-dat2[,c("FID","IID","stillbirth")]
write.table(dat2,"./paternal_stillbirth.txt",sep="\t",row.names=FALSE, quote=FALSE)


#####################
### GWAS IN LINUX ###
#####################

##########################
### CLEAN GWAS RESULTS ###
##########################

setwd("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/preg_outc/")
z<-qnorm(1-0.05/2)
freq<-fread("./results/maternal.afreq",header=TRUE,sep="\t",sep2="\t")
freq<-freq[,c("ID","ALT_FREQS")]
names(freq)<-c("SNP","eaf")

dat<-fread("./results/maternal_htapreg.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Hypertensive disorder of pregnancy, mothers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/htapreg_mom.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/maternal_eclampsia.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Preeclampsia and eclampsia, mothers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/eclampsia_mom.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/maternal_miscarriage_mbrn.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Miscarriage, mothers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/miscarriage_mbrn_mom.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/maternal_stillbirth.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Stillbirth, mothers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/stillbirth_mom.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/maternal_preterm.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Preterm birth, mothers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/preterm_mom.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/maternal_preterm_sp.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Spontaneous preterm birth, mothers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/preterm_sp_mom.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/maternal_sga.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Small for gestational age, mothers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/sga_mom.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/maternal_lga.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Large for gestational age, mothers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/lga_mom.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)


dat<-fread("./results/paternal_htapreg.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Hypertensive disorder of pregnancy, fathers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/htapreg_dad.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/paternal_eclampsia.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Preeclampsia and eclampsia, fathers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/eclampsia_dad.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/paternal_miscarriage_mbrn.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Miscarriage, fathers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/miscarriage_mbrn_dad.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/paternal_stillbirth.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Stillbirth, fathers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/stillbirth_dad.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/paternal_preterm.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Preterm birth, fathers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/preterm_dad.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/paternal_preterm_sp.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Spontaneous preterm birth, fathers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/preterm_sp_dad.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/paternal_sga.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Small for gestational age, fathers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/sga_dad.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/paternal_lga.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Large for gestational age, fathers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/lga_dad.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

dat<-fread("./results/paternal_gdm.txt",header=FALSE,sep="\t",sep2="\t")
names(dat)<-c("chrom","pos","SNP","other_allele","effect_allele","effect_allele2","test","samplesize","or","se","zstat","pval")
dat<-merge(dat,freq,by="SNP",all.x=TRUE)
dat$beta<-log(dat$or)
dat$Units<-c("Units")
dat$Phenotype<-c("Gestational diabetes, fathers")
dat<-dat[,c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dat,"./definitive/gdm_dad.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)

moba_all<-fread("N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_translate.txt",header=TRUE,sep="\t",sep2="\t")
hunt_all<-fread("N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/hunt_translate.txt",header=TRUE,sep="\t",sep2="\t")
alspac_all<-fread("N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/alspac_translate.txt",header=TRUE,sep="\t",sep2="\t")

dat<-fread("./results/maternal.gdm.fixedP",header=TRUE,sep="\t",sep2="\t")
dat<-dat[,c("MarkerName","Chromosome","Position","EA","NEA","EAF","beta_0","se_0","P.value_association","Nsample")]
dat<-dat[!is.na(dat$beta_0) & !is.na(dat$se_0) & !is.na(dat$P.value_association) & !is.na(dat$EAF),]
dat<-dat[dat$EAF>0.01,]
dat<-dat[dat$EAF<0.99,]
dat<-rename.vars(dat, from=c("MarkerName"), to=c("MARKERNAME"))
dat$MARKERNAME<-paste(paste(dat$Chromosome,dat$Position,sep=":"),
                      paste(dat$NEA,dat$EA,sep="/"),sep="_")
dat<-merge(dat,moba_all,by="MARKERNAME",all.x=TRUE)
dat$MARKERNAME<-with(dat,ifelse(!is.na(RSID),RSID,MARKERNAME))
dat$RSID<-NULL
moba_all<-NULL
dat<-merge(dat,hunt_all,by="MARKERNAME",all.x=TRUE)
dat$MARKERNAME<-with(dat,ifelse(!is.na(RSID),RSID,MARKERNAME))
dat$RSID<-NULL
hunt_all<-NULL
dat<-merge(dat,alspac_all,by="MARKERNAME",all.x=TRUE)
dat$MARKERNAME<-with(dat,ifelse(!is.na(RSID),RSID,MARKERNAME))
dat$RSID<-NULL
alspac_all<-NULL
dat<-dat[startsWith(dat$MARKERNAME,"rs"),]

names(dat)<-c("SNP","chrom","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize")
dat$Units<-c("Units")
dat$Phenotype<-c("Gestational diabetes, mothers")
fwrite(dat,"./definitive/gdm_mom.txt",append=FALSE,sep="\t",sep2=c("\t","|","\t"),row.names=FALSE, col.names=TRUE)


##########################################
### TWO-SAMPLE MENDELIAN RANDOMIZATION ###
##########################################

dir.create("N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr")

setwd("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/preg_outc/definitive")

# EXPOSURE : CAD #

cad<-read.csv2("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/cad/gwas_cad_varderharst2018.csv",
              header=TRUE,sep=";",dec=".")
cad$Phenotype<-c("Coronary heart disease")
cad$units<-c("units")
cad$id<-c("Van Der Harst P, 2018")
cad$samplesize<-547261
cad<-rename.vars(cad,
                from=c("chr","rsid"),
                to=c("chrom","SNP"))
cad<-cad[,c("Phenotype","SNP","beta","se","eaf","effect_allele","other_allele","pval",
          "units","samplesize","id")]
write.table(cad,"./cad.txt",
            sep="\t",row.names=FALSE, quote=FALSE)

cad_harm<-read_exposure_data(filename = "./cad.txt",
                            clump = FALSE,
                            sep="\t",
                            phenotype_col = "Phenotype",
                            snp_col = "SNP",
                            beta_col = "beta",
                            se_col = "se",
                            eaf_col = "eaf",
                            effect_allele_col = "effect_allele",
                            other_allele_col = "other_allele",
                            pval_col = "pval",
                            units_col = "units",
                            samplesize_col = "samplesize",
                            id_col = "id",
                            min_pval = 1e-200,
                            log_pval = FALSE)
save(cad_harm,file="./cad_harm.RData")


# OUTCOME : Miscarriage - Mothers #

miscarriage_mbrn_mom_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./miscarriage_mbrn_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, miscarriage_mbrn_mom_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_miscarriage_mbrn_mom_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_miscarriage_mbrn_mom_withpalindromic.RData")

dat<-harmonise_data(cad_harm, miscarriage_mbrn_mom_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_miscarriage_mbrn_mom_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_miscarriage_mbrn_mom_nopalindromic.RData")


# OUTCOME : stillbirth - Mothers #

stillbirth_mom_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./stillbirth_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, stillbirth_mom_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_stillbirth_mom_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_stillbirth_mom_withpalindromic.RData")

dat<-harmonise_data(cad_harm, stillbirth_mom_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_stillbirth_mom_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_stillbirth_mom_nopalindromic.RData")


# OUTCOME : HTA disorders in pregnancy - Mothers #

htapreg_mom_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./htapreg_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, htapreg_mom_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_htapreg_mom_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_htapreg_mom_withpalindromic.RData")

dat<-harmonise_data(cad_harm, htapreg_mom_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_htapreg_mom_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_htapreg_mom_nopalindromic.RData")


# OUTCOME : Preeclampsia and eclampsia - Mothers #

eclampsia_mom_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./eclampsia_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, eclampsia_mom_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_eclampsia_mom_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_eclampsia_mom_withpalindromic.RData")

dat<-harmonise_data(cad_harm, eclampsia_mom_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_eclampsia_mom_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_eclampsia_mom_nopalindromic.RData")


# OUTCOME : Preterm - Mothers #

preterm_mom_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./preterm_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, preterm_mom_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_preterm_mom_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_preterm_mom_withpalindromic.RData")

dat<-harmonise_data(cad_harm, preterm_mom_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_preterm_mom_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_preterm_mom_nopalindromic.RData")


# OUTCOME : Spontaneous preterm - Mothers #

preterm_sp_mom_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./preterm_sp_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, preterm_sp_mom_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_preterm_sp_mom_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_preterm_sp_mom_withpalindromic.RData")

dat<-harmonise_data(cad_harm, preterm_sp_mom_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_preterm_sp_mom_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_preterm_sp_mom_nopalindromic.RData")


# OUTCOME : SGA - Mothers #

sga_mom_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./sga_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, sga_mom_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_sga_mom_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_sga_mom_withpalindromic.RData")

dat<-harmonise_data(cad_harm, sga_mom_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_sga_mom_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_sga_mom_nopalindromic.RData")


# OUTCOME : LGA - Mothers #

lga_mom_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./lga_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, lga_mom_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_lga_mom_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_lga_mom_withpalindromic.RData")

dat<-harmonise_data(cad_harm, lga_mom_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_lga_mom_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_lga_mom_nopalindromic.RData")


# OUTCOME : Gestational diabetes - Mothers #

gdm_mom_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./gdm_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, gdm_mom_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_gdm_mom_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_gdm_mom_withpalindromic.RData")

dat<-harmonise_data(cad_harm, gdm_mom_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_gdm_mom_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_gdm_mom_nopalindromic.RData")


# OUTCOME : Miscarriage - Fathers #

miscarriage_mbrn_dad_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./miscarriage_mbrn_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, miscarriage_mbrn_dad_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_miscarriage_mbrn_dad_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_miscarriage_mbrn_dad_withpalindromic.RData")

dat<-harmonise_data(cad_harm, miscarriage_mbrn_dad_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_miscarriage_mbrn_dad_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_miscarriage_mbrn_dad_nopalindromic.RData")


# OUTCOME : stillbirth - Fathers #

stillbirth_dad_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./stillbirth_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, stillbirth_dad_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_stillbirth_dad_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_stillbirth_dad_withpalindromic.RData")

dat<-harmonise_data(cad_harm, stillbirth_dad_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_stillbirth_dad_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_stillbirth_dad_nopalindromic.RData")


# OUTCOME : HTA disorders in pregnancy - Fathers #

htapreg_dad_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./htapreg_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, htapreg_dad_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_htapreg_dad_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_htapreg_dad_withpalindromic.RData")

dat<-harmonise_data(cad_harm, htapreg_dad_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_htapreg_dad_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_htapreg_dad_nopalindromic.RData")


# OUTCOME : Preeclampsia and eclampsia - Fathers #

eclampsia_dad_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./eclampsia_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, eclampsia_dad_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_eclampsia_dad_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_eclampsia_dad_withpalindromic.RData")

dat<-harmonise_data(cad_harm, eclampsia_dad_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_eclampsia_dad_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_eclampsia_dad_nopalindromic.RData")


# OUTCOME : Preterm - Fathers #

preterm_dad_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./preterm_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, preterm_dad_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_preterm_dad_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_preterm_dad_withpalindromic.RData")

dat<-harmonise_data(cad_harm, preterm_dad_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_preterm_dad_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_preterm_dad_nopalindromic.RData")


# OUTCOME : Preterm - Fathers #

preterm_sp_dad_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./preterm_sp_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, preterm_sp_dad_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_preterm_sp_dad_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_preterm_sp_dad_withpalindromic.RData")

dat<-harmonise_data(cad_harm, preterm_sp_dad_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_preterm_sp_dad_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_preterm_sp_dad_nopalindromic.RData")


# OUTCOME : SGA - Fathers #

sga_dad_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./sga_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, sga_dad_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_sga_dad_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_sga_dad_withpalindromic.RData")

dat<-harmonise_data(cad_harm, sga_dad_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_sga_dad_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_sga_dad_nopalindromic.RData")


# OUTCOME : LGA - Fathers #

lga_dad_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./lga_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, lga_dad_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_lga_dad_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_lga_dad_withpalindromic.RData")

dat<-harmonise_data(cad_harm, lga_dad_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_lga_dad_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_lga_dad_nopalindromic.RData")


# OUTCOME : Gestational diabetes - Fathers #

gdm_dad_harm <- read_outcome_data(
  snps = cad$SNP,
  filename = "./gdm_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(cad_harm, gdm_dad_harm, action = 1)
write.table(dat$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snps_cad_gdm_dad_withpalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_gdm_dad_withpalindromic.RData")

dat<-harmonise_data(cad_harm, gdm_dad_harm, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr/snp_cad_gdm_dad_nopalindromic.csv",
            sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./cad_gdm_dad_nopalindromic.RData")


##############################
### TWO-SAMPLE MR ANALYSES ###
##############################

setwd("N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/2smr")

vars01<-c("miscarriage_mbrn","miscarriage_mbrn","stillbirth","stillbirth","eclampsia","eclampsia","htapreg","htapreg",
          "preterm","preterm","preterm_sp","preterm_sp","sga","sga","lga","lga","gdm","gdm",
          "miscarriage_mbrn","miscarriage_mbrn","stillbirth","stillbirth","eclampsia","eclampsia","htapreg","htapreg",
          "preterm","preterm","preterm_sp","preterm_sp","sga","sga","lga","lga","gdm","gdm")
vars02<-c("mom","mom","mom","mom","mom","mom","mom","mom","mom","mom","mom","mom","mom","mom","mom","mom","mom","mom",
          "dad","dad","dad","dad","dad","dad","dad","dad","dad","dad","dad","dad","dad","dad","dad","dad","dad","dad")
vars03<-c("_withpalindromic","_nopalindromic","_withpalindromic","_nopalindromic",
          "_withpalindromic","_nopalindromic","_withpalindromic","_nopalindromic",
          "_withpalindromic","_nopalindromic","_withpalindromic","_nopalindromic",
          "_withpalindromic","_nopalindromic","_withpalindromic","_nopalindromic",
          "_withpalindromic","_nopalindromic","_withpalindromic","_nopalindromic",
          "_withpalindromic","_nopalindromic","_withpalindromic","_nopalindromic",
          "_withpalindromic","_nopalindromic","_withpalindromic","_nopalindromic",
          "_withpalindromic","_nopalindromic","_withpalindromic","_nopalindromic",
          "_withpalindromic","_nopalindromic","_withpalindromic","_nopalindromic")

tab<-NULL
for(i in 1:length(vars01))
  
{
  namedat<-paste("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/preg_outc/definitive/cad_",vars01[i],
                 "_",vars02[i],vars03[i],".RData",sep="")
  load(namedat)
  
  mr_results<-mr(dat,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
  mr_results$pval<-pval_guapa(mr_results$pval)
  mr_results$beta<-paste(guapa(mr_results$b)," (",guapa(mr_results$se),")",sep="")
  mr_results$or<-risk_se_ic_guapa(mr_results$b,mr_results$se)
  forest<-risk_se_ic_guapa2(mr_results$b[1],mr_results$se[1])
  mr_raps<-mr.raps(b_exp=dat$beta.exposure,b_out=dat$beta.outcome,se_exp=dat$se.exposure,se_out=dat$se.outcome,diagnosis=FALSE)
  mr_raps_beta<-paste(guapa(mr_raps$beta.hat)," (",guapa(mr_raps$beta.se),")",sep="")
  mr_raps_or<-risk_se_ic_guapa(mr_raps$beta.hat,mr_raps$beta.se)
  mr_raps_pval<-pval_guapa(mr_raps$beta.p.value)
  q1<-paste(round(mr_heterogeneity(dat)$Q[2],2)," (P=",pval_guapa(mr_heterogeneity(dat)$Q_pval[2]),")",sep="")
  q2<-paste(round(mr_heterogeneity(dat)$Q[1],2)," (P=",pval_guapa(mr_heterogeneity(dat)$Q_pval[1]),")",sep="")
  plei<-pval_guapa(mr_pleiotropy_test(dat)$pval)
  
  tab<-rbind(tab,cbind(mr_results$nsnp[1],
                       mr_results$beta[1],mr_results$or[1],mr_results$pval[1],
                       mr_results$beta[2],mr_results$or[2],mr_results$pval[2],
                       mr_results$beta[3],mr_results$or[3],mr_results$pval[3],
                       mr_results$beta[4],mr_results$or[4],mr_results$pval[4],
                       mr_raps_beta,mr_raps_or,mr_raps_pval,plei,q1,q2,forest))
}

colnames(tab)<-c("SNPs","IVW_b","IVW_or","IVW_p","Egger_b","Egger_or","Egger_p",
                 "WMe_b","WMe_or","Wme_p","WMo_b","WMo_or","WMo_p",
                 "RAPS_b","RAPS_or","RAPS_p","Egger_pleio","Cochran_Q","Rucker_Q","forestplot")
rownames(tab)<-paste("cad_",vars01,"_",vars02,vars03,sep="")
write.table(tab,file="./2smr.csv",sep=";",col.names=NA)



### APPROXIMATION TO F-STATISTIC AND UNWEIGHTED I-SQUARED ###
#############################################################

# F>10: good IVW performance
# Isq close to 1: good MR-Egger performance

Isq = function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

load("./cad_harm.RData")
fstat<-mean((abs(cad_harm$beta.exposure))^2/cad_harm$se.exposure^2,na.rm=TRUE)
unIsq<-Isq(abs(cad_harm$beta.exposure),cad_harm$se.exposure)





####################################################################################################

### LOGISTIC REGRESSION: NON-LINEAR ASSOCIATIONS (SMOOTHED SPLINES) ###
#######################################################################

# EXECUTE MANUALLY: i=1, then the loop #

load("N:/data/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")
dat$adj<-0
dat$preg_plan2<-with(dat,ifelse(preg_plan==1,0,1))

vars00<-c("cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad",
          "cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad",
          "cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad",
          "stroke_grs_mom","stroke_grs_dad","stroke_grs_mom","stroke_grs_dad",
          "stroke_grs_mom","stroke_grs_dad","stroke_grs_mom","stroke_grs_dad",
          "stroke_grs_mom","stroke_grs_dad","stroke_grs_mom","stroke_grs_dad","stroke_grs_mom","stroke_grs_dad")
vars01<-c("cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z",
          "cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z",
          "cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z",
          "stroke_grs_mom_z","stroke_grs_dad_z","stroke_grs_mom_z","stroke_grs_dad_z",
          "stroke_grs_mom_z","stroke_grs_dad_z","stroke_grs_mom_z","stroke_grs_dad_z",
          "stroke_grs_mom_z","stroke_grs_dad_z","stroke_grs_mom_z","stroke_grs_dad_z","stroke_grs_mom_z","stroke_grs_dad_z")
vars02<-c("subf_mom","subf_dad","miscarriage_mom","miscarriage_dad",
          "eclampsia_mom","eclampsia_dad","hta_preg_mom","hta_preg_dad",
          "preterm_mom","preterm_dad","sga_mom","sga_dad","gdm_mom","gdm_dad",
          "subf_mom","subf_dad","miscarriage_mom","miscarriage_dad",
          "eclampsia_mom","eclampsia_dad","hta_preg_mom","hta_preg_dad",
          "preterm_mom","preterm_dad","sga_mom","sga_dad","gdm_mom","gdm_dad")
vars03<-c("preg_plan2","preg_plan2","adj","adj",
          "adj","adj","hta_chronic_mom","hta_chronic_mom",
          "adj","adj","adj","adj","adj","adj",
          "preg_plan2","preg_plan2","adj","adj",
          "adj","adj","hta_chronic_mom","hta_chronic_mom",
          "adj","adj","adj","adj","adj","adj")
vars04<-c("exclude_genotype_mom","exclude_genotype_dad","exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad","exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad","exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad","exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad","exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad","exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad")
vars08<-c("pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad",
          "pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad",
          "pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad")
vars09<-c("pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad",
          "pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad",
          "pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad")
vars10<-c("pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad",
          "pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad",
          "pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad")
vars11<-c("pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad",
          "pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad",
          "pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad")
vars12<-c("pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad",
          "pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad",
          "pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad")
vars13<-c("pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad",
          "pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad",
          "pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad")
vars14<-c("pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad",
          "pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad",
          "pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad")
vars15<-c("pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad",
          "pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad",
          "pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad")
vars16<-c("pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad",
          "pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad",
          "pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad")
vars17<-c("pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad",
          "pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad",
          "pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad")
vars18<-c("m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374",
          "m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374",
          "m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374")
vars20<-c("CAD-GRS, mothers","CAD-GRS, fathers",
          "CAD-GRS, mothers","CAD-GRS, fathers",
          "CAD-GRS, mothers","CAD-GRS, fathers",
          "CAD-GRS, mothers","CAD-GRS, fathers",
          "CAD-GRS, mothers","CAD-GRS, fathers",
          "CAD-GRS, mothers","CAD-GRS, fathers",
          "CAD-GRS, mothers","CAD-GRS, fathers",
          "Stroke-GRS, mothers","Stroke-GRS, fathers",
          "Stroke-GRS, mothers","Stroke-GRS, fathers",
          "Stroke-GRS, mothers","Stroke-GRS, fathers",
          "Stroke-GRS, mothers","Stroke-GRS, fathers",
          "Stroke-GRS, mothers","Stroke-GRS, fathers",
          "Stroke-GRS, mothers","Stroke-GRS, fathers",
          "Stroke-GRS, mothers","Stroke-GRS, fathers")
vars21<-c("Subfertility (odds ratio, adjusted)","Subfertility (odds ratio, adjusted)",
          "Miscarriage (odds ratio, adjusted)","Miscarriage (odds ratio, adjusted)",
          "Pre-/Eclampsia (odds ratio, adjusted)","Pre-/Eclampsia (odds ratio, adjusted)",
          "HT in pregnancy (odds ratio, adjusted)","HT in pregnancy (odds ratio, adjusted)",
          "Preterm birth (odds ratio, adjusted)","Preterm birth (odds ratio, adjusted)",
          "SGA (odds ratio, adjusted)","SGA (odds ratio, adjusted)",
          "Gestat. diabetes (odds ratio, adjusted)","Gestat. diabetes (odds ratio, adjusted)",
          "Subfertility (odds ratio, adjusted)","Subfertility (odds ratio, adjusted)",
          "Miscarriage (odds ratio, adjusted)","Miscarriage (odds ratio, adjusted)",
          "Pre-/Eclampsia (odds ratio, adjusted)","Pre-/Eclampsia (odds ratio, adjusted)",
          "HT in pregnancy (odds ratio, adjusted)","HT in pregnancy (odds ratio, adjusted)",
          "Preterm birth (odds ratio, adjusted)","Preterm birth (odds ratio, adjusted)",
          "SGA (odds ratio, adjusted)","SGA (odds ratio, adjusted)",
          "Gestat. diabetes (odds ratio, adjusted)","Gestat. diabetes (odds ratio, adjusted)")

tab<-NULL
for(i in 1:length(vars01))
  
{
  varstot<-c(vars00[i],vars02[i],vars03[i],vars04[i],vars08[i],vars09[i],vars10[i],vars11[i],vars12[i],vars13[i],
             vars14[i],vars15[i],vars16[i],vars17[i],vars18[i],"preg_plan")
  dat2<-na.omit(dat[,varstot])
  dat2<-subset2(dat2,"dat2[,vars03[i]]==0 & dat2[,vars04[i]]==0")
  dat2<-dat2[!duplicated(dat2[,vars18[i]]),]
  dat2[,vars01[i]]<-as.numeric(with(dat2,scale(dat2[,vars00[i]])))
  aaa<-dat2[,vars01[i]]
  minim<-mean(aaa,na.rm=TRUE)
  
  mod01<-gam(formula=as.factor(dat2[,vars02[i]])~bs(aaa,df=4)
             +dat2[,vars08[i]]+dat2[,vars09[i]]+dat2[,vars10[i]]+dat2[,vars11[i]]
             +dat2[,vars12[i]]+dat2[,vars13[i]]+dat2[,vars14[i]]
             +dat2[,vars15[i]]+dat2[,vars16[i]]+dat2[,vars17[i]],
             data=dat2, family="binomial")
  ptemp<-termplot(mod01,term=1,se=TRUE,plot=FALSE)
  temp<-ptemp$aaa
  value<-closest(temp$x,minim)
  #min_val<-temp$x[which(temp$y==min(temp$y,na.rm=TRUE))]
  center<-with(temp, y[x==value])
  z<-qnorm(1-0.05/2)
  ytemp<-temp$y+outer(temp$se,c(0,-z,z),'*')
  min_val<-guapa(temp$x[which(temp$y==min(temp$y,na.rm=TRUE))])
  ci<-exp(ytemp-center)
  name<-paste("./results/",vars01[i],"_",vars02[i],".jpg",sep="")
  labely<-vars21[i]
  
  plot.data<-as.data.frame(cbind(temp$x,ci))
  colnames(plot.data)<-c("xval","yest","lci","uci")
  ft<-as.data.frame(table(dat2[,vars01[i]]))
  names(ft)<-c("xval","freq")
  plot.data<-merge2(plot.data,ft,by.id=c("xval"),all.x=TRUE,sort=FALSE)
  infl01<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[1],]$xval
  infl02<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[2],]$xval
  infl03<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[3],]$xval
  infl04<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[4],]$xval
  infl05<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[5],]$xval
  infl06<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[6],]$xval
  infl07<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[7],]$xval
  infl08<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[8],]$xval
  infl09<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[9],]$xval
  infl10<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[10],]$xval
  
  figure<-ggplot() +
    geom_histogram(data=dat2, aes(x=dat2[,vars01[i]], y=(..density../max(..density..,na.rm=TRUE))*3.75), 
                   bins=20, color="grey80", fill="grey90") +
    #scale_x_continuous(limits = c(min(dat2[,vars01[i]],na.rm=TRUE),max(dat2[,vars01[i]],na.rm=TRUE))) +
    scale_y_continuous(limits = c(0,4), sec.axis = sec_axis(~., name = "Relative density")) +
    geom_hline(aes(yintercept=1), data=plot.data, colour="black", linetype=2) + 
    geom_line(aes(x=xval, y=yest), data=plot.data, color="black") + 
    geom_line(aes(x=xval, y=lci), data=plot.data, color="grey35") + 
    geom_line(aes(x=xval ,y=uci), data=plot.data, color="grey35") + 
    theme_bw() +
    labs(x=vars20[i],y=labely) +
    theme(axis.title.x = element_text(vjust=0.5, size=20), 
          axis.title.y = element_text(vjust=0.5, size=20),
          axis.text.x = element_text(size=16, colour = 'black'),
          axis.text.y = element_text(size=16, colour = 'black'),
          axis.ticks.x = element_line(colour = 'black'),
          axis.ticks.y = element_line(colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank())  
  #    geom_point(aes(x=infl01, y=with(plot.data, yest[xval==infl01])), data=plot.data, colour="black", size=2.5) +
  #    geom_point(aes(x=infl02, y=with(plot.data, yest[xval==infl02])), data=plot.data, colour="black", size=2.5)
  
  jpeg(filename=name,width = 8000, height = 8000, res=1200)
  par(las=1,cex=1.2,mar=c(6,6,2,0),bty="n",lheight=0.9)
  figure
  dev.off()
  
  tab<-rbind(tab,cbind(infl01,infl02,infl03,infl04,infl05,infl06,infl07,infl08,infl09,infl10))
}

rownames(tab)<-vars01
write.table(tab,file="./results/inflection_points.csv",sep=";",col.names=NA)


### MENDELIAN RANDOMIZATION IN STRATA ###
#########################################

load("N:/data/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")
dat$adj<-0
dat$preg_plan2<-with(dat,ifelse(preg_plan==1,0,1))

vars00<-c("cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad",
          "cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad",
          "cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad",
          "stroke_grs_mom","stroke_grs_dad","stroke_grs_mom","stroke_grs_dad",
          "stroke_grs_mom","stroke_grs_dad","stroke_grs_mom","stroke_grs_dad",
          "stroke_grs_mom","stroke_grs_dad","stroke_grs_mom","stroke_grs_dad","stroke_grs_mom","stroke_grs_dad")
vars01<-c("cad_grs_mom_q","cad_grs_dad_q","cad_grs_mom_q","cad_grs_dad_q",
          "cad_grs_mom_q","cad_grs_dad_q","cad_grs_mom_q","cad_grs_dad_q",
          "cad_grs_mom_q","cad_grs_dad_q","cad_grs_mom_q","cad_grs_dad_q","cad_grs_mom_q","cad_grs_dad_q",
          "stroke_grs_mom_q","stroke_grs_dad_q","stroke_grs_mom_q","stroke_grs_dad_q",
          "stroke_grs_mom_q","stroke_grs_dad_q","stroke_grs_mom_q","stroke_grs_dad_q",
          "stroke_grs_mom_q","stroke_grs_dad_q","stroke_grs_mom_q","stroke_grs_dad_q","stroke_grs_mom_q","stroke_grs_dad_q")
vars02<-c("cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z",
          "cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z",
          "cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z",
          "stroke_grs_mom_z","stroke_grs_dad_z","stroke_grs_mom_z","stroke_grs_dad_z",
          "stroke_grs_mom_z","stroke_grs_dad_z","stroke_grs_mom_z","stroke_grs_dad_z",
          "stroke_grs_mom_z","stroke_grs_dad_z","stroke_grs_mom_z","stroke_grs_dad_z","stroke_grs_mom_z","stroke_grs_dad_z")
vars03<-c("subf_mom","subf_dad","miscarriage_mom","miscarriage_dad",
          "eclampsia_mom","eclampsia_dad","hta_preg_mom","hta_preg_dad",
          "preterm_mom","preterm_dad","sga_mom","sga_dad","gdm_mom","gdm_dad",
          "subf_mom","subf_dad","miscarriage_mom","miscarriage_dad",
          "eclampsia_mom","eclampsia_dad","hta_preg_mom","hta_preg_dad",
          "preterm_mom","preterm_dad","sga_mom","sga_dad","gdm_mom","gdm_dad")
vars04<-c("preg_plan2","preg_plan2","adj","adj",
          "adj","adj","hta_chronic_mom","hta_chronic_mom",
          "adj","adj","adj","adj","adj","adj",
          "preg_plan2","preg_plan2","adj","adj",
          "adj","adj","hta_chronic_mom","hta_chronic_mom",
          "adj","adj","adj","adj","adj","adj")
vars05<-c("exclude_genotype_mom","exclude_genotype_dad","exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad","exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad","exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad","exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad","exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad","exclude_genotype_mom","exclude_genotype_dad",
          "exclude_genotype_mom","exclude_genotype_dad")
vars08<-c("pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad",
          "pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad",
          "pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad")
vars09<-c("pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad",
          "pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad",
          "pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad")
vars10<-c("pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad",
          "pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad",
          "pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad")
vars11<-c("pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad",
          "pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad",
          "pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad")
vars12<-c("pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad",
          "pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad",
          "pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad")
vars13<-c("pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad",
          "pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad",
          "pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad")
vars14<-c("pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad",
          "pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad",
          "pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad")
vars15<-c("pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad",
          "pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad",
          "pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad")
vars16<-c("pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad",
          "pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad",
          "pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad")
vars17<-c("pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad",
          "pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad",
          "pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad")
vars18<-c("m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374",
          "m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374",
          "m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374")
vars20<-c("CAD-GRS - Subfertility, mothers","CAD-GRS - Subfertility, fathers",
          "CAD-GRS - Miscarriage, mothers","CAD-GRS - Miscarriage, fathers",
          "CAD-GRS - Pre-/Eclampsia, mothers","CAD-GRS - Pre-/Eclampsia, fathers",
          "CAD-GRS - HT in pregnancy, mothers","CAD-GRS - HT in pregnancy, fathers",
          "CAD-GRS - Preterm, mothers","CAD-GRS - Preterm, fathers",
          "CAD-GRS - SGA, mothers","CAD-GRS - SGA, fathers",
          "CAD-GRS - Gestat. diabetes, mothers","CAD-GRS - Gestat. diabetes, fathers",
          "Stroke-GRS - Subfertility, mothers","Stroke-GRS - Subfertility, fathers",
          "Stroke-GRS - Miscarriage, mothers","Stroke-GRS - Miscarriage, fathers",
          "Stroke-GRS - Pre-/Eclampsia, mothers","Stroke-GRS - Pre-/Eclampsia, fathers",
          "Stroke-GRS - HT in pregnancy, mothers","Stroke-GRS - HT in pregnancy, fathers",
          "Stroke-GRS - Preterm, mothers","Stroke-GRS - Preterm, fathers",
          "Stroke-GRS - SGA, mothers","Stroke-GRS - SGA, fathers",
          "Stroke-GRS - Gestat. diabetes, mothers","Stroke-GRS - Gestat. diabetes, fathers")

tab<-NULL
for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"!is.na(dat[,vars00[i]]) & dat[,vars04[i]]==0 & dat[,vars05[i]]==0")
  datx<-datx[!duplicated(datx[,vars18[i]]),]
  datx[,vars01[i]]<-as.numeric(ntile(datx[,vars00[i]], 5))
  datx[,vars02[i]]<-as.numeric(with(datx,scale(datx[,vars00[i]])))
  
  dat01<-subset2(datx,"datx[,vars01[i]]==1")
  dat02<-subset2(datx,"datx[,vars01[i]]==2")
  dat03<-subset2(datx,"datx[,vars01[i]]==3")
  dat04<-subset2(datx,"datx[,vars01[i]]==4")
  dat05<-subset2(datx,"datx[,vars01[i]]==5")
  
  mod01<-miceadds::glm.cluster(formula=as.factor(dat01[,vars03[i]])~dat01[,vars02[i]]
                               +dat01[,vars08[i]]+dat01[,vars09[i]]+dat01[,vars10[i]]+dat01[,vars11[i]]
                               +dat01[,vars12[i]]+dat01[,vars13[i]]+dat01[,vars14[i]]
                               +dat01[,vars15[i]]+dat01[,vars16[i]]+dat01[,vars17[i]],
                               data=dat01, cluster=dat01[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod01)[2,1])
  se<-as.numeric(summary(mod01)[2,2])
  coef01<-risk_se_ic_guapa(estimate,se)
  pval01<-pval_guapa(as.numeric(summary(mod01)[2,4]))
  
  mod02<-miceadds::glm.cluster(formula=as.factor(dat02[,vars03[i]])~dat02[,vars02[i]]
                               +dat02[,vars08[i]]+dat02[,vars09[i]]+dat02[,vars10[i]]+dat02[,vars11[i]]
                               +dat02[,vars12[i]]+dat02[,vars13[i]]+dat02[,vars14[i]]
                               +dat02[,vars15[i]]+dat02[,vars16[i]]+dat02[,vars17[i]],
                               data=dat02, cluster=dat02[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod02)[2,1])
  se<-as.numeric(summary(mod02)[2,2])
  coef02<-risk_se_ic_guapa(estimate,se)
  pval02<-pval_guapa(as.numeric(summary(mod02)[2,4]))
  
  mod03<-miceadds::glm.cluster(formula=as.factor(dat03[,vars03[i]])~dat03[,vars02[i]]
                               +dat03[,vars08[i]]+dat03[,vars09[i]]+dat03[,vars10[i]]+dat03[,vars11[i]]
                               +dat03[,vars12[i]]+dat03[,vars13[i]]+dat03[,vars14[i]]
                               +dat03[,vars15[i]]+dat03[,vars16[i]]+dat03[,vars17[i]],
                               data=dat03, cluster=dat03[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod03)[2,1])
  se<-as.numeric(summary(mod03)[2,2])
  coef03<-risk_se_ic_guapa(estimate,se)
  pval03<-pval_guapa(as.numeric(summary(mod03)[2,4]))
  
  mod04<-miceadds::glm.cluster(formula=as.factor(dat04[,vars03[i]])~dat04[,vars02[i]]
                               +dat04[,vars08[i]]+dat04[,vars09[i]]+dat04[,vars10[i]]+dat04[,vars11[i]]
                               +dat04[,vars12[i]]+dat04[,vars13[i]]+dat04[,vars14[i]]
                               +dat04[,vars15[i]]+dat04[,vars16[i]]+dat04[,vars17[i]],
                               data=dat04, cluster=dat04[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod04)[2,1])
  se<-as.numeric(summary(mod04)[2,2])
  coef04<-risk_se_ic_guapa(estimate,se)
  pval04<-pval_guapa(as.numeric(summary(mod04)[2,4]))
  
  mod05<-miceadds::glm.cluster(formula=as.factor(dat05[,vars03[i]])~dat05[,vars02[i]]
                               +dat05[,vars08[i]]+dat05[,vars09[i]]+dat05[,vars10[i]]+dat05[,vars11[i]]
                               +dat05[,vars12[i]]+dat05[,vars13[i]]+dat05[,vars14[i]]
                               +dat05[,vars15[i]]+dat05[,vars16[i]]+dat05[,vars17[i]],
                               data=dat05, cluster=dat05[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod05)[2,1])
  se<-as.numeric(summary(mod05)[2,2])
  coef05<-risk_se_ic_guapa(estimate,se)
  pval05<-pval_guapa(as.numeric(summary(mod05)[2,4]))
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02,coef03,pval03,coef04,pval04,coef05,pval05))
}

colnames(tab)<-c("Q1 (OR)","Q1 (pval)","Q2 (OR)","Q2 (pval)","Q3 (OR)","Q3 (pval)","Q4 (OR)","Q4 (pval)","Q5 (OR)","Q5 (pval)")
rownames(tab)<-vars20
write.table(tab,file="./results/mr_strata_q5.csv",sep=";",col.names=NA)


######################################################
### STEIGER FILTERING: CARDIOVASCULAR RISK FACTORS ###
######################################################

library(devtools)
library(googleAuthR)
library(MendelianRandomization)
library(mr.raps)
library(meta)
library(MRPRESSO)
library(MRInstruments)
library(MRMix)
library(RadialMR)
library(ieugwasr)
library(TwoSampleMR)


### CAD ###

cad<-read.csv2("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/cad/gwas_cad_varderharst2018.csv",
              header=TRUE,sep=";",dec=".")
cad$Phenotype<-c("Coronary artery disease")
cad$units<-c("units")
cad$id<-c("van der Harst P, 2018")
cad$samplesize<-547261
cad<-rename.vars(cad,
                from=c("chr","rsid"),
                to=c("chrom","SNP"))
cad<-cad[,c("Phenotype","SNP","beta","se","eaf","effect_allele","other_allele","pval",
          "units","samplesize","id")]
write.table(cad,"N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad_summary.txt",sep="\t",row.names=FALSE, quote=FALSE)


### BMI ###

bmi<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_CVRF/2smr/bmi.txt",
                              header=TRUE,sep="\t"))
names(bmi)<-tolower(names(bmi))
bmi<-rename.vars(bmi,
                 from=c("snp","tested_allele","p","n","freq_tested_allele_in_hrs"),
                 to=c("SNP","effect_allele","pval","samplesize","eaf"))
bmi<-bmi[,c("SNP","effect_allele","other_allele","beta","se","eaf","pval","samplesize")]
bmi$Units<-c("Units")
bmi$Phenotype<-c("BMI")
write.table(bmi,"N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/bmi_summary.txt",sep="\t",row.names=FALSE, quote=FALSE)


### FG ###

fg<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/fg_lagou2021.txt",
                             header=TRUE,sep="\t"))
names(fg)<-tolower(names(fg))
fg<-rename.vars(fg,
                from=c("rsid","a2","a1","p.value","n","maf"),
                to=c("SNP","effect_allele","other_allele","pval","samplesize","eaf"))
fg<-fg[,c("SNP","effect_allele","other_allele","beta","se","eaf","pval","samplesize")]
fg$Units<-c("Units")
fg$Phenotype<-c("Fasting glucose")
write.table(fg,"N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/fg_summary.txt",sep="\t",row.names=FALSE, quote=FALSE)


### FINS ###

fins<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/fins_lagou2021.txt",
                               header=TRUE,sep="\t"))
names(fins)<-tolower(names(fins))
fins<-rename.vars(fins,
                  from=c("rsid","a2","a1","p.value","n","maf"),
                  to=c("SNP","effect_allele","other_allele","pval","samplesize","eaf"))
fins<-fins[,c("SNP","effect_allele","other_allele","beta","se","eaf","pval","samplesize")]
fins$Units<-c("Units")
fins$Phenotype<-c("Fasting insulin")
write.table(fins,"N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/fins_summary.txt",sep="\t",row.names=FALSE, quote=FALSE)


### SBP ###

sbp<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/sbp_evangelou2018.txt",
                              header=TRUE,sep="\t"))
names(sbp)<-tolower(names(sbp))
sbp<-rename.vars(sbp,
                 from=c("id","alt","ref","n","alt_freq"),
                 to=c("SNP","effect_allele","other_allele","samplesize","eaf"))
sbp<-sbp[,c("SNP","effect_allele","other_allele","beta","se","eaf","pval","samplesize")]
sbp$Units<-c("Units")
sbp$Phenotype<-c("Systolic blood pressure")
write.table(sbp,"N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/sbp_summary.txt",sep="\t",row.names=FALSE, quote=FALSE)


### HDL-C ###

hdlc<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/hdlc_graham2021.results",
                               header=TRUE,sep="\t"))
names(hdlc)<-tolower(names(hdlc))
hdlc<-rename.vars(hdlc,
                  from=c("rsid","alt","ref","n","pooled_alt_af","effect_size","pvalue"),
                  to=c("SNP","effect_allele","other_allele","samplesize","eaf","beta","pval"))
hdlc<-hdlc[,c("SNP","effect_allele","other_allele","beta","se","eaf","pval","samplesize")]
hdlc$Units<-c("Units")
hdlc$Phenotype<-c("HDL cholesterol")
write.table(hdlc,"N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/hdlc_summary.txt",sep="\t",row.names=FALSE, quote=FALSE)


### LDL-C ###

ldlc<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/ldlc_graham2021.results",
                               header=TRUE,sep="\t"))
names(ldlc)<-tolower(names(ldlc))
ldlc<-rename.vars(ldlc,
                  from=c("rsid","alt","ref","n","pooled_alt_af","effect_size","pvalue"),
                  to=c("SNP","effect_allele","other_allele","samplesize","eaf","beta","pval"))
ldlc<-ldlc[,c("SNP","effect_allele","other_allele","beta","se","eaf","pval","samplesize")]
ldlc$Units<-c("Units")
ldlc$Phenotype<-c("HDL cholesterol")
write.table(ldlc,"N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/ldlc_summary.txt",sep="\t",row.names=FALSE, quote=FALSE)


### TG ###

tg<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/tg_graham2021.results",
                               header=TRUE,sep="\t"))
names(tg)<-tolower(names(tg))
tg<-rename.vars(tg,
                  from=c("rsid","alt","ref","n","pooled_alt_af","effect_size","pvalue"),
                  to=c("SNP","effect_allele","other_allele","samplesize","eaf","beta","pval"))
tg<-tg[,c("SNP","effect_allele","other_allele","beta","se","eaf","pval","samplesize")]
tg$Units<-c("Units")
tg$Phenotype<-c("HDL cholesterol")
write.table(tg,"N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/tg_summary.txt",sep="\t",row.names=FALSE, quote=FALSE)


### STEIGER FILTERING OF ALL TRAITS ###
#######################################

cad_harm<-read_exposure_data(
  filename = "N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad_summary.txt",
  clump = FALSE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE
)


### BMI ###

bmi_harm <- read_outcome_data(
  snps = cad_harm$SNP,
  filename = "N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/bmi_summary.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

cad_bmi<-harmonise_data(cad_harm, bmi_harm, action = 1)
cad_bmi_steiger<-steiger_filtering(cad_bmi)
cad_bmi_steiger<-subset2(cad_bmi_steiger,"cad_bmi_steiger$steiger_dir==TRUE")
dim(cad_bmi)[1]-dim(cad_bmi_steiger)[1]
save(cad_bmi_steiger,file="N:/data/durable/Projects/Magnus_MR_CVRF/Outputs/steiger/source/cad_bmi_steiger.RData")

# 3 SNPs EXPLAIN MORE VARIATION IN BMI THAN IN CAD #

bmi_steiger_excl<-setdiff(cad_bmi$SNP,cad_bmi_steiger$SNP)
bmi_steiger_ok<-intersect(cad_bmi$SNP,cad_bmi_steiger$SNP)


### FG ###

fg_harm <- read_outcome_data(
  snps = cad_harm$SNP,
  filename = "N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/fg_summary.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

cad_fg<-harmonise_data(cad_harm, fg_harm, action = 1)
cad_fg_steiger<-steiger_filtering(cad_fg)
cad_fg_steiger<-subset2(cad_fg_steiger,"cad_fg_steiger$steiger_dir==TRUE")
dim(cad_fg)[1]-dim(cad_fg_steiger)[1]
save(cad_fg_steiger,file="N:/data/durable/Projects/Magnus_MR_CVRF/Outputs/steiger/source/cad_fg_steiger.RData")

# 3 SNPs EXPLAIN MORE VARIATION IN fg THAN IN CAD #

fg_steiger_excl<-setdiff(cad_fg$SNP,cad_fg_steiger$SNP)
fg_steiger_ok<-intersect(cad_fg$SNP,cad_fg_steiger$SNP)


### FINS ###

fins_harm <- read_outcome_data(
  snps = cad_harm$SNP,
  filename = "N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/fins_summary.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

cad_fins<-harmonise_data(cad_harm, fins_harm, action = 1)
cad_fins_steiger<-steiger_filtering(cad_fins)
cad_fins_steiger<-subset2(cad_fins_steiger,"cad_fins_steiger$steiger_dir==TRUE")
dim(cad_fins)[1]-dim(cad_fins_steiger)[1]
save(cad_fins_steiger,file="N:/data/durable/Projects/Magnus_MR_CVRF/Outputs/steiger/source/cad_fins_steiger.RData")

# 7 SNPs EXPLAIN MORE VARIATION IN fins THAN IN CAD #

fins_steiger_excl<-setdiff(cad_fins$SNP,cad_fins_steiger$SNP)
fins_steiger_ok<-intersect(cad_fins$SNP,cad_fins_steiger$SNP)


### SBP ###

sbp_harm <- read_outcome_data(
  snps = cad_harm$SNP,
  filename = "N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/sbp_summary.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

cad_sbp<-harmonise_data(cad_harm, sbp_harm, action = 1)
cad_sbp_steiger<-steiger_filtering(cad_sbp)
cad_sbp_steiger<-subset2(cad_sbp_steiger,"cad_sbp_steiger$steiger_dir==TRUE")
dim(cad_sbp)[1]-dim(cad_sbp_steiger)[1]
save(cad_sbp_steiger,file="N:/data/durable/Projects/Magnus_MR_CVRF/Outputs/steiger/source/cad_sbp_steiger.RData")

# 7 SNPs EXPLAIN MORE VARIATION IN sbp THAN IN CAD #

sbp_steiger_excl<-setdiff(cad_sbp$SNP,cad_sbp_steiger$SNP)
sbp_steiger_ok<-intersect(cad_sbp$SNP,cad_sbp_steiger$SNP)


### HDL-C ###

hdlc_harm <- read_outcome_data(
  snps = cad_harm$SNP,
  filename = "N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/hdlc_summary.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

cad_hdlc<-harmonise_data(cad_harm, hdlc_harm, action = 1)
cad_hdlc_steiger<-steiger_filtering(cad_hdlc)
cad_hdlc_steiger<-subset2(cad_hdlc_steiger,"cad_hdlc_steiger$steiger_dir==TRUE")
dim(cad_hdlc)[1]-dim(cad_hdlc_steiger)[1]
save(cad_hdlc_steiger,file="N:/data/durable/Projects/Magnus_MR_CVRF/Outputs/steiger/source/cad_hdlc_steiger.RData")

# 15 SNPs EXPLAIN MORE VARIATION IN hdlc THAN IN CAD #

hdlc_steiger_excl<-setdiff(cad_hdlc$SNP,cad_hdlc_steiger$SNP)
hdlc_steiger_ok<-intersect(cad_hdlc$SNP,cad_hdlc_steiger$SNP)


### LDL-C ###

ldlc_harm <- read_outcome_data(
  snps = cad_harm$SNP,
  filename = "N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/ldlc_summary.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

cad_ldlc<-harmonise_data(cad_harm, ldlc_harm, action = 2)
cad_ldlc_steiger<-steiger_filtering(cad_ldlc)
cad_ldlc_steiger<-subset2(cad_ldlc_steiger,"cad_ldlc_steiger$steiger_dir==TRUE")
dim(cad_ldlc)[1]-dim(cad_ldlc_steiger)[1]
save(cad_ldlc_steiger,file="N:/data/durable/Projects/Magnus_MR_CVRF/Outputs/steiger/source/cad_ldlc_steigerr.RData")

# 20 SNPs EXPLAIN MORE VARIATION IN ldlc THAN IN CAD #

ldlc_steiger_excl<-setdiff(cad_ldlc$SNP,cad_ldlc_steiger$SNP)
ldlc_steiger_ok<-intersect(cad_ldlc$SNP,cad_ldlc_steiger$SNP)


### TG ###

tg_harm <- read_outcome_data(
  snps = cad_harm$SNP,
  filename = "N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/tg_summary.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

cad_tg<-harmonise_data(cad_harm, tg_harm, action = 2)
cad_tg_steiger<-steiger_filtering(cad_tg)
cad_tg_steiger<-subset2(cad_tg_steiger,"cad_tg_steiger$steiger_dir==TRUE")
dim(cad_tg)[1]-dim(cad_tg_steiger)[1]
save(cad_tg_steiger,file="N:/data/durable/Projects/Magnus_MR_CVRF/Outputs/steiger/source/cad_tg_steigerr.RData")

# 13 SNPs EXPLAIN MORE VARIATION IN tg THAN IN CAD #

tg_steiger_excl<-setdiff(cad_tg$SNP,cad_tg_steiger$SNP)
tg_steiger_ok<-intersect(cad_tg$SNP,cad_tg_steiger$SNP)


### SNPs TO BE EXCLUDED ###

snp_steiger<-union(bmi_steiger_excl,fg_steiger_excl)
snp_steiger<-union(snp_steiger,fins_steiger_excl)
snp_steiger<-union(snp_steiger,sbp_steiger_excl)
snp_steiger<-union(snp_steiger,hdlc_steiger_excl)
snp_steiger<-union(snp_steiger,ldlc_steiger_excl)
snp_steiger<-union(snp_steiger,tg_steiger_excl)

indep_cad_snps<-setdiff(cad$SNP,snp_steiger)
write.table(indep_cad_snps,"N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/indep_cad_snps.txt",sep="\t",row.names=FALSE, quote=FALSE)


################################
### EXTRACT SNPs USING PLINK ###
################################


############################
### GENERATION OF DM GRS ###
############################

# MoBa SELECTED SNPs ON DM #

chr01<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr1.raw",header=TRUE,sep="\t"))
chr02<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr2.raw",header=TRUE,sep="\t"))
chr03<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr3.raw",header=TRUE,sep="\t"))
chr04<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr4.raw",header=TRUE,sep="\t"))
chr05<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr5.raw",header=TRUE,sep="\t"))
chr06<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr6.raw",header=TRUE,sep="\t"))
chr07<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr7.raw",header=TRUE,sep="\t"))
chr08<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr8.raw",header=TRUE,sep="\t"))
chr09<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr9.raw",header=TRUE,sep="\t"))
chr10<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr10.raw",header=TRUE,sep="\t"))
chr11<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr11.raw",header=TRUE,sep="\t"))
chr12<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr12.raw",header=TRUE,sep="\t"))
chr13<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr13.raw",header=TRUE,sep="\t"))
chr14<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr14.raw",header=TRUE,sep="\t"))
chr15<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr15.raw",header=TRUE,sep="\t"))
chr16<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr16.raw",header=TRUE,sep="\t"))
chr17<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr17.raw",header=TRUE,sep="\t"))
chr18<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr18.raw",header=TRUE,sep="\t"))
chr19<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr19.raw",header=TRUE,sep="\t"))
chr20<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr20.raw",header=TRUE,sep="\t"))
chr21<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr21.raw",header=TRUE,sep="\t"))
chr22<-as.data.frame(read.delim("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/chr22.raw",header=TRUE,sep="\t"))

chr01<-chr01[,c(1,7:dim(chr01)[2])]
chr02<-chr02[,c(1,7:dim(chr02)[2])]
chr03<-chr03[,c(1,7:dim(chr03)[2])]
chr04<-chr04[,c(1,7:dim(chr04)[2])]
chr05<-chr05[,c(1,7:dim(chr05)[2])]
chr06<-chr06[,c(1,7:dim(chr06)[2])]
chr07<-chr07[,c(1,7:dim(chr07)[2])]
chr08<-chr08[,c(1,7:dim(chr08)[2])]
chr09<-chr09[,c(1,7:dim(chr09)[2])]
chr10<-chr10[,c(1,7:dim(chr10)[2])]
chr11<-chr11[,c(1,7:dim(chr11)[2])]
chr12<-chr12[,c(1,7:dim(chr12)[2])]
chr13<-chr13[,c(1,7:dim(chr13)[2])]
chr14<-chr14[,c(1,7:dim(chr14)[2])]
chr15<-chr15[,c(1,7:dim(chr15)[2])]
chr16<-chr16[,c(1,7:dim(chr16)[2])]
chr17<-chr17[,c(1,7:dim(chr17)[2])]
#chr18<-chr18[,c(1,7:dim(chr18)[2])]
chr19<-chr19[,c(1,7:dim(chr19)[2])]
chr20<-chr20[,c(1,7:dim(chr20)[2])]
chr21<-chr21[,c(1,7:dim(chr21)[2])]
#chr22<-chr22[,c(1,7:dim(chr22)[2])]

dat<-merge2(chr01,chr02,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr03,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr04,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr05,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr06,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr07,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr08,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr09,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr10,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr11,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr12,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr13,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr14,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr15,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr16,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr17,by.id=c("FID"),all.x=TRUE,sort=FALSE)
#dat<-merge2(dat,chr18,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr19,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr20,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr21,by.id=c("FID"),all.x=TRUE,sort=FALSE)
#dat<-merge2(dat,chr22,by.id=c("FID"),all.x=TRUE,sort=FALSE)


# SNPs RELATED TO CAD #

cad_ma<-read.csv2("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/steiger/cad/gwas_cad_varderharst2018.csv",header=TRUE,
                  sep=";",dec=".")

names(cad_ma)<-tolower(names(cad_ma))
cad_ma<-rename.vars(cad_ma,
                    from=c("effect_allele","eaf"),
                    to=c("tested_allele_ma","maf_ma"))
cad_ma$tested_allele_ma<-tolower(cad_ma$tested_allele_ma)


### CHECK THAT THE 102 SNPs IN DAT BELONG TO THE META-ANALYSIS RESULTS ###

bbb<-as.character(sort(cad_ma$rsid,decreasing=FALSE))

dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))

length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND cad_ma MATCH


# TABLE WITH NUMBER OF MISSING VALUES OF ALL VARIABLES IN dat #

vars01<-names(dat)
tab_na<-NULL

for(i in 1:length(vars01))
  
{
  x<-length(which(is.na(dat[,vars01[i]])))
  tab_na<-as.data.frame(rbind(tab_na,cbind(x)))
}

colnames(tab_na)<-c("NAs")
rownames(tab_na)<-vars01
table(tab_na$NAs)


# EXTRACT THE INFORMATION ON THE ALLELE TESTED IN MOBA #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
cad_ma<-merge2(cad_ma,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
cad_ma<-na.omit(cad_ma)


# SNP FLIPPING IF THE ALLELE ORDER IN MoBa AND THE META-ANALYSIS DIFFER #

cad_ma$same_coding<-with(cad_ma,ifelse(tested_allele_moba==tested_allele_ma,1,
                                       ifelse(tested_allele_moba!=tested_allele_ma,0,NA)))

cad_ma$flip<-with(cad_ma,ifelse(beta>=0 & same_coding==1,0,
                                ifelse(beta<0 & same_coding==0,0,
                                       ifelse(beta>=0 & same_coding==0,1,
                                              ifelse(beta<0 & same_coding==1,1,NA)))))

cad_ma$maf<-with(cad_ma,ifelse(same_coding==1,maf_ma,
                               ifelse(same_coding==0,1-maf_ma,NA)))
cad_ma$coef<-with(cad_ma,ifelse(beta>0,beta,
                                ifelse(beta<0,-beta,NA)))


### CALCULATION OF CAD GRS ###

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("FID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-cad_ma[cad_ma$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-cad_ma[cad_ma$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$cad_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","cad_grs")]
dat<-rename.vars(dat,from=c("cad_grs"),to=c("cad_steiger_grs"))
save(dat,file="N:/data/durable/Projects/Hernaez_MR_CVD/source_files/cad_steiger_grs.RData")

cad_steiger_grs<-dat
load("N:/data/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")

cad_steiger_grs<-rename.vars(cad_steiger_grs,from=c("id","cad_steiger_grs"),to=c("sentrixid_mom","cad_steiger_grs_mom"))
dat<-merge2(dat,cad_steiger_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
cad_steiger_grs<-rename.vars(cad_steiger_grs,from=c("sentrixid_mom","cad_steiger_grs_mom"),to=c("sentrixid_dad","cad_steiger_grs_dad"))
dat<-merge2(dat,cad_steiger_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

save(dat,file="N:/data/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")


#####################
### MAIN ANALYSES ###
#####################

setwd("N:/data/durable/Projects/Hernaez_MR_CVD")

dir.create("./Outputs/sensitivity")
dir.create("./Outputs/sensitivity/steiger")

setwd("N:/data/durable/Projects/Hernaez_MR_CVD/Outputs/sensitivity/steiger")


### LOGISTIC REGRESSION / MENDELIAN RANDOMIZATION: LINEAR ASSOCIATIONS ###
##########################################################################

load("N:/data/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")

vars00<-c("cad_steiger_grs_mom","cad_steiger_grs_dad","cad_steiger_grs_mom","cad_steiger_grs_dad","cad_steiger_grs_mom","cad_steiger_grs_dad","cad_steiger_grs_mom","cad_steiger_grs_dad","cad_steiger_grs_mom","cad_steiger_grs_dad")
vars01<-c("cad_steiger_grs_mom_z","cad_steiger_grs_dad_z","cad_steiger_grs_mom_z","cad_steiger_grs_dad_z","cad_steiger_grs_mom_z","cad_steiger_grs_dad_z","cad_steiger_grs_mom_z","cad_steiger_grs_dad_z","cad_steiger_grs_mom_z","cad_steiger_grs_dad_z")
vars02<-c("subfertility_12plus","subfertility_12plus","miscarriage","miscarriage","eclampsia","eclampsia","preterm","preterm","sga","sga")
vars08<-c("pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad")
vars09<-c("pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad")
vars10<-c("pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad")
vars11<-c("pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad")
vars12<-c("pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad")
vars13<-c("pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad")
vars14<-c("pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad")
vars15<-c("pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad")
vars16<-c("pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad")
vars17<-c("pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad")
vars18<-c("m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374")
vars20<-c("CAD-GRS - Subfertility, mothers","CAD-GRS - Subfertility, fathers",
          "CAD-GRS - Miscarriage, mothers","CAD-GRS - Miscarriage, fathers",
          "CAD-GRS - Pre-/Eclampsia, mothers","CAD-GRS - Pre-/Eclampsia, fathers",
          "CAD-GRS - Preterm, mothers","CAD-GRS - Preterm, fathers",
          "CAD-GRS - SGA, mothers","CAD-GRS - SGA, fathers")

tab<-NULL
for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"!is.na(dat[,vars00[i]]) & dat$preg_plan==1")
  datx[,vars01[i]]<-as.numeric(with(datx,scale(datx[,vars00[i]])))
  
  mod03<-miceadds::glm.cluster(formula=as.factor(datx[,vars02[i]])~datx[,vars01[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod03)[2,1])
  se<-as.numeric(summary(mod03)[2,2])
  coef03<-risk_se_ic_guapa(estimate,se)
  pval03<-pval_guapa(as.numeric(summary(mod03)[2,4]))
  
  mod04<-miceadds::glm.cluster(formula=as.factor(datx[,vars02[i]])~datx[,vars01[i]]
                               +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
                               +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
                               +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod04)[2,1])
  se<-as.numeric(summary(mod04)[2,2])
  coef04<-risk_se_ic_guapa(estimate,se)
  pval04<-pval_guapa(as.numeric(summary(mod04)[2,4]))
  forest_mr<-risk_se_ic_guapa2(estimate,se)
  
  aaa<-datx[,vars01[i]]
  mod_base<-gam(formula=as.factor(datx[,vars02[i]])~datx[,vars08[i]]
                +datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
                +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
                +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]],
                data=datx, family="binomial")
  mod_lin<-gam(formula=as.factor(datx[,vars02[i]])~aaa
               +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
               +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
               +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]],
               data=datx, family="binomial")
  mod_nlin<-gam(formula=as.factor(datx[,vars02[i]])~bs(aaa,df=4)
                +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
                +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
                +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]],
                data=datx, family="binomial")
  p_lin<-pval_guapa(lrtest(mod_base,mod_lin)[2,5])
  p_nonlin<-pval_guapa(lrtest(mod_base,mod_nlin)[2,5])
  p_lrtest<-pval_guapa(lrtest(mod_lin,mod_nlin)[2,5])
  
  tab<-rbind(tab,cbind(coef03,pval03,coef04,pval04,p_lin,p_nonlin,p_lrtest,forest_mr))
}

colnames(tab)<-c("OR (raw)","pval (raw)","OR (adj.)","pval (adj.)",
                 "p_lin","p_nonlin","p_lrtest","forest")
rownames(tab)<-vars20
write.table(tab,file="./linear.csv",sep=";",col.names=NA)


### LOGISTIC REGRESSION: NON-LINEAR ASSOCIATIONS (SMOOTHED SPLINES) ###
#######################################################################

# EXECUTE MANUALLY: i=1, then the loop #

load("N:/data/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")

vars00<-c("cad_steiger_grs_mom","cad_steiger_grs_dad","cad_steiger_grs_mom","cad_steiger_grs_dad","cad_steiger_grs_mom","cad_steiger_grs_dad","cad_steiger_grs_mom","cad_steiger_grs_dad","cad_steiger_grs_mom","cad_steiger_grs_dad")
vars01<-c("cad_steiger_grs_mom_z","cad_steiger_grs_dad_z","cad_steiger_grs_mom_z","cad_steiger_grs_dad_z","cad_steiger_grs_mom_z","cad_steiger_grs_dad_z","cad_steiger_grs_mom_z","cad_steiger_grs_dad_z","cad_steiger_grs_mom_z","cad_steiger_grs_dad_z")
vars02<-c("subfertility_12plus","subfertility_12plus","miscarriage","miscarriage","eclampsia","eclampsia","preterm","preterm","sga","sga")
vars08<-c("pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad")
vars09<-c("pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad")
vars10<-c("pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad")
vars11<-c("pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad")
vars12<-c("pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad")
vars13<-c("pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad")
vars14<-c("pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad")
vars15<-c("pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad")
vars16<-c("pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad")
vars17<-c("pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad")
vars18<-c("m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374")
vars20<-c("CAD-GRS, mothers","CAD-GRS, fathers",
          "CAD-GRS, mothers","CAD-GRS, fathers",
          "CAD-GRS, mothers","CAD-GRS, fathers",
          "CAD-GRS, mothers","CAD-GRS, fathers",
          "CAD-GRS, mothers","CAD-GRS, fathers")
vars21<-c("Subfertility (odds ratio, adjusted)","Subfertility (odds ratio, adjusted)",
          "Miscarriage (odds ratio, adjusted)","Miscarriage (odds ratio, adjusted)",
          "Pre-/Eclampsia (odds ratio, adjusted)","Pre-/Eclampsia (odds ratio, adjusted)",
          "Preterm birth (odds ratio, adjusted)","Preterm birth (odds ratio, adjusted)",
          "SGA (odds ratio, adjusted)","SGA (odds ratio, adjusted)")

tab<-NULL
for(i in 1:length(vars01))
  
{
  varstot<-c(vars00[i],vars02[i],vars08[i],vars09[i],vars10[i],vars11[i],vars12[i],vars13[i],
             vars14[i],vars15[i],vars16[i],vars17[i],vars18[i],"preg_plan")
  dat2<-na.omit(dat[,varstot])
  dat2<-subset2(dat2,"dat2$preg_plan==1")
  dat2[,vars01[i]]<-as.numeric(with(dat2,scale(dat2[,vars00[i]])))
  aaa<-dat2[,vars01[i]]
  minim<-mean(aaa,na.rm=TRUE)
  
  mod01<-gam(formula=as.factor(dat2[,vars02[i]])~bs(aaa,df=4)
             +dat2[,vars08[i]]+dat2[,vars09[i]]+dat2[,vars10[i]]+dat2[,vars11[i]]
             +dat2[,vars12[i]]+dat2[,vars13[i]]+dat2[,vars14[i]]
             +dat2[,vars15[i]]+dat2[,vars16[i]]+dat2[,vars17[i]],
             data=dat2, family="binomial")
  ptemp<-termplot(mod01,term=1,se=TRUE,plot=FALSE)
  temp<-ptemp$aaa
  value<-closest(temp$x,minim)
  #min_val<-temp$x[which(temp$y==min(temp$y,na.rm=TRUE))]
  center<-with(temp, y[x==value])
  z<-qnorm(1-0.05/2)
  ytemp<-temp$y+outer(temp$se,c(0,-z,z),'*')
  min_val<-guapa(temp$x[which(temp$y==min(temp$y,na.rm=TRUE))])
  ci<-exp(ytemp-center)
  name<-paste("./",vars01[i],"_",vars02[i],".jpg",sep="")
  labely<-vars21[i]
  
  plot.data<-as.data.frame(cbind(temp$x,ci))
  colnames(plot.data)<-c("xval","yest","lci","uci")
  ft<-as.data.frame(table(dat2[,vars01[i]]))
  names(ft)<-c("xval","freq")
  plot.data<-merge2(plot.data,ft,by.id=c("xval"),all.x=TRUE,sort=FALSE)
  infl01<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[1],]$xval
  infl02<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[2],]$xval
  infl03<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[3],]$xval
  infl04<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[4],]$xval
  infl05<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[5],]$xval
  infl06<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[6],]$xval
  infl07<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[7],]$xval
  infl08<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[8],]$xval
  infl09<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[9],]$xval
  infl10<-plot.data[which(c(FALSE, diff(diff(plot.data$yest)>0)!=0)==TRUE)[10],]$xval
  
  figure<-ggplot() +
    geom_histogram(data=dat2, aes(x=dat2[,vars01[i]], y=(..density../max(..density..,na.rm=TRUE))*3.75), 
                   bins=20, color="grey80", fill="grey90") +
    #scale_x_continuous(limits = c(min(dat2[,vars01[i]],na.rm=TRUE),max(dat2[,vars01[i]],na.rm=TRUE))) +
    scale_y_continuous(limits = c(0,4), sec.axis = sec_axis(~., name = "Relative density")) +
    geom_hline(aes(yintercept=1), data=plot.data, colour="black", linetype=2) + 
    geom_line(aes(x=xval, y=yest), data=plot.data, color="black") + 
    geom_line(aes(x=xval, y=lci), data=plot.data, color="grey35") + 
    geom_line(aes(x=xval ,y=uci), data=plot.data, color="grey35") + 
    theme_bw() +
    labs(x=vars20[i],y=labely) +
    theme(axis.title.x = element_text(vjust=0.5, size=20), 
          axis.title.y = element_text(vjust=0.5, size=20),
          axis.text.x = element_text(size=16, colour = 'black'),
          axis.text.y = element_text(size=16, colour = 'black'),
          axis.ticks.x = element_line(colour = 'black'),
          axis.ticks.y = element_line(colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank())  
  #    geom_point(aes(x=infl01, y=with(plot.data, yest[xval==infl01])), data=plot.data, colour="black", size=2.5) +
  #    geom_point(aes(x=infl02, y=with(plot.data, yest[xval==infl02])), data=plot.data, colour="black", size=2.5)
  
  jpeg(filename=name,width = 8000, height = 8000, res=1200)
  par(las=1,cex=1.2,mar=c(6,6,2,0),bty="n",lheight=0.9)
  figure
  dev.off()
  
  tab<-rbind(tab,cbind(infl01,infl02,infl03,infl04,infl05,infl06,infl07,infl08,infl09,infl10))
}

rownames(tab)<-vars01
write.table(tab,file="./inflection_points.csv",sep=";",col.names=NA)


### MENDELIAN RANDOMIZATION IN STRATA ###
#########################################

load("N:/data/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")

vars00<-c("cad_steiger_grs_mom","cad_steiger_grs_dad","cad_steiger_grs_mom","cad_steiger_grs_dad","cad_steiger_grs_mom","cad_steiger_grs_dad","cad_steiger_grs_mom","cad_steiger_grs_dad","cad_steiger_grs_mom","cad_steiger_grs_dad")
vars01<-c("cad_steiger_grs_mom_q","cad_steiger_grs_dad_q","cad_steiger_grs_mom_q","cad_steiger_grs_dad_q","cad_steiger_grs_mom_q","cad_steiger_grs_dad_q","cad_steiger_grs_mom_q","cad_steiger_grs_dad_q","cad_steiger_grs_mom_q","cad_steiger_grs_dad_q")
vars02<-c("cad_steiger_grs_mom_z","cad_steiger_grs_dad_z","cad_steiger_grs_mom_z","cad_steiger_grs_dad_z","cad_steiger_grs_mom_z","cad_steiger_grs_dad_z","cad_steiger_grs_mom_z","cad_steiger_grs_dad_z","cad_steiger_grs_mom_z","cad_steiger_grs_dad_z")
vars03<-c("subfertility_12plus","subfertility_12plus","miscarriage","miscarriage","eclampsia","eclampsia","preterm","preterm","sga","sga")
vars08<-c("pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad","pc1_mom","pc1_dad")
vars09<-c("pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad","pc2_mom","pc2_dad")
vars10<-c("pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad","pc3_mom","pc3_dad")
vars11<-c("pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad","pc4_mom","pc4_dad")
vars12<-c("pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad","pc5_mom","pc5_dad")
vars13<-c("pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad","pc6_mom","pc6_dad")
vars14<-c("pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad","pc7_mom","pc7_dad")
vars15<-c("pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad","pc8_mom","pc8_dad")
vars16<-c("pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad","pc9_mom","pc9_dad")
vars17<-c("pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad")
vars18<-c("m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374","m_id_2374","f_id_2374")
vars20<-c("CAD-GRS - Subfertility, mothers","CAD-GRS - Subfertility, fathers",
          "CAD-GRS - Miscarriage, mothers","CAD-GRS - Miscarriage, fathers",
          "CAD-GRS - Pre-/Eclampsia, mothers","CAD-GRS - Pre-/Eclampsia, fathers",
          "CAD-GRS - Preterm, mothers","CAD-GRS - Preterm, fathers",
          "CAD-GRS - SGA, mothers","CAD-GRS - SGA, fathers")

tab<-NULL
for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"!is.na(dat[,vars00[i]]) & dat$preg_plan==1")
  datx[,vars01[i]]<-as.numeric(ntile(datx[,vars00[i]], 5))
  datx[,vars02[i]]<-as.numeric(with(datx,scale(datx[,vars00[i]])))
  
  dat01<-subset2(datx,"datx[,vars01[i]]==1")
  dat02<-subset2(datx,"datx[,vars01[i]]==2")
  dat03<-subset2(datx,"datx[,vars01[i]]==3")
  dat04<-subset2(datx,"datx[,vars01[i]]==4")
  dat05<-subset2(datx,"datx[,vars01[i]]==5")
  
  mod01<-miceadds::glm.cluster(formula=as.factor(dat01[,vars03[i]])~dat01[,vars02[i]]
                               +dat01[,vars08[i]]+dat01[,vars09[i]]+dat01[,vars10[i]]+dat01[,vars11[i]]
                               +dat01[,vars12[i]]+dat01[,vars13[i]]+dat01[,vars14[i]]
                               +dat01[,vars15[i]]+dat01[,vars16[i]]+dat01[,vars17[i]],
                               data=dat01, cluster=dat01[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod01)[2,1])
  se<-as.numeric(summary(mod01)[2,2])
  coef01<-risk_se_ic_guapa(estimate,se)
  pval01<-pval_guapa(as.numeric(summary(mod01)[2,4]))
  
  mod02<-miceadds::glm.cluster(formula=as.factor(dat02[,vars03[i]])~dat02[,vars02[i]]
                               +dat02[,vars08[i]]+dat02[,vars09[i]]+dat02[,vars10[i]]+dat02[,vars11[i]]
                               +dat02[,vars12[i]]+dat02[,vars13[i]]+dat02[,vars14[i]]
                               +dat02[,vars15[i]]+dat02[,vars16[i]]+dat02[,vars17[i]],
                               data=dat02, cluster=dat02[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod02)[2,1])
  se<-as.numeric(summary(mod02)[2,2])
  coef02<-risk_se_ic_guapa(estimate,se)
  pval02<-pval_guapa(as.numeric(summary(mod02)[2,4]))
  
  mod03<-miceadds::glm.cluster(formula=as.factor(dat03[,vars03[i]])~dat03[,vars02[i]]
                               +dat03[,vars08[i]]+dat03[,vars09[i]]+dat03[,vars10[i]]+dat03[,vars11[i]]
                               +dat03[,vars12[i]]+dat03[,vars13[i]]+dat03[,vars14[i]]
                               +dat03[,vars15[i]]+dat03[,vars16[i]]+dat03[,vars17[i]],
                               data=dat03, cluster=dat03[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod03)[2,1])
  se<-as.numeric(summary(mod03)[2,2])
  coef03<-risk_se_ic_guapa(estimate,se)
  pval03<-pval_guapa(as.numeric(summary(mod03)[2,4]))
  
  mod04<-miceadds::glm.cluster(formula=as.factor(dat04[,vars03[i]])~dat04[,vars02[i]]
                               +dat04[,vars08[i]]+dat04[,vars09[i]]+dat04[,vars10[i]]+dat04[,vars11[i]]
                               +dat04[,vars12[i]]+dat04[,vars13[i]]+dat04[,vars14[i]]
                               +dat04[,vars15[i]]+dat04[,vars16[i]]+dat04[,vars17[i]],
                               data=dat04, cluster=dat04[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod04)[2,1])
  se<-as.numeric(summary(mod04)[2,2])
  coef04<-risk_se_ic_guapa(estimate,se)
  pval04<-pval_guapa(as.numeric(summary(mod04)[2,4]))
  
  mod05<-miceadds::glm.cluster(formula=as.factor(dat05[,vars03[i]])~dat05[,vars02[i]]
                               +dat05[,vars08[i]]+dat05[,vars09[i]]+dat05[,vars10[i]]+dat05[,vars11[i]]
                               +dat05[,vars12[i]]+dat05[,vars13[i]]+dat05[,vars14[i]]
                               +dat05[,vars15[i]]+dat05[,vars16[i]]+dat05[,vars17[i]],
                               data=dat05, cluster=dat05[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod05)[2,1])
  se<-as.numeric(summary(mod05)[2,2])
  coef05<-risk_se_ic_guapa(estimate,se)
  pval05<-pval_guapa(as.numeric(summary(mod05)[2,4]))
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02,coef03,pval03,coef04,pval04,coef05,pval05))
}

colnames(tab)<-c("Q1 (OR)","Q1 (pval)","Q2 (OR)","Q2 (pval)","Q3 (OR)","Q3 (pval)","Q4 (OR)","Q4 (pval)","Q5 (OR)","Q5 (pval)")
rownames(tab)<-vars20
write.table(tab,file="./mr_strata_q5.csv",sep=";",col.names=NA)


##############################################################################################################

#######################################################################
### ADDITIONAL STEPS IN GWAS RUN BY CHROMOSOMES (OLD GENOTYPE DATA) ###
#######################################################################


### stillbirth ###
##################

### HARVEST - MOTHERS - stillbirth ###

setwd("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/preg_outc/results")
z<-qnorm(1-0.05/2)

h_m_chr_1<-as.data.frame(read.delim("./h_m_chr_1_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_2<-as.data.frame(read.delim("./h_m_chr_2_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_3<-as.data.frame(read.delim("./h_m_chr_3_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_4<-as.data.frame(read.delim("./h_m_chr_4_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_5<-as.data.frame(read.delim("./h_m_chr_5_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_6<-as.data.frame(read.delim("./h_m_chr_6_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_7<-as.data.frame(read.delim("./h_m_chr_7_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_8<-as.data.frame(read.delim("./h_m_chr_8_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_9<-as.data.frame(read.delim("./h_m_chr_9_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_10<-as.data.frame(read.delim("./h_m_chr_10_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_11<-as.data.frame(read.delim("./h_m_chr_11_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_12<-as.data.frame(read.delim("./h_m_chr_12_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_13<-as.data.frame(read.delim("./h_m_chr_13_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_14<-as.data.frame(read.delim("./h_m_chr_14_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_15<-as.data.frame(read.delim("./h_m_chr_15_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_16<-as.data.frame(read.delim("./h_m_chr_16_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_17<-as.data.frame(read.delim("./h_m_chr_17_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_18<-as.data.frame(read.delim("./h_m_chr_18_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_19<-as.data.frame(read.delim("./h_m_chr_19_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_20<-as.data.frame(read.delim("./h_m_chr_20_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_21<-as.data.frame(read.delim("./h_m_chr_21_stillbirth.txt",header=FALSE,sep="\t"))
h_m_chr_22<-as.data.frame(read.delim("./h_m_chr_22_stillbirth.txt",header=FALSE,sep="\t"))

harvest_mom<-rbind(h_m_chr_1,h_m_chr_2,h_m_chr_3,h_m_chr_4,h_m_chr_5,h_m_chr_6,h_m_chr_7,h_m_chr_8,
                   h_m_chr_9,h_m_chr_10,h_m_chr_11,h_m_chr_12,h_m_chr_13,h_m_chr_1,h_m_chr_15,
                   h_m_chr_16,h_m_chr_17,h_m_chr_18,h_m_chr_19,h_m_chr_20,h_m_chr_21,h_m_chr_22)
names(harvest_mom)<-c("CHROM","POS","MARKERNAME","NEA","EA","A1","TEST","N","OR","SE","Z_STAT","PVAL")
harvest_mom$BETA<-log(harvest_mom$OR)
harvest_mom$OR_SE<-harvest_mom$OR*harvest_mom$SE
harvest_mom$OR_95L<-harvest_mom$OR-(z*harvest_mom$OR_SE)
harvest_mom$OR_95U<-harvest_mom$OR+(z*harvest_mom$OR_SE)
write.table(harvest_mom,"./harvest_mom_stillbirth.txt",sep="\t",row.names=FALSE, quote=FALSE)


### ROTTERDAM1 - MOTHERS - stillbirth ###

r1_m_chr_1<-as.data.frame(read.delim("./r1_m_chr_1_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_2<-as.data.frame(read.delim("./r1_m_chr_2_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_3<-as.data.frame(read.delim("./r1_m_chr_3_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_4<-as.data.frame(read.delim("./r1_m_chr_4_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_5<-as.data.frame(read.delim("./r1_m_chr_5_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_6<-as.data.frame(read.delim("./r1_m_chr_6_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_7<-as.data.frame(read.delim("./r1_m_chr_7_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_8<-as.data.frame(read.delim("./r1_m_chr_8_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_9<-as.data.frame(read.delim("./r1_m_chr_9_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_10<-as.data.frame(read.delim("./r1_m_chr_10_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_11<-as.data.frame(read.delim("./r1_m_chr_11_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_12<-as.data.frame(read.delim("./r1_m_chr_12_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_13<-as.data.frame(read.delim("./r1_m_chr_13_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_14<-as.data.frame(read.delim("./r1_m_chr_14_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_15<-as.data.frame(read.delim("./r1_m_chr_15_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_16<-as.data.frame(read.delim("./r1_m_chr_16_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_17<-as.data.frame(read.delim("./r1_m_chr_17_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_18<-as.data.frame(read.delim("./r1_m_chr_18_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_19<-as.data.frame(read.delim("./r1_m_chr_19_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_20<-as.data.frame(read.delim("./r1_m_chr_20_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_21<-as.data.frame(read.delim("./r1_m_chr_21_stillbirth.txt",header=FALSE,sep="\t"))
r1_m_chr_22<-as.data.frame(read.delim("./r1_m_chr_22_stillbirth.txt",header=FALSE,sep="\t"))

rotterdam1_mom<-rbind(r1_m_chr_1,r1_m_chr_2,r1_m_chr_3,r1_m_chr_4,r1_m_chr_5,r1_m_chr_6,r1_m_chr_7,r1_m_chr_8,
                      r1_m_chr_9,r1_m_chr_10,r1_m_chr_11,r1_m_chr_12,r1_m_chr_13,r1_m_chr_1,r1_m_chr_15,
                      r1_m_chr_16,r1_m_chr_17,r1_m_chr_18,r1_m_chr_19,r1_m_chr_20,r1_m_chr_21,r1_m_chr_22)
names(rotterdam1_mom)<-c("CHROM","POS","MARKERNAME","NEA","EA","A1","TEST","N","OR","SE","Z_STAT","PVAL")
rotterdam1_mom$BETA<-log(rotterdam1_mom$OR)
rotterdam1_mom$OR_SE<-rotterdam1_mom$OR*rotterdam1_mom$SE
rotterdam1_mom$OR_95L<-rotterdam1_mom$OR-(z*rotterdam1_mom$OR_SE)
rotterdam1_mom$OR_95U<-rotterdam1_mom$OR+(z*rotterdam1_mom$OR_SE)
write.table(rotterdam1_mom,"./rotterdam1_mom_stillbirth.txt",sep="\t",row.names=FALSE, quote=FALSE)


### ROTTERDAM2 - MOTHERS - stillbirth ###

r2_m_chr_1<-as.data.frame(read.delim("./r2_m_chr_1_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_2<-as.data.frame(read.delim("./r2_m_chr_2_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_3<-as.data.frame(read.delim("./r2_m_chr_3_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_4<-as.data.frame(read.delim("./r2_m_chr_4_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_5<-as.data.frame(read.delim("./r2_m_chr_5_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_6<-as.data.frame(read.delim("./r2_m_chr_6_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_7<-as.data.frame(read.delim("./r2_m_chr_7_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_8<-as.data.frame(read.delim("./r2_m_chr_8_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_9<-as.data.frame(read.delim("./r2_m_chr_9_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_10<-as.data.frame(read.delim("./r2_m_chr_10_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_11<-as.data.frame(read.delim("./r2_m_chr_11_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_12<-as.data.frame(read.delim("./r2_m_chr_12_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_13<-as.data.frame(read.delim("./r2_m_chr_13_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_14<-as.data.frame(read.delim("./r2_m_chr_14_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_15<-as.data.frame(read.delim("./r2_m_chr_15_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_16<-as.data.frame(read.delim("./r2_m_chr_16_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_17<-as.data.frame(read.delim("./r2_m_chr_17_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_18<-as.data.frame(read.delim("./r2_m_chr_18_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_19<-as.data.frame(read.delim("./r2_m_chr_19_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_20<-as.data.frame(read.delim("./r2_m_chr_20_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_21<-as.data.frame(read.delim("./r2_m_chr_21_stillbirth.txt",header=FALSE,sep="\t"))
r2_m_chr_22<-as.data.frame(read.delim("./r2_m_chr_22_stillbirth.txt",header=FALSE,sep="\t"))

rotterdam2_mom<-rbind(r2_m_chr_1,r2_m_chr_2,r2_m_chr_3,r2_m_chr_4,r2_m_chr_5,r2_m_chr_6,r2_m_chr_7,r2_m_chr_8,
                      r2_m_chr_9,r2_m_chr_10,r2_m_chr_11,r2_m_chr_12,r2_m_chr_13,r2_m_chr_1,r2_m_chr_15,
                      r2_m_chr_16,r2_m_chr_17,r2_m_chr_18,r2_m_chr_19,r2_m_chr_20,r2_m_chr_21,r2_m_chr_22)
names(rotterdam2_mom)<-c("CHROM","POS","MARKERNAME","NEA","EA","A1","TEST","N","OR","SE","Z_STAT","PVAL")
rotterdam2_mom$BETA<-log(rotterdam2_mom$OR)
rotterdam2_mom$OR_SE<-rotterdam2_mom$OR*rotterdam2_mom$SE
rotterdam2_mom$OR_95L<-rotterdam2_mom$OR-(z*rotterdam2_mom$OR_SE)
rotterdam2_mom$OR_95U<-rotterdam2_mom$OR+(z*rotterdam2_mom$OR_SE)
write.table(rotterdam2_mom,"./rotterdam2_mom_stillbirth.txt",sep="\t",row.names=FALSE, quote=FALSE)


### NORMENT_FEB18 - MOTHERS - stillbirth ###

nf18_m_chr_1<-as.data.frame(read.delim("./nf18_m_chr_1_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_2<-as.data.frame(read.delim("./nf18_m_chr_2_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_3<-as.data.frame(read.delim("./nf18_m_chr_3_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_4<-as.data.frame(read.delim("./nf18_m_chr_4_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_5<-as.data.frame(read.delim("./nf18_m_chr_5_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_6<-as.data.frame(read.delim("./nf18_m_chr_6_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_7<-as.data.frame(read.delim("./nf18_m_chr_7_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_8<-as.data.frame(read.delim("./nf18_m_chr_8_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_9<-as.data.frame(read.delim("./nf18_m_chr_9_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_10<-as.data.frame(read.delim("./nf18_m_chr_10_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_11<-as.data.frame(read.delim("./nf18_m_chr_11_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_12<-as.data.frame(read.delim("./nf18_m_chr_12_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_13<-as.data.frame(read.delim("./nf18_m_chr_13_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_14<-as.data.frame(read.delim("./nf18_m_chr_14_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_15<-as.data.frame(read.delim("./nf18_m_chr_15_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_16<-as.data.frame(read.delim("./nf18_m_chr_16_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_17<-as.data.frame(read.delim("./nf18_m_chr_17_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_18<-as.data.frame(read.delim("./nf18_m_chr_18_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_19<-as.data.frame(read.delim("./nf18_m_chr_19_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_20<-as.data.frame(read.delim("./nf18_m_chr_20_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_21<-as.data.frame(read.delim("./nf18_m_chr_21_stillbirth.txt",header=FALSE,sep="\t"))
nf18_m_chr_22<-as.data.frame(read.delim("./nf18_m_chr_22_stillbirth.txt",header=FALSE,sep="\t"))

norment_feb18_mom<-rbind(nf18_m_chr_1,nf18_m_chr_2,nf18_m_chr_3,nf18_m_chr_4,nf18_m_chr_5,nf18_m_chr_6,nf18_m_chr_7,nf18_m_chr_8,
                         nf18_m_chr_9,nf18_m_chr_10,nf18_m_chr_11,nf18_m_chr_12,nf18_m_chr_13,nf18_m_chr_1,nf18_m_chr_15,
                         nf18_m_chr_16,nf18_m_chr_17,nf18_m_chr_18,nf18_m_chr_19,nf18_m_chr_20,nf18_m_chr_21,nf18_m_chr_22)
names(norment_feb18_mom)<-c("CHROM","POS","MARKERNAME","NEA","EA","A1","TEST","N","OR","SE","Z_STAT","PVAL")
norment_feb18_mom$BETA<-log(norment_feb18_mom$OR)
norment_feb18_mom$OR_SE<-norment_feb18_mom$OR*norment_feb18_mom$SE
norment_feb18_mom$OR_95L<-norment_feb18_mom$OR-(z*norment_feb18_mom$OR_SE)
norment_feb18_mom$OR_95U<-norment_feb18_mom$OR+(z*norment_feb18_mom$OR_SE)
write.table(norment_feb18_mom,"./norment_feb18_mom_stillbirth.txt",sep="\t",row.names=FALSE, quote=FALSE)


### NORMENT_MAY16 - MOTHERS - stillbirth ###

nm16_m_chr_1<-as.data.frame(read.delim("./nm16_m_chr_1_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_2<-as.data.frame(read.delim("./nm16_m_chr_2_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_3<-as.data.frame(read.delim("./nm16_m_chr_3_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_4<-as.data.frame(read.delim("./nm16_m_chr_4_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_5<-as.data.frame(read.delim("./nm16_m_chr_5_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_6<-as.data.frame(read.delim("./nm16_m_chr_6_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_7<-as.data.frame(read.delim("./nm16_m_chr_7_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_8<-as.data.frame(read.delim("./nm16_m_chr_8_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_9<-as.data.frame(read.delim("./nm16_m_chr_9_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_10<-as.data.frame(read.delim("./nm16_m_chr_10_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_11<-as.data.frame(read.delim("./nm16_m_chr_11_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_12<-as.data.frame(read.delim("./nm16_m_chr_12_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_13<-as.data.frame(read.delim("./nm16_m_chr_13_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_14<-as.data.frame(read.delim("./nm16_m_chr_14_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_15<-as.data.frame(read.delim("./nm16_m_chr_15_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_16<-as.data.frame(read.delim("./nm16_m_chr_16_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_17<-as.data.frame(read.delim("./nm16_m_chr_17_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_18<-as.data.frame(read.delim("./nm16_m_chr_18_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_19<-as.data.frame(read.delim("./nm16_m_chr_19_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_20<-as.data.frame(read.delim("./nm16_m_chr_20_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_21<-as.data.frame(read.delim("./nm16_m_chr_21_stillbirth.txt",header=FALSE,sep="\t"))
nm16_m_chr_22<-as.data.frame(read.delim("./nm16_m_chr_22_stillbirth.txt",header=FALSE,sep="\t"))

norment_may16_mom<-rbind(nm16_m_chr_1,nm16_m_chr_2,nm16_m_chr_3,nm16_m_chr_4,nm16_m_chr_5,nm16_m_chr_6,nm16_m_chr_7,nm16_m_chr_8,
                         nm16_m_chr_9,nm16_m_chr_10,nm16_m_chr_11,nm16_m_chr_12,nm16_m_chr_13,nm16_m_chr_1,nm16_m_chr_15,
                         nm16_m_chr_16,nm16_m_chr_17,nm16_m_chr_18,nm16_m_chr_19,nm16_m_chr_20,nm16_m_chr_21,nm16_m_chr_22)
names(norment_may16_mom)<-c("CHROM","POS","MARKERNAME","NEA","EA","A1","TEST","N","OR","SE","Z_STAT","PVAL")
norment_may16_mom$BETA<-log(norment_may16_mom$OR)
norment_may16_mom$OR_SE<-norment_may16_mom$OR*norment_may16_mom$SE
norment_may16_mom$OR_95L<-norment_may16_mom$OR-(z*norment_may16_mom$OR_SE)
norment_may16_mom$OR_95U<-norment_may16_mom$OR+(z*norment_may16_mom$OR_SE)
write.table(norment_may16_mom,"./norment_may16_mom_stillbirth.txt",sep="\t",row.names=FALSE, quote=FALSE)


### HARVEST - FATHERS - stillbirth ###

h_p_chr_1<-as.data.frame(read.delim("./h_p_chr_1_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_2<-as.data.frame(read.delim("./h_p_chr_2_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_3<-as.data.frame(read.delim("./h_p_chr_3_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_4<-as.data.frame(read.delim("./h_p_chr_4_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_5<-as.data.frame(read.delim("./h_p_chr_5_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_6<-as.data.frame(read.delim("./h_p_chr_6_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_7<-as.data.frame(read.delim("./h_p_chr_7_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_8<-as.data.frame(read.delim("./h_p_chr_8_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_9<-as.data.frame(read.delim("./h_p_chr_9_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_10<-as.data.frame(read.delim("./h_p_chr_10_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_11<-as.data.frame(read.delim("./h_p_chr_11_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_12<-as.data.frame(read.delim("./h_p_chr_12_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_13<-as.data.frame(read.delim("./h_p_chr_13_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_14<-as.data.frame(read.delim("./h_p_chr_14_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_15<-as.data.frame(read.delim("./h_p_chr_15_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_16<-as.data.frame(read.delim("./h_p_chr_16_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_17<-as.data.frame(read.delim("./h_p_chr_17_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_18<-as.data.frame(read.delim("./h_p_chr_18_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_19<-as.data.frame(read.delim("./h_p_chr_19_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_20<-as.data.frame(read.delim("./h_p_chr_20_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_21<-as.data.frame(read.delim("./h_p_chr_21_stillbirth.txt",header=FALSE,sep="\t"))
h_p_chr_22<-as.data.frame(read.delim("./h_p_chr_22_stillbirth.txt",header=FALSE,sep="\t"))

harvest_dad<-rbind(h_p_chr_1,h_p_chr_2,h_p_chr_3,h_p_chr_4,h_p_chr_5,h_p_chr_6,h_p_chr_7,h_p_chr_8,
                   h_p_chr_9,h_p_chr_10,h_p_chr_11,h_p_chr_12,h_p_chr_13,h_p_chr_1,h_p_chr_15,
                   h_p_chr_16,h_p_chr_17,h_p_chr_18,h_p_chr_19,h_p_chr_20,h_p_chr_21,h_p_chr_22)
names(harvest_dad)<-c("CHROM","POS","MARKERNAME","NEA","EA","A1","TEST","N","OR","SE","Z_STAT","PVAL")
harvest_dad$BETA<-log(harvest_dad$OR)
harvest_dad$OR_SE<-harvest_dad$OR*harvest_dad$SE
harvest_dad$OR_95L<-harvest_dad$OR-(z*harvest_dad$OR_SE)
harvest_dad$OR_95U<-harvest_dad$OR+(z*harvest_dad$OR_SE)
write.table(harvest_dad,"./harvest_dad_stillbirth.txt",sep="\t",row.names=FALSE, quote=FALSE)


### ROTTERDAM1 - FATHERS - stillbirth ###

r1_p_chr_1<-as.data.frame(read.delim("./r1_p_chr_1_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_2<-as.data.frame(read.delim("./r1_p_chr_2_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_3<-as.data.frame(read.delim("./r1_p_chr_3_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_4<-as.data.frame(read.delim("./r1_p_chr_4_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_5<-as.data.frame(read.delim("./r1_p_chr_5_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_6<-as.data.frame(read.delim("./r1_p_chr_6_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_7<-as.data.frame(read.delim("./r1_p_chr_7_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_8<-as.data.frame(read.delim("./r1_p_chr_8_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_9<-as.data.frame(read.delim("./r1_p_chr_9_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_10<-as.data.frame(read.delim("./r1_p_chr_10_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_11<-as.data.frame(read.delim("./r1_p_chr_11_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_12<-as.data.frame(read.delim("./r1_p_chr_12_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_13<-as.data.frame(read.delim("./r1_p_chr_13_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_14<-as.data.frame(read.delim("./r1_p_chr_14_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_15<-as.data.frame(read.delim("./r1_p_chr_15_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_16<-as.data.frame(read.delim("./r1_p_chr_16_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_17<-as.data.frame(read.delim("./r1_p_chr_17_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_18<-as.data.frame(read.delim("./r1_p_chr_18_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_19<-as.data.frame(read.delim("./r1_p_chr_19_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_20<-as.data.frame(read.delim("./r1_p_chr_20_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_21<-as.data.frame(read.delim("./r1_p_chr_21_stillbirth.txt",header=FALSE,sep="\t"))
r1_p_chr_22<-as.data.frame(read.delim("./r1_p_chr_22_stillbirth.txt",header=FALSE,sep="\t"))

rotterdam1_dad<-rbind(r1_p_chr_1,r1_p_chr_2,r1_p_chr_3,r1_p_chr_4,r1_p_chr_5,r1_p_chr_6,r1_p_chr_7,r1_p_chr_8,
                      r1_p_chr_9,r1_p_chr_10,r1_p_chr_11,r1_p_chr_12,r1_p_chr_13,r1_p_chr_1,r1_p_chr_15,
                      r1_p_chr_16,r1_p_chr_17,r1_p_chr_18,r1_p_chr_19,r1_p_chr_20,r1_p_chr_21,r1_p_chr_22)
names(rotterdam1_dad)<-c("CHROM","POS","MARKERNAME","NEA","EA","A1","TEST","N","OR","SE","Z_STAT","PVAL")
rotterdam1_dad$BETA<-log(rotterdam1_dad$OR)
rotterdam1_dad$OR_SE<-rotterdam1_dad$OR*rotterdam1_dad$SE
rotterdam1_dad$OR_95L<-rotterdam1_dad$OR-(z*rotterdam1_dad$OR_SE)
rotterdam1_dad$OR_95U<-rotterdam1_dad$OR+(z*rotterdam1_dad$OR_SE)
write.table(rotterdam1_dad,"./rotterdam1_dad_stillbirth.txt",sep="\t",row.names=FALSE, quote=FALSE)


### ROTTERDAM2 - FATHERS - stillbirth ###

r2_p_chr_1<-as.data.frame(read.delim("./r2_p_chr_1_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_2<-as.data.frame(read.delim("./r2_p_chr_2_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_3<-as.data.frame(read.delim("./r2_p_chr_3_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_4<-as.data.frame(read.delim("./r2_p_chr_4_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_5<-as.data.frame(read.delim("./r2_p_chr_5_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_6<-as.data.frame(read.delim("./r2_p_chr_6_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_7<-as.data.frame(read.delim("./r2_p_chr_7_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_8<-as.data.frame(read.delim("./r2_p_chr_8_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_9<-as.data.frame(read.delim("./r2_p_chr_9_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_10<-as.data.frame(read.delim("./r2_p_chr_10_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_11<-as.data.frame(read.delim("./r2_p_chr_11_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_12<-as.data.frame(read.delim("./r2_p_chr_12_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_13<-as.data.frame(read.delim("./r2_p_chr_13_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_14<-as.data.frame(read.delim("./r2_p_chr_14_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_15<-as.data.frame(read.delim("./r2_p_chr_15_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_16<-as.data.frame(read.delim("./r2_p_chr_16_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_17<-as.data.frame(read.delim("./r2_p_chr_17_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_18<-as.data.frame(read.delim("./r2_p_chr_18_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_19<-as.data.frame(read.delim("./r2_p_chr_19_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_20<-as.data.frame(read.delim("./r2_p_chr_20_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_21<-as.data.frame(read.delim("./r2_p_chr_21_stillbirth.txt",header=FALSE,sep="\t"))
r2_p_chr_22<-as.data.frame(read.delim("./r2_p_chr_22_stillbirth.txt",header=FALSE,sep="\t"))

rotterdam2_dad<-rbind(r2_p_chr_1,r2_p_chr_2,r2_p_chr_3,r2_p_chr_4,r2_p_chr_5,r2_p_chr_6,r2_p_chr_7,r2_p_chr_8,
                      r2_p_chr_9,r2_p_chr_10,r2_p_chr_11,r2_p_chr_12,r2_p_chr_13,r2_p_chr_1,r2_p_chr_15,
                      r2_p_chr_16,r2_p_chr_17,r2_p_chr_18,r2_p_chr_19,r2_p_chr_20,r2_p_chr_21,r2_p_chr_22)
names(rotterdam2_dad)<-c("CHROM","POS","MARKERNAME","NEA","EA","A1","TEST","N","OR","SE","Z_STAT","PVAL")
rotterdam2_dad$BETA<-log(rotterdam2_dad$OR)
rotterdam2_dad$OR_SE<-rotterdam2_dad$OR*rotterdam2_dad$SE
rotterdam2_dad$OR_95L<-rotterdam2_dad$OR-(z*rotterdam2_dad$OR_SE)
rotterdam2_dad$OR_95U<-rotterdam2_dad$OR+(z*rotterdam2_dad$OR_SE)
write.table(rotterdam2_dad,"./rotterdam2_dad_stillbirth.txt",sep="\t",row.names=FALSE, quote=FALSE)


### NORMENT_FEB18 - FATHERS - stillbirth ###

nf18_p_chr_1<-as.data.frame(read.delim("./nf18_p_chr_1_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_2<-as.data.frame(read.delim("./nf18_p_chr_2_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_3<-as.data.frame(read.delim("./nf18_p_chr_3_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_4<-as.data.frame(read.delim("./nf18_p_chr_4_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_5<-as.data.frame(read.delim("./nf18_p_chr_5_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_6<-as.data.frame(read.delim("./nf18_p_chr_6_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_7<-as.data.frame(read.delim("./nf18_p_chr_7_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_8<-as.data.frame(read.delim("./nf18_p_chr_8_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_9<-as.data.frame(read.delim("./nf18_p_chr_9_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_10<-as.data.frame(read.delim("./nf18_p_chr_10_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_11<-as.data.frame(read.delim("./nf18_p_chr_11_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_12<-as.data.frame(read.delim("./nf18_p_chr_12_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_13<-as.data.frame(read.delim("./nf18_p_chr_13_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_14<-as.data.frame(read.delim("./nf18_p_chr_14_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_15<-as.data.frame(read.delim("./nf18_p_chr_15_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_16<-as.data.frame(read.delim("./nf18_p_chr_16_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_17<-as.data.frame(read.delim("./nf18_p_chr_17_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_18<-as.data.frame(read.delim("./nf18_p_chr_18_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_19<-as.data.frame(read.delim("./nf18_p_chr_19_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_20<-as.data.frame(read.delim("./nf18_p_chr_20_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_21<-as.data.frame(read.delim("./nf18_p_chr_21_stillbirth.txt",header=FALSE,sep="\t"))
nf18_p_chr_22<-as.data.frame(read.delim("./nf18_p_chr_22_stillbirth.txt",header=FALSE,sep="\t"))

norment_feb18_dad<-rbind(nf18_p_chr_1,nf18_p_chr_2,nf18_p_chr_3,nf18_p_chr_4,nf18_p_chr_5,nf18_p_chr_6,nf18_p_chr_7,nf18_p_chr_8,
                         nf18_p_chr_9,nf18_p_chr_10,nf18_p_chr_11,nf18_p_chr_12,nf18_p_chr_13,nf18_p_chr_1,nf18_p_chr_15,
                         nf18_p_chr_16,nf18_p_chr_17,nf18_p_chr_18,nf18_p_chr_19,nf18_p_chr_20,nf18_p_chr_21,nf18_p_chr_22)
names(norment_feb18_dad)<-c("CHROM","POS","MARKERNAME","NEA","EA","A1","TEST","N","OR","SE","Z_STAT","PVAL")
norment_feb18_dad$BETA<-log(norment_feb18_dad$OR)
norment_feb18_dad$OR_SE<-norment_feb18_dad$OR*norment_feb18_dad$SE
norment_feb18_dad$OR_95L<-norment_feb18_dad$OR-(z*norment_feb18_dad$OR_SE)
norment_feb18_dad$OR_95U<-norment_feb18_dad$OR+(z*norment_feb18_dad$OR_SE)
write.table(norment_feb18_dad,"./norment_feb18_dad_stillbirth.txt",sep="\t",row.names=FALSE, quote=FALSE)


### NORMENT_MAY16 - FATHERS - stillbirth ###

nm16_p_chr_1<-as.data.frame(read.delim("./nm16_p_chr_1_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_2<-as.data.frame(read.delim("./nm16_p_chr_2_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_3<-as.data.frame(read.delim("./nm16_p_chr_3_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_4<-as.data.frame(read.delim("./nm16_p_chr_4_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_5<-as.data.frame(read.delim("./nm16_p_chr_5_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_6<-as.data.frame(read.delim("./nm16_p_chr_6_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_7<-as.data.frame(read.delim("./nm16_p_chr_7_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_8<-as.data.frame(read.delim("./nm16_p_chr_8_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_9<-as.data.frame(read.delim("./nm16_p_chr_9_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_10<-as.data.frame(read.delim("./nm16_p_chr_10_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_11<-as.data.frame(read.delim("./nm16_p_chr_11_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_12<-as.data.frame(read.delim("./nm16_p_chr_12_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_13<-as.data.frame(read.delim("./nm16_p_chr_13_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_14<-as.data.frame(read.delim("./nm16_p_chr_14_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_15<-as.data.frame(read.delim("./nm16_p_chr_15_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_16<-as.data.frame(read.delim("./nm16_p_chr_16_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_17<-as.data.frame(read.delim("./nm16_p_chr_17_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_18<-as.data.frame(read.delim("./nm16_p_chr_18_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_19<-as.data.frame(read.delim("./nm16_p_chr_19_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_20<-as.data.frame(read.delim("./nm16_p_chr_20_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_21<-as.data.frame(read.delim("./nm16_p_chr_21_stillbirth.txt",header=FALSE,sep="\t"))
nm16_p_chr_22<-as.data.frame(read.delim("./nm16_p_chr_22_stillbirth.txt",header=FALSE,sep="\t"))

norment_may16_dad<-rbind(nm16_p_chr_1,nm16_p_chr_2,nm16_p_chr_3,nm16_p_chr_4,nm16_p_chr_5,nm16_p_chr_6,nm16_p_chr_7,nm16_p_chr_8,
                         nm16_p_chr_9,nm16_p_chr_10,nm16_p_chr_11,nm16_p_chr_12,nm16_p_chr_13,nm16_p_chr_1,nm16_p_chr_15,
                         nm16_p_chr_16,nm16_p_chr_17,nm16_p_chr_18,nm16_p_chr_19,nm16_p_chr_20,nm16_p_chr_21,nm16_p_chr_22)
names(norment_may16_dad)<-c("CHROM","POS","MARKERNAME","NEA","EA","A1","TEST","N","OR","SE","Z_STAT","PVAL")
norment_may16_dad$BETA<-log(norment_may16_dad$OR)
norment_may16_dad$OR_SE<-norment_may16_dad$OR*norment_may16_dad$SE
norment_may16_dad$OR_95L<-norment_may16_dad$OR-(z*norment_may16_dad$OR_SE)
norment_may16_dad$OR_95U<-norment_may16_dad$OR+(z*norment_may16_dad$OR_SE)
write.table(norment_may16_dad,"./norment_may16_dad_stillbirth.txt",sep="\t",row.names=FALSE, quote=FALSE)


z<-qnorm(1-0.05/2)
vars01<-c("harvest_mom_stillbirth","rotterdam1_mom_stillbirth","rotterdam2_mom_stillbirth",
          "norment_feb18_mom_stillbirth","norment_may16_mom_stillbirth",
          "harvest_dad_stillbirth","rotterdam1_dad_stillbirth","rotterdam2_dad_stillbirth",
          "norment_feb18_dad_stillbirth","norment_may16_dad_stillbirth")

for(i in 1:length(vars01))
{
  name<-paste("./",vars01[i],".txt",sep="")
  name2<-paste("./",vars01[i],"_cat_gwama.txt",sep="")
  name3<-paste("./",vars01[i],"_cont_gwama.txt",sep="")
  dat<-as.data.frame(read.delim(name,header=TRUE,sep="\t"))
  dat_cat<-dat[,c("CHROM","POS","MARKERNAME","EA","NEA","OR","OR_95L","OR_95U","N")]
  write.table(dat_cat,name2,sep="\t",row.names=FALSE, quote=FALSE)
  dat_cont<-dat[,c("CHROM","POS","MARKERNAME","EA","NEA","BETA","SE","N")]
  write.table(dat_cont,name3,sep="\t",row.names=FALSE, quote=FALSE)
}

#############################################################################
### Copy "gwama_mom.txt" and "gwama_dad.txt" to the "gwas_results" folder ###
#############################################################################

dir.create("N:/data/durable/Projects/Hernaez_MR_CVD/gwas/preg_outc/definitive")

#################################
### RUN GWAMA - META-ANALYSIS ###
#################################






