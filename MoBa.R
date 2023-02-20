rm(list=ls())

### SELF-MADE FUNCTIONS TO OBTAIN CLEAN ESTIMATES ###
#####################################################

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

dat<-as.data.frame(fread("...",header=TRUE,sep="\t",sep2="\t"))
dat<-dat[,c(2,7:dim(dat)[2])]


# SNPs RELATED TO CAD #

cad_ma<-read.csv2("./gwas_cad_varderharst2018.csv",header=TRUE,
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

save(dat,file="./source_files/cad_grs.RData")


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


### SMALL FOR GESTATIONAL AGE (SGA, LGA) ### 

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

# SGA in boys #

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
mbrn_boys<-mbrn_boys[,c("lopenr_barn","sga_boys")]


# SGA in girls #

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

mbrn_girls<-mbrn_girls[,c("lopenr_barn","sga_girls")]

mbrn<-merge2(mbrn,mbrn_boys,by.id=c("lopenr_barn"),all.x=TRUE,sort=FALSE)
mbrn<-merge2(mbrn,mbrn_girls,by.id=c("lopenr_barn"),all.x=TRUE,sort=FALSE)

vars<-c("sga_boys","sga_girls")
for(i in 1:length(vars))
{
  mbrn[,vars[i]]<-with(mbrn,ifelse(is.na(mbrn[,vars[i]]),0,mbrn[,vars[i]]))
}

mbrn$sga<-with(mbrn,ifelse(sga_boys==1,1,
                           ifelse(sga_girls==1,1,0)))


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

mom_miscarriage_mbrn<-mbrn_mom[,c("lopenr_mor","miscarriage_mbrn")]
mom_miscarriage_mbrn<-mom_miscarriage_mbrn[order(mom_miscarriage_mbrn$lopenr_mor,-abs(mom_miscarriage_mbrn$miscarriage_mbrn)),]
mom_miscarriage_mbrn<-mom_miscarriage_mbrn[!duplicated(mom_miscarriage_mbrn$lopenr_mor),]

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
mbrn_mom<-merge2(mbrn_mom,mom_miscarriage_mbrn,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_stillbirth,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_agemax,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_yearpregmax,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_birthyear,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_paritymax,by.id=c("lopenr_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-subset2(mbrn_mom,"!is.na(mbrn_mom$mobaid_mfr_2374)")
names(mbrn_mom)<-c("lopenr_mor","mobaid_mfr_2374","flerfodsel","eclampsia_mom","hta_preg_mom","hta_chronic_mom",
                   "preterm_mom","preterm_sp_mom","gdm_mom","sga_mom","miscarriage_mbrn_mom",
                   "stillbirth_mom","agemax_mom","yearpregmax_mom","birthyear_mom","paritymax_mom")
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

dad_miscarriage_mbrn<-mbrn_dad[,c("lopenr_far","miscarriage_mbrn")]
dad_miscarriage_mbrn<-dad_miscarriage_mbrn[order(dad_miscarriage_mbrn$lopenr_far,-abs(dad_miscarriage_mbrn$miscarriage_mbrn)),]
dad_miscarriage_mbrn<-dad_miscarriage_mbrn[!duplicated(dad_miscarriage_mbrn$lopenr_far),]

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
mbrn_dad<-merge2(mbrn_dad,dad_miscarriage_mbrn,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_stillbirth,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_agemax,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_yearpregmax,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_birthyear,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_paritymax,by.id=c("lopenr_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-subset2(mbrn_dad,"!is.na(mbrn_dad$mobaid_mfr_2374)")
names(mbrn_dad)<-c("lopenr_far","mobaid_mfr_2374","flerfodsel","eclampsia_dad","hta_preg_dad",
                   "preterm_dad","preterm_sp_dad","gdm_dad","sga_dad","miscarriage_mbrn_dad",
                   "stillbirth_dad","agemax_dad","yearpregmax_dad","birthyear_dad","paritymax_dad")
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
            "miscarriage","cad_grs_mom","cad_grs_dad",
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

moba_mom<-merge2(mom_misc,mom_bmimax,by.id=c("m_id_2374"),all.x=TRUE,sort=FALSE)
moba_mom<-merge2(moba_mom,mom_eduyearsmax,by.id=c("m_id_2374"),all.x=TRUE,sort=FALSE)
moba_mom<-merge2(moba_mom,mom_smkinitmax,by.id=c("m_id_2374"),all.x=TRUE,sort=FALSE)


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

moba_dad<-merge2(dad_misc,dad_bmimax,by.id=c("f_id_2374"),all.x=TRUE,sort=FALSE)
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
             "hta_preg_mom","eclampsia_mom","gdm_mom","preterm_mom","preterm_sp_mom","sga_mom","miscarriage_mbrn_mom","stillbirth_mom")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx, method=c("birthyear_mom"=2,"agemax_mom"=2,"bmimax_mom"=2,"smkinitmax_mom"=3,"paritymax_mom"=3,
                                             "miscarriage_mbrn_mom"=3,"stillbirth_mom"=3,"eclampsia_mom"=3,
                                             "hta_preg_mom"=3,"preterm_mom"=3,"preterm_sp_mom"=3,"sga_mom"=3,"gdm_mom"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(cbind(all$descr[,1]))
colnames(tab1)<-c("Mothers-All")
write.table(tab1,file="./descriptive/descriptive_mothers.csv",sep=";",col.names=NA)


datx<-subset2(dat,"!is.na(dat$cad_grs_dad)")
datx<-datx[!duplicated(datx$f_id_2374),]

xxx<-datx[,c("birthyear_dad","agemax_dad","eduyearsmax_dad","bmimax_dad","smkinitmax_dad","paritymax_dad",
             "hta_preg_dad","eclampsia_dad","gdm_dad","preterm_dad","preterm_sp_dad","sga_dad","miscarriage_mbrn_dad","stillbirth_dad")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx, method=c("birthyear_dad"=2,"agemax_dad"=2,"bmimax_dad"=2,"smkinitmax_dad"=3,"paritymax_dad"=3,
                                             "miscarriage_mbrn_dad"=3,"stillbirth_dad"=3,"eclampsia_dad"=3,
                                             "hta_preg_dad"=3,"preterm_dad"=3,"preterm_sp_dad"=3,"sga_dad"=3,"gdm_dad"=3)),
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
             "hta_preg_mom","eclampsia_mom","gdm_mom","preterm_mom","preterm_sp_mom","sga_mom","miscarriage_mbrn_mom","stillbirth_mom")]

all<-NULL
all<-createTable(compareGroups(exclusion~.,
                               xxx, method=c("birthyear_mom"=2,"agemax_mom"=2,"bmimax_mom"=2,"smkinitmax_mom"=3,"paritymax_mom"=3,
                                             "miscarriage_mbrn_mom"=3,"stillbirth_mom"=3,"eclampsia_mom"=3,
                                             "hta_preg_mom"=3,"preterm_mom"=3,"preterm_sp_mom"=3,"sga_mom"=3,"gdm_mom"=3)),
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
             "hta_preg_dad","eclampsia_dad","gdm_dad","preterm_dad","preterm_sp_dad","sga_dad","miscarriage_mbrn_dad","stillbirth_dad")]

all<-NULL
all<-createTable(compareGroups(exclusion~.,
                               xxx, method=c("birthyear_dad"=2,"agemax_dad"=2,"bmimax_dad"=2,"smkinitmax_dad"=3,"paritymax_dad"=3,
                                             "miscarriage_mbrn_dad"=3,"stillbirth_dad"=3,"eclampsia_dad"=3,
                                             "hta_preg_dad"=3,"preterm_dad"=3,"preterm_sp_dad"=3,"sga_dad"=3,"gdm_dad"=3)),
                 show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Included","Non-included","P-value","N")
write.table(tab1,file="./descriptive/selectionbias_fathers.csv",sep=";",col.names=NA)


### LOGISTIC REGRESSION - LINEAR ASSOCIATIONS ###
#################################################

load("N:/durable/Projects/Hernaez_MR_CVD/MoBa_cvd.RData")
dat$adj<-0

vars00<-c("cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad",
          "cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad",
          "cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad",
          "cad_grs_mom","cad_grs_dad","cad_grs_mom","cad_grs_dad")
vars01<-c("cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z",
          "cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z",
          "cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z",
          "cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z")
vars02<-c("hta_preg_mom","hta_preg_dad","eclampsia_mom","eclampsia_dad","gdm_mom","gdm_dad",
          "preterm_mom","preterm_dad","preterm_sp_mom","preterm_sp_dad","sga_mom","sga_dad",
          "miscarriage_mbrn_mom","miscarriage_mbrn_dad","stillbirth_mom","stillbirth_dad")
vars03<-c("hta_chronic_mom","hta_chronic_mom","hta_chronic_mom","hta_chronic_mom",
          "adj","adj","adj","adj","adj","adj","adj","adj","adj","adj","adj","adj")
vars04<-c("cad_grs_dad","cad_grs_mom","cad_grs_dad","cad_grs_mom",
          "cad_grs_dad","cad_grs_mom","cad_grs_dad","cad_grs_mom",
          "cad_grs_dad","cad_grs_mom","cad_grs_dad","cad_grs_mom",
          "cad_grs_dad","cad_grs_mom","cad_grs_dad","cad_grs_mom")
vars05<-c("cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z",
          "cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z",
          "cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z",
          "cad_grs_dad_z","cad_grs_mom_z","cad_grs_dad_z","cad_grs_mom_z")
vars08<-c("pc01_mom","pc01_dad","pc01_mom","pc01_dad","pc01_mom","pc01_dad","pc01_mom","pc01_dad","pc01_mom","pc01_dad",
          "pc01_mom","pc01_dad","pc01_mom","pc01_dad","pc01_mom","pc01_dad")
vars09<-c("pc02_mom","pc02_dad","pc02_mom","pc02_dad","pc02_mom","pc02_dad","pc02_mom","pc02_dad","pc02_mom","pc02_dad",
          "pc02_mom","pc02_dad","pc02_mom","pc02_dad","pc02_mom","pc02_dad")
vars10<-c("pc03_mom","pc03_dad","pc03_mom","pc03_dad","pc03_mom","pc03_dad","pc03_mom","pc03_dad","pc03_mom","pc03_dad",
          "pc03_mom","pc03_dad","pc03_mom","pc03_dad","pc03_mom","pc03_dad")
vars11<-c("pc04_mom","pc04_dad","pc04_mom","pc04_dad","pc04_mom","pc04_dad","pc04_mom","pc04_dad","pc04_mom","pc04_dad",
          "pc04_mom","pc04_dad","pc04_mom","pc04_dad","pc04_mom","pc04_dad")
vars12<-c("pc05_mom","pc05_dad","pc05_mom","pc05_dad","pc05_mom","pc05_dad","pc05_mom","pc05_dad","pc05_mom","pc05_dad",
          "pc05_mom","pc05_dad","pc05_mom","pc05_dad","pc05_mom","pc05_dad")
vars13<-c("pc06_mom","pc06_dad","pc06_mom","pc06_dad","pc06_mom","pc06_dad","pc06_mom","pc06_dad","pc06_mom","pc06_dad",
          "pc06_mom","pc06_dad","pc06_mom","pc06_dad","pc06_mom","pc06_dad")
vars14<-c("pc07_mom","pc07_dad","pc07_mom","pc07_dad","pc07_mom","pc07_dad","pc07_mom","pc07_dad","pc07_mom","pc07_dad",
          "pc07_mom","pc07_dad","pc07_mom","pc07_dad","pc07_mom","pc07_dad")
vars15<-c("pc08_mom","pc08_dad","pc08_mom","pc08_dad","pc08_mom","pc08_dad","pc08_mom","pc08_dad","pc08_mom","pc08_dad",
          "pc08_mom","pc08_dad","pc08_mom","pc08_dad","pc08_mom","pc08_dad")
vars16<-c("pc09_mom","pc09_dad","pc09_mom","pc09_dad","pc09_mom","pc09_dad","pc09_mom","pc09_dad","pc09_mom","pc09_dad",
          "pc09_mom","pc09_dad","pc09_mom","pc09_dad","pc09_mom","pc09_dad")
vars17<-c("pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad",
          "pc10_mom","pc10_dad","pc10_mom","pc10_dad","pc10_mom","pc10_dad")
vars18<-c("genotype_batch_mom","genotype_batch_dad","genotype_batch_mom","genotype_batch_dad","genotype_batch_mom","genotype_batch_dad",
          "genotype_batch_mom","genotype_batch_dad","genotype_batch_mom","genotype_batch_dad","genotype_batch_mom","genotype_batch_dad",
          "genotype_batch_mom","genotype_batch_dad","genotype_batch_mom","genotype_batch_dad")
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

