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
library(Hmisc)
library(survival)
library(foreign)
library(bitops)
library(caTools)
library(gplots)
library(ROCR)
library(mice)
library(writexl)
library(officer)
library(uuid)
library(HardyWeinberg)
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
library(fmsb)
library(remotes)
library(sessioninfo)
library(usethis)
library(iterpc)
library(CompQuadForm)
library(arrangements)
library(rjson)

library(devtools)
library(googleAuthR)
library(MendelianRandomization)
library(mr.raps)
library(meta)

#devtools::install("N:/durable/projects/ALHE_MR_CHD/R_packages/MR-PRESSO-master/MR-PRESSO-master/")
#devtools::install("N:/durable/projects/ALHE_MR_CHD/R_packages/MRInstruments-master/MRInstruments-master/")
#devtools::install("N:/durable/projects/ALHE_MR_CHD/R_packages/MRMix-master/MRMix-master/")
#utils::install.packages("N:/durable/projects/ALHE_MR_CHD/R_packages/RadialMR-master/RadialMR-master/", repos=NULL, type="source")
#utils::install.packages("N:/durable/projects/ALHE_MR_CHD/R_packages/ieugwasr-0.1.5.tar.gz", repos=NULL, type="source")
#utils::install.packages("N:/durable/projects/ALHE_MR_CHD/R_packages/TwoSampleMR-0.5.6.tar.gz", repos=NULL, type="source")

library(MRPRESSO)
library(MRInstruments)
library(MRMix)
library(RadialMR)
library(ieugwasr)
library(TwoSampleMR)


RutinesLocals<- "N:/durable/projects/ALHE_MR_CHD/R_packages/Routines"
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
  redondeo<-ifelse(abs(x)<0.0001,signif(x,1),
                   ifelse(abs(x)<0.001,signif(x,1),
                          ifelse(abs(x)<0.1,round(x,3),
                                 ifelse(abs(x)<1,round(x,2),signif(x,3)))))
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


dir.create("N:/durable/projects/ALHE_MR_CHD/Outputs")
dir.create("N:/durable/projects/ALHE_MR_CHD/Outputs/descriptive")
dir.create("N:/durable/projects/ALHE_MR_CHD/Outputs/results")

setwd("N:/durable/projects/ALHE_MR_CHD")


### GENERATION OF DATABASE ###
##############################

dat<-as.data.frame(read.delim("./Data/pop_hunt_220817.txt",header=TRUE,sep=" "))
names(dat)<-tolower(names(dat))

dat$pid<-as.character(dat$pid110334)
dat$sex<-with(dat,ifelse(sex=="M",0,
                         ifelse(sex=="F",1,NA)))
dat$fert<-with(dat,ifelse(fert=="fertile",0,
                          ifelse(fert=="subfertile",1,
                                 ifelse(fert=="infertile",2,NA))))
dat$subf<-with(dat,ifelse(fert==0,0,
                          ifelse(fert==1,1,
                                 ifelse(fert==2,1,NA))))
dat$smok<-with(dat,ifelse(smok=="non-smoker",0,
                          ifelse(smok=="smoker",1,NA)))
dat$educ<-with(dat,ifelse(educ=="secondary school",1,
                          ifelse(educ=="upper secondary school",2,
                                 ifelse(educ=="higher education",3,NA))))
dat$eduyears<-with(dat,ifelse(educ==1,10,
                              ifelse(educ==2,13,
                                     ifelse(educ==3,19,NA))))
dat$children_questionnaire<-with(dat,ifelse(children_questionnaire=="8<",8,children_questionnaire))
dat$children_questionnaire<-as.numeric(dat$children_questionnaire)
dat$children_mbrn<-with(dat,ifelse(children_mbrn=="8<",8,children_mbrn))
dat$children_mbrn<-as.numeric(dat$children_mbrn)

#dat$children<-with(dat,ifelse(is.na(children_questionnaire) & is.na(children_mbrn),0,
#                              ifelse(is.na(children_questionnaire) & !is.na(children_mbrn),children_mbrn,
#                                     ifelse(!is.na(children_questionnaire) & is.na(children_mbrn),children_questionnaire,
#                                            ifelse(!is.na(children_questionnaire) & !is.na(children_mbrn),max(dat$children_questionnaire,dat$children_mbrn,na.rm=T),NA)))))
dat$children <- apply(cbind(dat$children_questionnaire,dat$children_mbrn),1,max, na.rm=TRUE)
dat$children<-with(dat,ifelse(children=="-Inf",0,children))

#dat$children<-with(dat,ifelse(children>4,4,children))

dat$partner<-with(dat,ifelse(is.na(partner),NA,
                             ifelse(!is.na(partner),as.character(substr(partner,4,16)),NA)))
dat$chd_prior<-with(dat,ifelse(chd_prior=="FALSE",0,
                               ifelse(chd_prior=="TRUE",1,NA))) 
dat$chd_lifetime<-with(dat,ifelse(chd_lifetime=="FALSE",0,
                                  ifelse(chd_lifetime=="TRUE",1,NA)))  
dat<-rename.vars(dat,from=c("survey_age","smok","educ"),to=c("age_survey","smoking","education"))

dat<-dat[,c("pid","partner","couple","survey","age_survey","birthyear",
            "sex","bmi","eduyears","smoking","children","subf","fert",
            "chd_lifetime","chd_prior","chd_noangina_lifetime","chd_noangina_prior","chd_noangina_icd10_lifetime",
            "chd_noangina_icd10_prior","infarct_lifetime","infarct_prior","angina_lifetime","angina_prior")]


grs<-as.data.frame(read.delim("./Data/grs_chd_220803.txt",header=TRUE,sep=" "))
names(grs)<-tolower(names(grs))
grs$pid<-as.character(substr(grs$pid,4,16))
grs<-rename.vars(grs,from=c("score"),to=c("cad_grs"))
grs$score_std<-NULL
grs$sex<-NULL
grs$birthyear<-NULL
dat<-merge2(dat,grs,by.id=c("pid"),all.x=TRUE,sort=FALSE)
save(dat,file="N:/durable/projects/ALHE_MR_CHD/all.RData")


### UNIQUE VARIABLES IN HUNT ###
################################

# Unique variables in MoBa #

dat2<-subset2(dat,"dat$sex==1")

mom_subf<-dat2[,c("pid","subf")]
mom_subf<-mom_subf[order(mom_subf$pid,-abs(mom_subf$subf)),]
mom_subf<-mom_subf[!duplicated(mom_subf$pid),]

mom_bmimax<-dat2[,c("pid","bmi")]
mom_bmimax<-mom_bmimax[order(mom_bmimax$pid,-abs(mom_bmimax$bmi)),]
mom_bmimax<-mom_bmimax[!duplicated(mom_bmimax$pid),]
mom_bmimax<-rename.vars(mom_bmimax,from=c("bmi"),to=c("bmimax"))

mom_eduyearsmax<-dat2[,c("pid","eduyears")]
mom_eduyearsmax<-mom_eduyearsmax[order(mom_eduyearsmax$pid,-abs(mom_eduyearsmax$eduyears)),]
mom_eduyearsmax<-mom_eduyearsmax[!duplicated(mom_eduyearsmax$pid),]
mom_eduyearsmax<-rename.vars(mom_eduyearsmax,from=c("eduyears"),to=c("eduyearsmax"))

mom_smokingmax<-dat2[,c("pid","smoking")]
mom_smokingmax<-mom_smokingmax[order(mom_smokingmax$pid,-abs(mom_smokingmax$smoking)),]
mom_smokingmax<-mom_smokingmax[!duplicated(mom_smokingmax$pid),]
mom_smokingmax<-rename.vars(mom_smokingmax,from=c("smoking"),to=c("smokingmax"))

mom_childrenmax<-dat2[,c("pid","children")]
mom_childrenmax<-mom_childrenmax[order(mom_childrenmax$pid,-abs(mom_childrenmax$children)),]
mom_childrenmax<-mom_childrenmax[!duplicated(mom_childrenmax$pid),]
mom_childrenmax<-rename.vars(mom_childrenmax,from=c("children"),to=c("childrenmax"))


dat2<-dat2[,c("pid","sex","partner","couple","survey","children","age_survey","birthyear",
              "chd_lifetime","chd_prior","chd_noangina_lifetime","chd_noangina_prior","chd_noangina_icd10_lifetime",
              "chd_noangina_icd10_prior","infarct_lifetime","infarct_prior","angina_lifetime","angina_prior",
              "cad_grs","batch","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",
              "pc11","pc12","pc13","pc14","pc15","pc16","pc17","pc18","pc19","pc20")]

hunt_mom<-merge2(mom_subf,mom_bmimax,by.id=c("pid"),all.x=TRUE,sort=FALSE)
hunt_mom<-merge2(hunt_mom,mom_eduyearsmax,by.id=c("pid"),all.x=TRUE,sort=FALSE)
hunt_mom<-merge2(hunt_mom,mom_smokingmax,by.id=c("pid"),all.x=TRUE,sort=FALSE)
hunt_mom<-merge2(hunt_mom,mom_childrenmax,by.id=c("pid"),all.x=TRUE,sort=FALSE)
hunt_mom<-merge2(hunt_mom,dat2,by.id=c("pid"),all.x=TRUE,sort=FALSE)


dat2<-subset2(dat,"dat$sex==0")

dad_subf<-dat2[,c("pid","subf")]
dad_subf<-dad_subf[order(dad_subf$pid,-abs(dad_subf$subf)),]
dad_subf<-dad_subf[!duplicated(dad_subf$pid),]

dad_bmimax<-dat2[,c("pid","bmi")]
dad_bmimax<-dad_bmimax[order(dad_bmimax$pid,-abs(dad_bmimax$bmi)),]
dad_bmimax<-dad_bmimax[!duplicated(dad_bmimax$pid),]
dad_bmimax<-rename.vars(dad_bmimax,from=c("bmi"),to=c("bmimax"))

dad_eduyearsmax<-dat2[,c("pid","eduyears")]
dad_eduyearsmax<-dad_eduyearsmax[order(dad_eduyearsmax$pid,-abs(dad_eduyearsmax$eduyears)),]
dad_eduyearsmax<-dad_eduyearsmax[!duplicated(dad_eduyearsmax$pid),]
dad_eduyearsmax<-rename.vars(dad_eduyearsmax,from=c("eduyears"),to=c("eduyearsmax"))

dad_smokingmax<-dat2[,c("pid","smoking")]
dad_smokingmax<-dad_smokingmax[order(dad_smokingmax$pid,-abs(dad_smokingmax$smoking)),]
dad_smokingmax<-dad_smokingmax[!duplicated(dad_smokingmax$pid),]
dad_smokingmax<-rename.vars(dad_smokingmax,from=c("smoking"),to=c("smokingmax"))

dad_childrenmax<-dat2[,c("pid","children")]
dad_childrenmax<-dad_childrenmax[order(dad_childrenmax$pid,-abs(dad_childrenmax$children)),]
dad_childrenmax<-dad_childrenmax[!duplicated(dad_childrenmax$pid),]
dad_childrenmax<-rename.vars(dad_childrenmax,from=c("children"),to=c("childrenmax"))


dat2<-dat2[,c("pid","sex","partner","couple","survey","children","age_survey","birthyear",
              "chd_lifetime","chd_prior","chd_noangina_lifetime","chd_noangina_prior","chd_noangina_icd10_lifetime",
              "chd_noangina_icd10_prior","infarct_lifetime","infarct_prior","angina_lifetime","angina_prior",
              "cad_grs","batch","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",
              "pc11","pc12","pc13","pc14","pc15","pc16","pc17","pc18","pc19","pc20")]

hunt_dad<-merge2(dad_subf,dad_bmimax,by.id=c("pid"),all.x=TRUE,sort=FALSE)
hunt_dad<-merge2(hunt_dad,dad_eduyearsmax,by.id=c("pid"),all.x=TRUE,sort=FALSE)
hunt_dad<-merge2(hunt_dad,dad_smokingmax,by.id=c("pid"),all.x=TRUE,sort=FALSE)
hunt_dad<-merge2(hunt_dad,dad_childrenmax,by.id=c("pid"),all.x=TRUE,sort=FALSE)
hunt_dad<-merge2(hunt_dad,dat2,by.id=c("pid"),all.x=TRUE,sort=FALSE)


### BIRTH REGISTRY ###
######################

mbrn<-spss.get("N:/durable/RAW/MFR/MFRdata_202250_2020Q4_versjon2.sav",
               use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(mbrn)<-tolower(names(mbrn))
mbrn$pid110334_barn<-gsub(" ","",mbrn$pid110334_barn)
mbrn$pid110334_barn2<-with(mbrn,ifelse(pid110334_barn=="",NA,pid110334_barn))


### PREECLAMPSIA / ECLAMPSIA ###

vars<-c("preekl","preekltidl","eklampsi","hellp","hypertensjon_alene","hypertensjon_kronisk","svlen",
        "diabetes_mellitus","spabort_12","spabort_12_5","spabort_23","spabort_23_5","dodkat","faar","fstart")
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
                                  ifelse(preterm==1 & fstart>1,NA,
                                         ifelse(preterm==0,0,NA))))

### GESTATIONAL DIABETES ###                                               

mbrn$gdm<-with(mbrn,ifelse(diabetes_mellitus==4,1,0))


### STILLBIRTH ###                                               

mbrn$stillbirth<-with(mbrn,ifelse(dodkat==7,1,
                                  ifelse(dodkat==8,1,
                                         ifelse(dodkat==9,1,0))))


### MISCARRIAGE_MBRN ###                                               

mbrn$miscarriage12_mbrn<-with(mbrn,ifelse(spabort_12>0 | spabort_12_5>0,1,0))
mbrn$miscarriage24_mbrn<-with(mbrn,ifelse(spabort_23>0 | spabort_23_5>0,1,0))
mbrn$miscarriage_mbrn<-with(mbrn,ifelse(miscarriage12_mbrn>0 | miscarriage24_mbrn>0,1,0))


### SMALL AND LARGE FOR GESTATIONAL AGE (SGA, LGA) ### 

# Erase children with sex not specified (0), uncertain (3) or missing (9) #

mbrn$kjonn<-with(mbrn,ifelse(kjonn==1,1,
                             ifelse(kjonn==2,2,NA)))
mbrn<-subset2(mbrn,"!is.na(mbrn$kjonn)")

mbrn$svlen2<-with(mbrn,ifelse(svlen<25,24,
                              ifelse(svlen>42,43,svlen)))
mbrn$vekt<-with(mbrn,ifelse(vekt<500,NA,
                            ifelse(vekt>6500,NA,vekt)))

mbrn_boys<-subset2(mbrn,"mbrn$kjonn==1 & !is.na(mbrn$pid110334_barn2)")
mbrn_boys<-mbrn_boys[,c("pid110334_barn","svlen2","vekt")]
mbrn_girls<-subset2(mbrn,"mbrn$kjonn==2 & !is.na(mbrn$pid110334_barn2)")
mbrn_girls<-mbrn_girls[,c("pid110334_barn","svlen2","vekt")]

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
mbrn_boys<-mbrn_boys[,c("pid110334_barn","sga_boys","lga_boys")]


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

mbrn_girls<-mbrn_girls[,c("pid110334_barn","sga_girls","lga_girls")]

mbrn<-merge2(mbrn,mbrn_boys,by.id=c("pid110334_barn"),all.x=TRUE,sort=FALSE)
mbrn<-merge2(mbrn,mbrn_girls,by.id=c("pid110334_barn"),all.x=TRUE,sort=FALSE)

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

mbrn$flerfodsel<-with(mbrn,ifelse(is.na(flerfodsel),0,flerfodsel))

# Mothers #

mbrn$pid110334_mor<-gsub(" ","",mbrn$pid110334_mor)
mbrn$pid110334_mor2<-with(mbrn,ifelse(pid110334_mor=="",NA,pid110334_mor))
mbrn_mom<-subset2(mbrn,"!is.na(mbrn$pid110334_mor2)")

mom_agemax<-mbrn_mom[,c("pid110334_mor","mors_alder")]
mom_agemax<-mom_agemax[order(mom_agemax$pid110334_mor,-abs(mom_agemax$mors_alder)),]
mom_agemax<-mom_agemax[!duplicated(mom_agemax$pid110334_mor),]

mom_yearpregmax<-mbrn_mom[,c("pid110334_mor","faar")]
mom_yearpregmax<-mom_yearpregmax[order(mom_yearpregmax$pid110334_mor,-abs(mom_yearpregmax$faar)),]
mom_yearpregmax<-mom_yearpregmax[!duplicated(mom_yearpregmax$pid110334_mor),]

mom_paritymax<-mbrn_mom[,c("pid110334_mor","paritet_5")]
mom_paritymax<-mom_paritymax[order(mom_paritymax$pid110334_mor,-abs(mom_paritymax$paritet_5)),]
mom_paritymax<-mom_paritymax[!duplicated(mom_paritymax$pid110334_mor),]

mom_eclampsia<-mbrn_mom[,c("pid110334_mor","eclampsia")]
mom_eclampsia<-mom_eclampsia[order(mom_eclampsia$pid110334_mor,-abs(mom_eclampsia$eclampsia)),]
mom_eclampsia<-mom_eclampsia[!duplicated(mom_eclampsia$pid110334_mor),]

mom_hta_preg<-mbrn_mom[,c("pid110334_mor","hta_preg")]
mom_hta_preg<-mom_hta_preg[order(mom_hta_preg$pid110334_mor,-abs(mom_hta_preg$hta_preg)),]
mom_hta_preg<-mom_hta_preg[!duplicated(mom_hta_preg$pid110334_mor),]

mom_hta_chronic<-mbrn_mom[,c("pid110334_mor","hta_chronic")]
mom_hta_chronic<-mom_hta_chronic[order(mom_hta_chronic$pid110334_mor,-abs(mom_hta_chronic$hta_chronic)),]
mom_hta_chronic<-mom_hta_chronic[!duplicated(mom_hta_chronic$pid110334_mor),]

mom_preterm<-mbrn_mom[,c("pid110334_mor","preterm")]
mom_preterm<-mom_preterm[order(mom_preterm$pid110334_mor,-abs(mom_preterm$preterm)),]
mom_preterm<-mom_preterm[!duplicated(mom_preterm$pid110334_mor),]

mom_preterm_sp<-mbrn_mom[,c("pid110334_mor","preterm_sp")]
mom_preterm_sp<-mom_preterm_sp[order(mom_preterm_sp$pid110334_mor,-abs(mom_preterm_sp$preterm_sp)),]
mom_preterm_sp<-mom_preterm_sp[!duplicated(mom_preterm_sp$pid110334_mor),]

mom_gdm<-mbrn_mom[,c("pid110334_mor","gdm")]
mom_gdm<-mom_gdm[order(mom_gdm$pid110334_mor,-abs(mom_gdm$gdm)),]
mom_gdm<-mom_gdm[!duplicated(mom_gdm$pid110334_mor),]

mom_sga<-mbrn_mom[,c("pid110334_mor","sga")]
mom_sga<-mom_sga[order(mom_sga$pid110334_mor,-abs(mom_sga$sga)),]
mom_sga<-mom_sga[!duplicated(mom_sga$pid110334_mor),]

mom_lga<-mbrn_mom[,c("pid110334_mor","lga")]
mom_lga<-mom_lga[order(mom_lga$pid110334_mor,-abs(mom_lga$lga)),]
mom_lga<-mom_lga[!duplicated(mom_lga$pid110334_mor),]

mom_miscarriage12_mbrn<-mbrn_mom[,c("pid110334_mor","miscarriage12_mbrn")]
mom_miscarriage12_mbrn<-mom_miscarriage12_mbrn[order(mom_miscarriage12_mbrn$pid110334_mor,-abs(mom_miscarriage12_mbrn$miscarriage12_mbrn)),]
mom_miscarriage12_mbrn<-mom_miscarriage12_mbrn[!duplicated(mom_miscarriage12_mbrn$pid110334_mor),]

mom_miscarriage24_mbrn<-mbrn_mom[,c("pid110334_mor","miscarriage24_mbrn")]
mom_miscarriage24_mbrn<-mom_miscarriage24_mbrn[order(mom_miscarriage24_mbrn$pid110334_mor,-abs(mom_miscarriage24_mbrn$miscarriage24_mbrn)),]
mom_miscarriage24_mbrn<-mom_miscarriage24_mbrn[!duplicated(mom_miscarriage24_mbrn$pid110334_mor),]

mom_miscarriage_mbrn<-mbrn_mom[,c("pid110334_mor","miscarriage_mbrn")]
mom_miscarriage_mbrn<-mom_miscarriage_mbrn[order(mom_miscarriage_mbrn$pid110334_mor,-abs(mom_miscarriage_mbrn$miscarriage_mbrn)),]
mom_miscarriage_mbrn<-mom_miscarriage_mbrn[!duplicated(mom_miscarriage_mbrn$pid110334_mor),]

mom_stillbirth<-mbrn_mom[,c("pid110334_mor","stillbirth")]
mom_stillbirth<-mom_stillbirth[order(mom_stillbirth$pid110334_mor,-abs(mom_stillbirth$stillbirth)),]
mom_stillbirth<-mom_stillbirth[!duplicated(mom_stillbirth$pid110334_mor),]


mbrn_mom<-subset2(mbrn,"!is.na(mbrn$pid110334_mor2) & mbrn$flerfodsel==0")
mbrn_mom<-as.data.frame(mbrn_mom[,c("pid110334_mor")])
names(mbrn_mom)<-c("pid110334_mor")
mbrn_mom<-as.data.frame(mbrn_mom[!duplicated(mbrn_mom$pid110334_mor),])
names(mbrn_mom)<-c("pid110334_mor")

mbrn_mom<-merge2(mbrn_mom,mom_eclampsia,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_hta_preg,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_hta_chronic,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_preterm,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_preterm_sp,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_gdm,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_sga,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_lga,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_miscarriage12_mbrn,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_miscarriage24_mbrn,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_miscarriage_mbrn,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_stillbirth,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_agemax,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_yearpregmax,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_paritymax,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
names(mbrn_mom)<-c("pid","eclampsia","hta_preg","hta_chronic","preterm","preterm_sp","gdm","sga","lga",
                   "miscarriage12_mbrn","miscarriage24_mbrn","miscarriage_mbrn","stillbirth",
                   "agemax","yearpregmax","paritymax")

mbrn_mom$pid<-with(mbrn_mom,ifelse(pid=="",NA,pid))
mbrn_mom<-subset2(mbrn_mom,"!is.na(mbrn_mom$pid)")


# Fathers #

mbrn$pid110334_far<-gsub(" ","",mbrn$pid110334_far)
mbrn$pid110334_far2<-with(mbrn,ifelse(pid110334_far=="",NA,pid110334_far))
mbrn_dad<-subset2(mbrn,"!is.na(mbrn$pid110334_far2)")

dad_agemax<-mbrn_dad[,c("pid110334_far","fars_alder")]
dad_agemax<-dad_agemax[order(dad_agemax$pid110334_far,-abs(dad_agemax$fars_alder)),]
dad_agemax<-dad_agemax[!duplicated(dad_agemax$pid110334_far),]

dad_yearpregmax<-mbrn_dad[,c("pid110334_far","faar")]
dad_yearpregmax<-dad_yearpregmax[order(dad_yearpregmax$pid110334_far,-abs(dad_yearpregmax$faar)),]
dad_yearpregmax<-dad_yearpregmax[!duplicated(dad_yearpregmax$pid110334_far),]

dad_paritymax<-mbrn_dad[,c("pid110334_far","paritet_5")]
dad_paritymax<-dad_paritymax[order(dad_paritymax$pid110334_far,-abs(dad_paritymax$paritet_5)),]
dad_paritymax<-dad_paritymax[!duplicated(dad_paritymax$pid110334_far),]

dad_eclampsia<-mbrn_dad[,c("pid110334_far","eclampsia")]
dad_eclampsia<-dad_eclampsia[order(dad_eclampsia$pid110334_far,-abs(dad_eclampsia$eclampsia)),]
dad_eclampsia<-dad_eclampsia[!duplicated(dad_eclampsia$pid110334_far),]

dad_hta_preg<-mbrn_dad[,c("pid110334_far","hta_preg")]
dad_hta_preg<-dad_hta_preg[order(dad_hta_preg$pid110334_far,-abs(dad_hta_preg$hta_preg)),]
dad_hta_preg<-dad_hta_preg[!duplicated(dad_hta_preg$pid110334_far),]

dad_hta_chronic<-mbrn_dad[,c("pid110334_far","hta_chronic")]
dad_hta_chronic<-dad_hta_chronic[order(dad_hta_chronic$pid110334_far,-abs(dad_hta_chronic$hta_chronic)),]
dad_hta_chronic<-dad_hta_chronic[!duplicated(dad_hta_chronic$pid110334_far),]

dad_preterm<-mbrn_dad[,c("pid110334_far","preterm")]
dad_preterm<-dad_preterm[order(dad_preterm$pid110334_far,-abs(dad_preterm$preterm)),]
dad_preterm<-dad_preterm[!duplicated(dad_preterm$pid110334_far),]

dad_preterm_sp<-mbrn_dad[,c("pid110334_far","preterm_sp")]
dad_preterm_sp<-dad_preterm_sp[order(dad_preterm_sp$pid110334_far,-abs(dad_preterm_sp$preterm_sp)),]
dad_preterm_sp<-dad_preterm_sp[!duplicated(dad_preterm_sp$pid110334_far),]

dad_gdm<-mbrn_dad[,c("pid110334_far","gdm")]
dad_gdm<-dad_gdm[order(dad_gdm$pid110334_far,-abs(dad_gdm$gdm)),]
dad_gdm<-dad_gdm[!duplicated(dad_gdm$pid110334_far),]

dad_sga<-mbrn_dad[,c("pid110334_far","sga")]
dad_sga<-dad_sga[order(dad_sga$pid110334_far,-abs(dad_sga$sga)),]
dad_sga<-dad_sga[!duplicated(dad_sga$pid110334_far),]

dad_lga<-mbrn_dad[,c("pid110334_far","lga")]
dad_lga<-dad_lga[order(dad_lga$pid110334_far,-abs(dad_lga$lga)),]
dad_lga<-dad_lga[!duplicated(dad_lga$pid110334_far),]

dad_miscarriage12_mbrn<-mbrn_dad[,c("pid110334_far","miscarriage12_mbrn")]
dad_miscarriage12_mbrn<-dad_miscarriage12_mbrn[order(dad_miscarriage12_mbrn$pid110334_far,-abs(dad_miscarriage12_mbrn$miscarriage12_mbrn)),]
dad_miscarriage12_mbrn<-dad_miscarriage12_mbrn[!duplicated(dad_miscarriage12_mbrn$pid110334_far),]

dad_miscarriage24_mbrn<-mbrn_dad[,c("pid110334_far","miscarriage24_mbrn")]
dad_miscarriage24_mbrn<-dad_miscarriage24_mbrn[order(dad_miscarriage24_mbrn$pid110334_far,-abs(dad_miscarriage24_mbrn$miscarriage24_mbrn)),]
dad_miscarriage24_mbrn<-dad_miscarriage24_mbrn[!duplicated(dad_miscarriage24_mbrn$pid110334_far),]

dad_miscarriage_mbrn<-mbrn_dad[,c("pid110334_far","miscarriage_mbrn")]
dad_miscarriage_mbrn<-dad_miscarriage_mbrn[order(dad_miscarriage_mbrn$pid110334_far,-abs(dad_miscarriage_mbrn$miscarriage_mbrn)),]
dad_miscarriage_mbrn<-dad_miscarriage_mbrn[!duplicated(dad_miscarriage_mbrn$pid110334_far),]

dad_stillbirth<-mbrn_dad[,c("pid110334_far","stillbirth")]
dad_stillbirth<-dad_stillbirth[order(dad_stillbirth$pid110334_far,-abs(dad_stillbirth$stillbirth)),]
dad_stillbirth<-dad_stillbirth[!duplicated(dad_stillbirth$pid110334_far),]


mbrn_dad<-subset2(mbrn,"!is.na(mbrn$pid110334_far2) & mbrn$flerfodsel==0")
mbrn_dad<-as.data.frame(mbrn_dad[,c("pid110334_far")])
names(mbrn_dad)<-c("pid110334_far")
mbrn_dad<-as.data.frame(mbrn_dad[!duplicated(mbrn_dad$pid110334_far),])
names(mbrn_dad)<-c("pid110334_far")

mbrn_dad<-merge2(mbrn_dad,dad_eclampsia,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_hta_preg,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_hta_chronic,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_preterm,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_preterm_sp,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_gdm,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_sga,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_lga,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_miscarriage12_mbrn,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_miscarriage24_mbrn,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_miscarriage_mbrn,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_stillbirth,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_agemax,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_yearpregmax,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_paritymax,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
names(mbrn_dad)<-c("pid","eclampsia","hta_preg","hta_chronic","preterm","preterm_sp","gdm","sga","lga",
                   "miscarriage12_mbrn","miscarriage24_mbrn","miscarriage_mbrn","stillbirth",
                   "agemax","yearpregmax","paritymax")

mbrn_dad$pid<-with(mbrn_dad,ifelse(pid=="",NA,pid))
mbrn_dad<-subset2(mbrn_dad,"!is.na(mbrn_dad$pid)")


mom<-merge2(hunt_mom,mbrn_mom,by.id=c("pid"),all.x=TRUE,sort=FALSE)
mom<-mom[order(mom$pid,-abs(mom$couple)),]
mom<-mom[!duplicated(mom$pid),]
mom$agemax<-mom$yearpregmax-mom$birthyear
save(mom,file="N:/durable/projects/ALHE_MR_CHD/mom.RData")

dad<-merge2(hunt_dad,mbrn_dad,by.id=c("pid"),all.x=TRUE,sort=FALSE)
dad<-dad[order(dad$pid,-abs(dad$couple)),]
dad<-dad[!duplicated(dad$pid),]
dad$agemax<-dad$yearpregmax-dad$birthyear
save(dad,file="N:/durable/projects/ALHE_MR_CHD/dad.RData")

dat<-rbind(mom,dad)
save(dat,file="N:/durable/projects/ALHE_MR_CHD/all.RData")


