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
                                                             ifelse(abs(x)>=100,round(x,0),round(x,0)))))))))
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
dat$chd_lifetime<-with(dat,ifelse(chd_lifetime=="FALSE",0,
                                  ifelse(chd_lifetime=="TRUE",1,NA)))  
dat<-rename.vars(dat,from=c("survey_age","smok","educ"),to=c("age_survey","smoking","education"))

dat<-dat[,c("pid","partner","couple","survey","age_survey","birthyear",
            "sex","bmi","eduyears","smoking","children",
            "chd_lifetime")]


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
              "chd_lifetime",
              "cad_grs","batch","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",
              "pc11","pc12","pc13","pc14","pc15","pc16","pc17","pc18","pc19","pc20")]

hunt_mom<-merge2(mom_bmimax,mom_eduyearsmax,by.id=c("pid"),all.x=TRUE,sort=FALSE)
hunt_mom<-merge2(hunt_mom,mom_smokingmax,by.id=c("pid"),all.x=TRUE,sort=FALSE)
hunt_mom<-merge2(hunt_mom,mom_childrenmax,by.id=c("pid"),all.x=TRUE,sort=FALSE)
hunt_mom<-merge2(hunt_mom,dat2,by.id=c("pid"),all.x=TRUE,sort=FALSE)


dat2<-subset2(dat,"dat$sex==0")

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
              "chd_lifetime",
              "cad_grs","batch","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",
              "pc11","pc12","pc13","pc14","pc15","pc16","pc17","pc18","pc19","pc20")]

hunt_dad<-merge2(dad_bmimax,dad_eduyearsmax,by.id=c("pid"),all.x=TRUE,sort=FALSE)
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


### SMALL FOR GESTATIONAL AGE (SGA, LGA) ### 

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
mbrn_boys<-mbrn_boys[,c("pid110334_barn","sga_boys")]


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

mbrn_girls<-mbrn_girls[,c("pid110334_barn","sga_girls")]

mbrn<-merge2(mbrn,mbrn_boys,by.id=c("pid110334_barn"),all.x=TRUE,sort=FALSE)
mbrn<-merge2(mbrn,mbrn_girls,by.id=c("pid110334_barn"),all.x=TRUE,sort=FALSE)

vars<-c("sga_boys","sga_girls")
for(i in 1:length(vars))
{
  mbrn[,vars[i]]<-with(mbrn,ifelse(is.na(mbrn[,vars[i]]),0,mbrn[,vars[i]]))
}

mbrn$sga<-with(mbrn,ifelse(sga_boys==1,1,
                           ifelse(sga_girls==1,1,0)))


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
mbrn_mom<-merge2(mbrn_mom,mom_miscarriage_mbrn,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_stillbirth,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_agemax,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_yearpregmax,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
mbrn_mom<-merge2(mbrn_mom,mom_paritymax,by.id=c("pid110334_mor"),all.x=TRUE,sort=FALSE)
names(mbrn_mom)<-c("pid","eclampsia","hta_preg","hta_chronic","preterm","preterm_sp","gdm","sga",
                   "miscarriage_mbrn","stillbirth",
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
mbrn_dad<-merge2(mbrn_dad,dad_miscarriage_mbrn,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_stillbirth,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_agemax,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_yearpregmax,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
mbrn_dad<-merge2(mbrn_dad,dad_paritymax,by.id=c("pid110334_far"),all.x=TRUE,sort=FALSE)
names(mbrn_dad)<-c("pid","eclampsia","hta_preg","hta_chronic","preterm","preterm_sp","gdm","sga",
                   "miscarriage_mbrn","stillbirth",
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


### POPULATION DESCRIPTION ###
##############################

load("N:/durable/projects/ALHE_MR_CHD/mom.RData")
datx<-subset2(mom,"!is.na(mom$cad_grs) & mom$birthyear>1952")
datx<-datx[!duplicated(datx$pid),]
datx$miscarriage_mbrn<-with(datx,ifelse(yearpregmax<1998,NA,miscarriage_mbrn))
datx$gdm<-with(datx,ifelse(yearpregmax<2006,NA,gdm))
datx$childrenmax<-with(datx,ifelse(childrenmax>=4,4,childrenmax))

xxx<-datx[,c("birthyear","agemax","eduyearsmax","bmimax","smokingmax","childrenmax","paritymax",
             "hta_preg","eclampsia","gdm","preterm","preterm_sp","sga","miscarriage_mbrn","stillbirth")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx, method=c("birthyear"=2,"agemax"=2,"bmimax"=2,"smokingmax"=3,"childrenmax"=3,"paritymax"=3,
                                             "subf"=3,"miscarriage_mbrn"=3,"stillbirth"=3,"eclampsia"=3,
                                             "hta_preg"=3,"preterm"=3,"preterm_sp"=3,"sga"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(cbind(all$descr[,1]))
colnames(tab1)<-c("Mothers-All")
write.table(tab1,file="./descriptive/descriptive_mothers.csv",sep=";",col.names=NA)


load("N:/durable/projects/ALHE_MR_CHD/dad.RData")
datx<-subset2(dad,"!is.na(dad$cad_grs) & dad$birthyear>1952")
datx<-datx[!duplicated(datx$pid),]
datx$miscarriage_mbrn<-with(datx,ifelse(yearpregmax<1998,NA,miscarriage_mbrn))
datx$gdm<-with(datx,ifelse(yearpregmax<2006,NA,gdm))
datx$childrenmax<-with(datx,ifelse(childrenmax>=4,4,childrenmax))

xxx<-datx[,c("birthyear","agemax","eduyearsmax","bmimax","smokingmax","childrenmax","paritymax",
             "hta_preg","eclampsia","gdm","preterm","preterm_sp","sga","miscarriage_mbrn","stillbirth")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.,
                               xxx, method=c("birthyear"=2,"agemax"=2,"bmimax"=2,"smokingmax"=3,"childrenmax"=3,"paritymax"=3,
                                             "subf"=3,"miscarriage_mbrn"=3,"stillbirth"=3,"eclampsia"=3,
                                             "hta_preg"=3,"preterm"=3,"preterm_sp"=3,"sga"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(cbind(all$descr[,1]))
colnames(tab1)<-c("Fathers-All")
write.table(tab1,file="./descriptive/descriptive_fathers.csv",sep=";",col.names=NA)


### SELECTION BIAS - INCLUDED vs NON-INCLUDED ###

load("N:/durable/projects/ALHE_MR_CHD/all.RData")
datx<-subset2(dat,"dat$birthyear>1952")
datx$exclusion<-with(datx,ifelse(!is.na(datx$cad_grs),0,1))
datx<-datx[!duplicated(datx$pid),]
datx$miscarriage_mbrn<-with(datx,ifelse(yearpregmax<1998,NA,miscarriage_mbrn))
datx$gdm<-with(datx,ifelse(yearpregmax<2006,NA,gdm))
datx$childrenmax<-with(datx,ifelse(childrenmax>=4,4,childrenmax))

xxx<-datx[,c("exclusion","sex","birthyear","agemax","eduyearsmax","bmimax","smokingmax","childrenmax","paritymax",
             "hta_preg","eclampsia","gdm","preterm","preterm_sp","sga","miscarriage_mbrn","stillbirth")]
all<-NULL
all<-createTable(compareGroups(exclusion~.,
                               xxx, method=c("sex"=3,"birthyear"=2,"agemax"=2,"bmimax"=2,"smokingmax"=3,"childrenmax"=3,"paritymax"=3,
                                             "miscarriage_mbrn"=3,"stillbirth"=3,"eclampsia"=3,
                                             "hta_preg"=3,"preterm"=3,"preterm_sp"=3,"sga"=3,"gdm"=3)),
                 show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Included","Non-included","P-value","N")
write.table(tab1,file="./descriptive/selectionbias.csv",sep=";",col.names=NA)

xxx2<-subset2(xxx,"xxx$sex==1")
all<-NULL
all<-createTable(compareGroups(exclusion~.-sex,
                               xxx2, method=c("birthyear"=2,"agemax"=2,"bmimax"=2,"smokingmax"=3,"childrenmax"=3,"paritymax"=3,
                                              "miscarriage_mbrn"=3,"stillbirth"=3,"eclampsia"=3,
                                              "hta_preg"=3,"preterm"=3,"preterm_sp"=3,"sga"=3,"gdm"=3)),
                 show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)
tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Included","Non-included","P-value","N")
write.table(tab1,file="./descriptive/selectionbias_mothers.csv",sep=";",col.names=NA)

xxx2<-subset2(xxx,"xxx$sex==0")
all<-NULL
all<-createTable(compareGroups(exclusion~.-sex,
                               xxx2, method=c("birthyear"=2,"agemax"=2,"bmimax"=2,"smokingmax"=3,"childrenmax"=3,"paritymax"=3,
                                              "miscarriage_mbrn"=3,"stillbirth"=3,"eclampsia"=3,
                                              "hta_preg"=3,"preterm"=3,"preterm_sp"=3,"sga"=3,"gdm"=3)),
                 show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)
tab1<-NULL
tab1<-as.data.frame(all$descr)
colnames(tab1)<-c("Included","Non-included","P-value","N")
write.table(tab1,file="./descriptive/selectionbias_fathers.csv",sep=";",col.names=NA)


### ROBUSTNESS OF CAD GRS ###
#############################

load("N:/durable/projects/ALHE_MR_CHD/all.RData")
library(pROC)

# Number of SNPs used for each GRS (CAD: n=146)

vars00<-c("cad_grs","cad_grs","cad_grs")
vars01<-c("cad_grs_z","cad_grs_z","cad_grs_z")
vars02<-c("chd_lifetime","chd_lifetime","chd_lifetime")
vars03<-c(0,1,2)
vars04<-c("women_","men_","all_")
vars05<-c("CHD (odds ratio)","CHD (odds ratio)","CHD (odds ratio)"
vars06<-c("CHD-GRS, women","CHD-GRS, men","CHD-GRS, all")

tab<-NULL

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat$sex!=vars03[i] & !is.na(dat[,vars00[i]]) & !is.na(dat[,vars02[i]])")
  datx[,vars01[i]]<-as.numeric(with(datx,scale(datx[,vars00[i]])))
  mod01<-glm(formula=as.factor(datx[,vars02[i]])~datx[,vars01[i]],
             data=datx, family="binomial")
  estimate<-as.numeric(summary(mod01)$coefficients[2,1])
  se<-as.numeric(summary(mod01)$coefficients[2,2])
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(estimate))
  ic95a<-round(exp(estimate-(z*se)),3)
  ic95b<-round(exp(estimate+(z*se)),3)
  coef<-ic_guapa(hr,ic95a,ic95b)
  llrtest<-round(as.numeric(anova(mod01,test=c("F"))[2,5]),0) # log-likelihood ratio test (~F in binomial GLM)
  pseudor2<-paste("'",guapa(NagelkerkeR2(mod01)$R2*100),sep="") # pseudoR2 by the Nagelkerke method
  n_snps<-c("146")
  datx$prob<-predict.glm(mod01,type=c("response"))
  roc_auc<-guapa(roc(as.factor(datx[,vars02[i]])~prob,data=datx)$auc)
  mean_grs<-paste(guapa(mean(datx[,vars00[i]],na.rm=TRUE))," (",guapa(sd(datx[,vars00[i]],na.rm=TRUE)),")",sep="")
  size<-dim(datx)[1]
  cases<-paste(table(datx[,vars02[i]])[2]," (",guapa(table(datx[,vars02[i]])[2]/dim(datx)[1]*100)," %)",sep="")
  
  aaa<-datx[,vars00[i]]
  mod_lin<-gam(formula=as.factor(datx[,vars02[i]])~aaa,
               data=datx, family="binomial")
  mod_nlin<-gam(formula=as.factor(datx[,vars02[i]])~bs(aaa,df=4),
                data=datx, family="binomial")
  p_lrtest<-pval_guapa(lrtest(mod_lin,mod_nlin)[2,5])
  
  aaa<-datx[,vars00[i]]
  minim<-mean(aaa,na.rm=TRUE)
  mod01<-gam(formula=as.factor(datx[,vars02[i]])~bs(aaa,df=4),
             data=datx, family="binomial")
  ptemp<-termplot(mod01,term=1,se=TRUE,plot=FALSE)
  temp<-ptemp$aaa
  value<-closest(temp$x,minim)
  #min_val<-temp$x[which(temp$y==min(temp$y,na.rm=TRUE))]
  center<-with(temp, y[x==value])
  z<-qnorm(1-0.05/2)
  ytemp<-temp$y+outer(temp$se,c(0,-z,z),'*')
  min_val<-guapa(temp$x[which(temp$y==min(temp$y,na.rm=TRUE))])
  ci<-exp(ytemp-center)
  name<-paste("./descriptive/robustness_",vars04[i],vars02[i],".jpg",sep="")
  labely<-vars05[i]
  
  plot.data<-as.data.frame(cbind(temp$x,ci))
  colnames(plot.data)<-c("xval","yest","lci","uci")
  figure<-ggplot() +
    geom_histogram(data=datx, aes(x=datx[,vars00[i]], y=(..density../max(..density..,na.rm=TRUE))*3.75), 
                   bins=20, color="grey80", fill="grey90") +
    #scale_x_continuous(limits = c(min(dat2[,vars01[i]],na.rm=TRUE),max(dat2[,vars01[i]],na.rm=TRUE))) +
    scale_y_continuous(limits = c(0,4), sec.axis = sec_axis(~., name = "Relative density")) +
    geom_hline(aes(yintercept=1), data=plot.data, colour="black", linetype=2) + 
    geom_line(aes(x=xval, y=yest), data=plot.data, color="black") + 
    geom_line(aes(x=xval, y=lci), data=plot.data, color="grey35") + 
    geom_line(aes(x=xval ,y=uci), data=plot.data, color="grey35") + 
    theme_bw() +
    labs(x=vars06[i],y=labely) +
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
  
  tab<-rbind(tab,cbind(size,cases,n_snps,mean_grs,coef,llrtest,roc_auc,pseudor2,p_lrtest))
}

rownames(tab)<-paste(vars04,vars02,sep="")
write.table(tab,file="./descriptive/robustness_grs.csv",sep=";",col.names=NA)


### LOGISTIC REGRESSION: LINEAR ASSOCIATIONS ###
################################################

load("N:/durable/projects/ALHE_MR_CHD/all.RData")
dat$adj<-0
dat$excl_misc<-with(dat,ifelse(yearpregmax<1998,1,0))
dat$excl_gdm<-with(dat,ifelse(yearpregmax<2006,1,0))

vars01<-c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0)
vars02<-c("miscarriage_mbrn","miscarriage_mbrn","stillbirth","stillbirth",
          "eclampsia","eclampsia","hta_preg","hta_preg",
          "preterm","preterm","preterm_sp","preterm_sp","sga","sga","gdm","gdm")
vars03<-c("excl_misc","excl_misc","adj","adj",
          "hta_chronic","hta_chronic","hta_chronic","hta_chronic",
          "adj","adj","adj","adj","adj","adj","excl_gdm","excl_gdm")
vars08<-c("pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1")
vars09<-c("pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2")
vars10<-c("pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3")
vars11<-c("pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4")
vars12<-c("pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5")
vars13<-c("pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6")
vars14<-c("pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7")
vars15<-c("pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8")
vars16<-c("pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9")
vars17<-c("pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10")
vars28<-c("CAD-GRS - Miscarriage-MBRN, mothers","CAD-GRS - Miscarriage-MBRN, fathers",
          "CAD-GRS - Stillbirth, mothers","CAD-GRS - Stillbirth, fathers",
          "CAD-GRS - Pre-/Eclampsia, mothers","CAD-GRS - Pre-/Eclampsia, fathers",
          "CAD-GRS - HT in pregnancy, mothers","CAD-GRS - HT in pregnancy, fathers",
          "CAD-GRS - Preterm, mothers","CAD-GRS - Preterm, fathers",
          "CAD-GRS - Spontaneous preterm, mothers","CAD-GRS - Spontaneous preterm, fathers",
          "CAD-GRS - SGA, mothers","CAD-GRS - SGA, fathers",
          "CAD-GRS - Gestat. diabetes, mothers","CAD-GRS - Gestat. diabetes, fathers")

tab<-NULL
z<-qnorm(1-0.05/2)

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat$sex==vars01[i] & dat[,vars03[i]]==0 & !is.na(dat$cad_grs) & !is.na(dat[,vars02[i]]) & dat$birthyear>1952")
  datx$cad_grs_z<-as.numeric(with(datx,scale(datx$cad_grs)))
  sample<-dim(datx)[1]
  
  mod01<-glm(formula=as.factor(datx[,vars02[i]])~cad_grs_z,
             data=datx, family="binomial")
  estimate01<-as.numeric(summary(mod01)$coefficients[2,1])
  se01<-as.numeric(summary(mod01)$coefficients[2,2])
  coef01<-risk_se_ic_guapa(estimate01,se01)
  pval01<-pval_guapa(as.numeric(summary(mod01)$coefficients[2,4]))
  
  mod02<-glm(formula=as.factor(datx[,vars02[i]])~cad_grs_z
             +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]+datx[,vars12[i]]
             +datx[,vars13[i]]+datx[,vars14[i]]+datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]]+as.factor(batch),
             data=datx, family="binomial")
  estimate02<-as.numeric(summary(mod02)$coefficients[2,1])
  se02<-as.numeric(summary(mod02)$coefficients[2,2])
  coef02<-risk_se_ic_guapa(estimate02,se02)
  pval02<-pval_guapa(as.numeric(summary(mod02)$coefficients[2,4]))
  or02<-round(exp(estimate02),5)
  ic95lo02<-round(exp(estimate02-(z*se02)),5)
  ic95hi02<-round(exp(estimate02+(z*se02)),5)
  
  aaa<-datx$cad_grs_z
  mod_base<-glm(formula=as.factor(datx[,vars02[i]])~datx[,vars08[i]]
                +datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]+datx[,vars12[i]]
                +datx[,vars13[i]]+datx[,vars14[i]]+datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]]+as.factor(batch),
                data=datx, family="binomial")
  mod_lin<-glm(formula=as.factor(datx[,vars02[i]])~aaa
               +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]+datx[,vars12[i]]
               +datx[,vars13[i]]+datx[,vars14[i]]+datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]]+as.factor(batch),
               data=datx, family="binomial")
  mod_nlin<-glm(formula=as.factor(datx[,vars02[i]])~bs(aaa,df=4)
                +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]+datx[,vars12[i]]
                +datx[,vars13[i]]+datx[,vars14[i]]+datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]]+as.factor(batch),
                data=datx, family="binomial")
  p_lrtest02<-pval_guapa(lrtest(mod_lin,mod_nlin)[2,5])
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02,p_lrtest02,
                       estimate02,se02,or02,ic95lo02,ic95hi02,sample))
}

colnames(tab)<-c("coef (raw)","pval (raw)","coef (adj.)","pval (adj.)","p_lrtest",
                 "beta","se","or","ic95lo","ic95hi","sample")
rownames(tab)<-vars28
write.table(tab,file="./results/linear.csv",sep=";",col.names=NA)


### ANALYSES ADJUSTED FOR GRS OF THE PARTNER ###
################################################

load("N:/durable/projects/ALHE_MR_CHD/mom.RData")
load("N:/durable/projects/ALHE_MR_CHD/dad.RData")
mom2<-subset2(mom,"!is.na(mom$partner)")
dad2<-subset2(dad,"!is.na(dad$partner)")
dad2<-dad2[,c("pid","cad_grs")]
dad2<-rename.vars(dad2,from=c("pid","cad_grs"),to=c("partner","cad_grs_partner"))
mom_merge<-merge2(mom2,dad2,by.id=c("partner"),all.x=TRUE,sort=FALSE)

load("N:/durable/projects/ALHE_MR_CHD/mom.RData")
load("N:/durable/projects/ALHE_MR_CHD/dad.RData")
dad2<-subset2(dad,"!is.na(dad$partner)")
mom2<-subset2(mom,"!is.na(mom$partner)")
mom2<-mom2[,c("pid","cad_grs")]
mom2<-rename.vars(mom2,from=c("pid","cad_grs"),to=c("partner","cad_grs_partner"))
dad_merge<-merge2(dad2,mom2,by.id=c("partner"),all.x=TRUE,sort=FALSE)

dat<-rbind(mom_merge,dad_merge)
dat$adj<-0
dat$excl_misc<-with(dat,ifelse(yearpregmax<1998,1,0))
dat$excl_gdm<-with(dat,ifelse(yearpregmax<2006,1,0))

vars01<-c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0)
vars02<-c("miscarriage_mbrn","miscarriage_mbrn","stillbirth","stillbirth",
          "eclampsia","eclampsia","hta_preg","hta_preg",
          "preterm","preterm","preterm_sp","preterm_sp","sga","sga","gdm","gdm")
vars03<-c("excl_misc","excl_misc","adj","adj",
          "hta_chronic","hta_chronic","hta_chronic","hta_chronic",
          "adj","adj","adj","adj","adj","adj","excl_gdm","excl_gdm")
vars08<-c("pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1","pc1")
vars09<-c("pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2","pc2")
vars10<-c("pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3","pc3")
vars11<-c("pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4","pc4")
vars12<-c("pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5","pc5")
vars13<-c("pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6","pc6")
vars14<-c("pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7","pc7")
vars15<-c("pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8","pc8")
vars16<-c("pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9","pc9")
vars17<-c("pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10","pc10")
vars28<-c("CAD-GRS - Miscarriage-MBRN, mothers","CAD-GRS - Miscarriage-MBRN, fathers",
          "CAD-GRS - Stillbirth, mothers","CAD-GRS - Stillbirth, fathers",
          "CAD-GRS - Pre-/Eclampsia, mothers","CAD-GRS - Pre-/Eclampsia, fathers",
          "CAD-GRS - HT in pregnancy, mothers","CAD-GRS - HT in pregnancy, fathers",
          "CAD-GRS - Preterm, mothers","CAD-GRS - Preterm, fathers",
          "CAD-GRS - Spontaneous preterm, mothers","CAD-GRS - Spontaneous preterm, fathers",
          "CAD-GRS - SGA, mothers","CAD-GRS - SGA, fathers",
          "CAD-GRS - Gestat. diabetes, mothers","CAD-GRS - Gestat. diabetes, fathers")

tab<-NULL
z<-qnorm(1-0.05/2)

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat$sex==vars01[i] & dat[,vars03[i]]==0 & !is.na(dat$cad_grs) & !is.na(dat[,vars02[i]]) & dat$birthyear>1952")
  datx$cad_grs_z<-as.numeric(with(datx,scale(datx$cad_grs)))
  datx$cad_grs_partner_z<-as.numeric(with(datx,scale(datx$cad_grs_partner)))
  sample<-dim(datx)[1]
  
  mod02<-glm(formula=as.factor(datx[,vars02[i]])~cad_grs_z
             +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]+datx[,vars12[i]]
             +datx[,vars13[i]]+datx[,vars14[i]]+datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]]+as.factor(batch)+cad_grs_partner_z,
             data=datx, family="binomial")
  estimate02<-as.numeric(summary(mod02)$coefficients[2,1])
  se02<-as.numeric(summary(mod02)$coefficients[2,2])
  coef02<-risk_se_ic_guapa(estimate02,se02)
  pval02<-pval_guapa(as.numeric(summary(mod02)$coefficients[2,4]))
  or02<-round(exp(estimate02),5)
  ic95lo02<-round(exp(estimate02-(z*se02)),5)
  ic95hi02<-round(exp(estimate02+(z*se02)),5)
  
  aaa<-datx$cad_grs_z
  mod_base<-glm(formula=as.factor(datx[,vars02[i]])~datx[,vars08[i]]
                +datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]+datx[,vars12[i]]
                +datx[,vars13[i]]+datx[,vars14[i]]+datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]]+as.factor(batch)+cad_grs_partner_z,
                data=datx, family="binomial")
  mod_lin<-glm(formula=as.factor(datx[,vars02[i]])~aaa
               +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]+datx[,vars12[i]]
               +datx[,vars13[i]]+datx[,vars14[i]]+datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]]+as.factor(batch)+cad_grs_partner_z,
               data=datx, family="binomial")
  mod_nlin<-glm(formula=as.factor(datx[,vars02[i]])~bs(aaa,df=4)
                +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]+datx[,vars12[i]]
                +datx[,vars13[i]]+datx[,vars14[i]]+datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]]+as.factor(batch)+cad_grs_partner_z,
                data=datx, family="binomial")
  p_lrtest02<-pval_guapa(lrtest(mod_lin,mod_nlin)[2,5])
  
  tab<-rbind(tab,cbind(coef02,pval02,p_lrtest02,
                       estimate02,se02,or02,ic95lo02,ic95hi02,sample))
}

colnames(tab)<-c("coef (adj.)","pval (adj.)","p_lrtest",
                 "beta","se","or","ic95lo","ic95hi","sample")
rownames(tab)<-vars28
write.table(tab,file="./results/linear_adj_grs_partner.csv",sep=";",col.names=NA)
