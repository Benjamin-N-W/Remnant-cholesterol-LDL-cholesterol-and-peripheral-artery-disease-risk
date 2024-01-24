#start
rm(list=ls())
detach(cgps)
detach(CcHS)

#set language
Sys.setenv(LANG = "en")

#Setting working directory
setwd("H:/Projekt remnant cholesterol/Artikel 4 - PAD and metabolomics/R working folder - PAD and metabolomics")

# Loading packages --------------------------------------------------------
library(Hmisc)
library(sjlabelled)
library(matrixStats)
library(haven)

library(pastecs)
library(xtable)
library(dplyr)
library(ggplot2)
library(ggpubr)

#Survival analysis

library(survival)
library(cmprsk)
library(survminer)
library(CMAverse)

#Plotting
library(eoffice)
library(stringr)

#Using ggplot
options(grType='RDevice')
library(broom)
library(flextable)
library(webshot)
library(magick)
library(grid)
library(cowplot)
library(CMAverse)

library(gridExtra)
library(gridGraphics)
library(grid)
library(haven)

library(rms)


#Function for P values in string format
signif.p<-function(x, exponent=F){
  
  result<-c()
  
  for(i in 1:length(x)){
    
    if(exponent==T){
      if (x[i]==0){
        result[i]<-as.character("-16")
        
        } else {
      if(log10(x[i])>(-4)){ 
        result[i]<-NA
      }  else {
        
        result[i]<-as.character(floor(log10(signif(x[i],1))))
      } 
    } 
    }else if (exponent==F) {
      
      
      if(log10(x[i])>(-4)){
        result[i]<-as.character(format(signif(x[i],2),scientific = F)
        )
      } else {
        
        if (x[i]==0){
          q<-"<2"
          
        } else {
        
        q<-signif(x[i],1)*10^-round(log10(x[i]))
       
        if (q<1) {
          
          q<-q*10
        }
        }
        
        result[i]<-paste(as.character(q),"x 10")
      }
    } 
  }
  return(result)
}



# Setting up data ---------------------------------------------------------



#load stata file

cgps<-read_dta("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstr?m/Dataset_obush.dta")
cgps<- cgps[-c(69:72)] 


op_koder<-read_dta("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstr?m/Fra_anders/PAD_op_clean_CGPS.dta")

pad_diagnr<-read_dta("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstrom/PAD_obush.dta")
pad_bdiagnr<-read_dta("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstrom/PAD_obush_medBdiag.dta")
lpa<-read_dta("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstr?m/lpa_imputed_validated.dta")
cgps_stroke<-read_dta("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstr?m/stroke_obush.dta")
cgps_hopkins<-read.csv("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstr?m/CGPS_ldl_hopkins_17022021.csv",header=TRUE)
cgps_metabolomics<-read.csv("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstr?m/Fra_shoaib/metabolomics_30_Mar_2020_obush1.csv")
cgps_metabolomics<- cgps_metabolomics[-2] 
cgps_apob<-read.csv2("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstr?m/Fra_anders/apob_mgdl_OBUSH.csv",header=TRUE)
cgps_koag<-read.csv("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstr?m/Fra_anders/OBUSH_aptt_koag_inr_fibrinog.csv",header=TRUE)


cgps<-merge(cgps,lpa,by="obushnr", all=TRUE)

cgps<- merge(cgps,cgps_stroke,by="obushnr", all=TRUE)

cgps<- merge(cgps,pad_bdiagnr,by="obushnr", all=TRUE)
cgps<- merge(cgps,op_koder,by="obushnr", all=TRUE)
cgps<- merge(cgps,cgps_hopkins,by="obushnr", all=TRUE)
cgps<- merge(cgps,cgps_metabolomics,by="obushnr", all=TRUE)
cgps<- merge(cgps,cgps_apob,by="obushnr", all=TRUE)
cgps<- merge(cgps,cgps_koag,by="obushnr", all=TRUE)


rm("cgps_metabolomics", "cgps_stroke", "lpa", "cgps_apob","cgps_koag", "op_koder","pad_bdiagnr","p_diagnr","cgps_hopkins")


# Converting data ---------------------------------------------------------

cgps$apob_mgdl<-as.numeric(cgps$apob_mgdl)

#As factor
cgps$P1methyl2pyrrolidone<-factor(cgps$P1methyl2pyrrolidone)


#Calculate Anbkle-brachial index
cgps$abileft<-cgps$OBJ22/cgps$OBJ18

cgps$abiright<-cgps$OBJ20/cgps$OBJ18

cgps$abilowest<-pmin(cgps$abileft,cgps$abiright)

cgps$abi0.9<-cut(cgps$abilowest, breaks=c(-Inf, 0.9 , Inf), right=FALSE, labels=c("<0.9",">0.9"), include.lowest = TRUE)
cgps$abi0.9<- relevel(cgps$abi0.9, ">0.9")

#NA values are converted to 0
cgps$PAD[is.na(cgps$PAD)] <- 0

# Removing renal artery stenosis (n=14)
cgps$PAD[cgps$c_diag=="44019"]<-0
cgps$PAD[cgps$c_diag=="DI701"]<-0

#Removing other and unspecified atherosclerosis (n=115)
cgps$PAD[cgps$c_diag=="44599"]<-0
cgps$PAD[cgps$c_diag=="44099"]<-0
cgps$PAD[cgps$c_diag=="44039"]<-0
cgps$PAD[cgps$c_diag=="DI708"]<-0
cgps$PAD[cgps$c_diag=="DI709"]<-0

#Questionnaire PAD
cgps$v25d<-with(cgps, ifelse(cgps$v25b==0 & is.na(v25d)==TRUE| v25d==0 , FALSE, TRUE)) 

#Creating prevalent PAD variable
cgps$PAD_pre<-with(cgps, ifelse(PAD==1 & PAD_dato<usdato| v25d==TRUE |abi0.9=="<0.9", TRUE, FALSE)) 
cgps$PAD_pre[is.na(cgps$PAD_pre)] <- FALSE

#Registry diagnosis at baseline
cgps$PAD_pre_diag<-with(cgps, ifelse(PAD==1 & PAD_dato<usdato, TRUE, FALSE)) 

#Remnant cholesterol variable
cgps$remncholcalc<-cgps$chol-cgps$ldl-cgps$hdl


#filling values for those with missing values for ldl but values for triglycerides. N=282
cgps$remncholcalc[is.na(cgps$remncholcalc)]<-(cgps$trig*5/11)[is.na(cgps$remncholcalc)]

#ropping observations with NA values for remnant cholesterol (mostly no blood samples taken). N=437
cgps<-cgps[!is.na(cgps$remncholcalc), ]


#96 subjects with remncholcalc <0.00001
cgps$remncholcalc[cgps$remncholcalc<0.00001]<-(cgps$trig*5/11)[cgps$remncholcalc<0.00001]

#prevalent diabetes type 1


#Creating variable for triglycerides >10
cgps$trig10<-cut(cgps$trig, breaks=c(-Inf, 10 , Inf), right=FALSE, labels=c("<10",">10"), include.lowest = TRUE)

#prevalent diabetes
cgps$dmall_pre<-with(cgps, ifelse(dmiidato<=usdato  & dmii==1 |dmidato<=usdato  & dmi==1 | glucose>11 | v27==1 | m11==1|m12==1, 1, 0)) 
cgps$dmi_pre<-with(cgps, ifelse(dmidato<=usdato  & dmi==1, TRUE, FALSE))
cgps$dmii_pre<-with(cgps, ifelse(dmiidato<=usdato  & dmii==1, TRUE, FALSE)) 


#Prevalent ami
cgps$ihd_pre<-with(cgps, ifelse(ihddato<=usdato  & ihd==1 , 1, 0)) 
cgps$ihd_pre<-as.logical(cgps$ihd_pre)

#Prevalent ischemic stroke
cgps$str_pre<-with(cgps, ifelse(isdato_uval<=usdato  & is_uval==1 , 1, 0)) 
cgps$str_pre<-as.logical(cgps$str_pre)

#Prevalent ami or stroke (ASCVD)
cgps$ascvd_pre<- with(cgps, ifelse(str_pre==TRUE | ihd_pre==TRUE , TRUE, FALSE)) 

#Start age
cgps$startage<-(cgps$usdato-cgps$fodedato)/365.25
cgps$startage<-as.numeric(cgps$startage)

#Stop age PAD
cgps$PAD_dato<-pmin(cgps$PAD_dato,cgps$doddato, na.rm = TRUE)
cgps$stopagePAD<-(cgps$PAD_dato-cgps$fodedato)/365.25
cgps$stopagePAD<-as.numeric(cgps$stopagePAD)

#Stop age ami
cgps$stopageIHD<-(cgps$amidato-cgps$fodedato)/365.25
cgps$stopageIHD<-as.numeric(cgps$stopageIHD)

#Stop age Ischemic stroke
cgps$stopageIS<-(cgps$isdato_uval-cgps$fodedato)/365.25
cgps$stopageIS<-as.numeric(cgps$stopageIS)

#Stop age mortality
cgps$stopagemort<-(cgps$doddato-cgps$fodedato)/365.25
cgps$stopagemort<-as.numeric(cgps$stopagemort)

#Birth date adjustment variable
cgps$birthdate<-as.numeric(cgps$fodedato)/365.25+1970


#Dropping 271 rows with fodedato=NA
cgps<-cgps[!is.na(cgps$fodedato), ]

#Dropping 3 rows with date of death before examination
cgps<-cgps[!cgps$doddato-cgps$usdato<0, ]

#Dropping subjects duplicated in cgps
cgps<-cgps[!cgps$duplicate_within_cgps==1, ]
cgps<- cgps[!duplicated(cgps$obushnr)==TRUE, ]

#Dropping subjects of another ethnicity (n=255)
cgps<-cgps[cgps$herkomst==1, ]
cgps<-cgps[cgps$oprindelse==1, ]

# Blanks are converted to NA values
cgps$diagnr1[cgps$diagnr1==""] <- NA
cgps$diagnr2[cgps$diagnr2==""] <- NA
cgps$diagnr3[cgps$diagnr3==""] <- NA
cgps$diagnr4[cgps$diagnr4==""] <- NA
cgps$diagnr5[cgps$diagnr5==""] <- NA
cgps$diagnr6[cgps$diagnr6==""] <- NA
cgps$diagnr7[cgps$diagnr7==""] <- NA
cgps$diagnr8[cgps$diagnr8==""] <- NA
cgps$diagnr9[cgps$diagnr9==""] <- NA


#Sex varaible is recoded to factor
cgps$sex <- factor(cgps$sex, 0:1, labels=c("Female", "Male"))

#Recoding  variables
cgps$status_smoking<-factor(cgps$status_smoking, levels=c(0,2,1,NULL),labels=c("Never smoker", "Ex smoker", "Current smoker"), ordered=TRUE)
cgps$v75<-factor(cgps$v75, levels=c(1,2,3,4, NULL),labels=c("Light exercise <2 hours/week", "Light exercise 2-4 hours/week", "Light exercise >4 hours/week", "Moderate-heavy exercise >4 hours/week"), ordered=TRUE)

#Converting endpoint variables to logical format
cgps$dod<-as.logical(cgps$dod)
cgps$PAD<-as.logical(cgps$PAD)
cgps$ami<-as.logical(cgps$ami)
cgps$is_uval<-as.logical(cgps$is_uval)

#Converting covariates to logical format
cgps$m5<-as.logical(cgps$m5)
cgps$m2a<-as.logical(cgps$m2a)
cgps$v25d<-as.factor(cgps$v25d)

#Creating multi-state event for competing risk
cgps$dod_comprisk<-cgps$dod
cgps$dod_comprisk[grepl("I",cgps$diagnr1)==TRUE]<-2
cgps$dod_comprisk[grepl("I",cgps$diagnr2)==TRUE]<-2
cgps$dod_comprisk[grepl("I",cgps$diagnr3)==TRUE]<-2
cgps$dod_comprisk[grepl("C",cgps$diagnr1)==TRUE]<-3
cgps$dod_comprisk[grepl("C",cgps$diagnr2)==TRUE]<-3 
cgps$dod_comprisk[grepl("C",cgps$diagnr3)==TRUE]<-3 
cgps$dod_comprisk<-factor(cgps$dod_comprisk, labels = c("Censored", "Other", "CVD", "Cancer"))


#PAD competing tisk variable
cgps$PAD_comprisk<-as.integer(cgps$PAD)
cgps$PAD_comprisk[cgps$dod==TRUE]<-1
cgps$PAD_comprisk[cgps$st131218==70]<-1
cgps$PAD_comprisk[cgps$st131218==80]<-1
cgps$PAD_comprisk[cgps$PAD==TRUE]<-2
cgps$PAD_comprisk<-factor(cgps$PAD_comprisk, labels = c("Censored", "Death", "PAD"))


library(mosaic)

#Current smoking status variable
cgps$current_smoking<-with(cgps, ifelse(status_smoking=="Current smoker", 1, 0)) 
cgps$current_smoking[is.na(cgps$status_smoking)]<-NA

#Ever smoking status variable
cgps$ever_smoking<-with(cgps, ifelse(status_smoking=="Ex smoker" | status_smoking=="Current smoker", TRUE, FALSE)) 
cgps$ever_smoking[is.na(cgps$status_smoking)]<-NA


cgps$cum_smoking[cgps$status_smoking=="Never smoker"]<-0

#NA for influential observations of smoking and systolic BP (>99.9th percentile)
cgps$systolic[cgps$systolic>200]<-NA
cgps$systolic[cgps$systolic<90]<-NA
cgps$diastolic[cgps$diastolic<55]<-NA
cgps$diastolic[cgps$diastolic>120]<-NA

#Creating hypertension variable
cgps$hypertens<-with(cgps, ifelse(systolic>=140 | diastolic>=90 | m2a==TRUE, 1, 0)) 

#Creating categories of remnant cholesterol
cgps$remncholcalc_cut4<-cut(cgps$remncholcalc, breaks=c(min(cgps$remncholcalc, na.rm=TRUE), quantile(cgps$remncholcalc, 0.25, na.rm=TRUE), quantile(cgps$remncholcalc, 0.5, na.rm=TRUE),
                                                        quantile(cgps$remncholcalc, 0.75, na.rm=TRUE), max(cgps$remncholcalc, na.rm=TRUE)), right=FALSE, include.lowest = TRUE, labels = c("1st quartile", "2nd quartile","3rd quartile","4th quartile"))


cgps$remncholcalc_cut05<-cut(cgps$remncholcalc, breaks=c(min(cgps$remncholcalc, na.rm=TRUE), 0.5, 1,
                                                         1.5, max(cgps$remncholcalc, na.rm=TRUE)), right=FALSE, include.lowest = TRUE, labels = c("<0.5", "0.5-1.0","1.0-1.5",">1.5"), ordered_result = TRUE)

#Creating categorical variables of continous variables
cgps$alder_cat<-cut2(cgps$alder, g=4)
cgps$ldl_cat<-cut2(cgps$ldl, g=4)
cgps$chol_cat<-cut2(cgps$chol, g=4)
cgps$apoa1_cat<-cut2(cgps$apoa1, g=4)
cgps$systolic_cat<-cut2(cgps$systolic, g=4)
cgps$diastolic_cat<-cut2(cgps$diastolic, g=4)
cgps$bmi_cat<-cut2(cgps$bmi, g=4)
cgps$glucose_cat<-cut2(cgps$glucose, g=4)
cgps$hdl_cat<-cut2(cgps$hdl, g=4)
cgps$trig_cat<-cut2(cgps$trig, g=4)
cgps$cum_smoking_cat<-cut2(cgps$cum_smoking, g=4)
cgps$alcohol_cat<-cut2(cgps$alcohol, g=4)
cgps$hscrp_cat<-cut2(cgps$hscrp, g=4)
cgps$creatini_cat<-cut2(cgps$creatini, g=4)
cgps$hdl_cat<-cut2(cgps$hdl, g=4)


#Categorical varaibles by biological cut-offs
cgps$bmi_cat25<-cut2(cgps$bmi, cuts=25)
cgps$alder_cat25<-cut2(cgps$alder, cuts=70)
cgps$ldl_cat25<-cut2(cgps$ldl, cuts=2.5)
cgps$birthdate_cat25<-cut2(cgps$birthdate, cuts=median(cgps$birthdate))
cgps$systolic_cat25<-as.numeric(cut2(cgps$systolic, cuts=140))
cgps$systolic_cat25[cgps$systolic_cat25==1]<-0
cgps$systolic_cat25[cgps$systolic_cat25==2]<-1
cgps$systolic_cat25<-as.logical(cgps$systolic_cat25)


# Imputing missing values, first using aregImpute for multiple imp --------

#Multiple imputation 
#impute_arg <- aregImpute(~ alder+sex+ihd_pre+status_smoking+current_smoking+ldl+hypertens+diastolic+ascvd_pre+cum_smoking+bmi+systolic+m2a+birthdate, data = cgps, n.impute = 1)

#cgps_imp<-cgps
#tbl_imp_areg <- impute.transcan(impute_arg,imputation = 1,  data = cgps,  list.out = TRUE, pr = FALSE,  check = FALSE)
#cgps_imp[names(tbl_imp_areg)] <- tbl_imp_areg
#saveRDS(cgps_imp, file="L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstr?m/cgps_imputed_metabolomics_PAD.Rda")

cgps_imp<-readRDS("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstrom/cgps_imputed_metabolomics_PAD.Rda")

cgps<-cgps_imp[!is.na(cgps_imp$fodedato), ]

cgps_clti<-read_dta("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstrom/Fra Peter/MALE_tilBenjamin_171023.dta")

cgps_clti<-rename(cgps_clti, CLTI  = Amp_Proc_smal,CLTI_dato  = Amp_Proc_smal_dato, CLTI_op=c_opr)

cgps<- merge(cgps,cgps_clti,by="obushnr", all=TRUE)

cgps$CLTI<-ifelse(is.na(cgps$CLTI),F,T)
cgps$CLTI_dato<-with(cgps,pmin(doddato,CLTI_dato,na.rm=T))
cgps<-subset(cgps,!is.na(startage))

#Adding gangrene and rest pain
cgps[with(cgps,(c_diag=="DI739C"|c_diag=="DI702A") & !is.na(c_diag)),"CLTI"]<-T

cgps[with(cgps,(c_diag=="DI739C"|c_diag=="DI702A") & !is.na(c_diag)& PAD_dato<CLTI_dato) ,"CLTI_op"]<-cgps[with(cgps,(c_diag=="DI739C"|c_diag=="DI702A") & !is.na(c_diag)& PAD_dato<CLTI_dato) ,"c_diag"]

cgps[with(cgps,(c_diag=="DI739C"|c_diag=="DI702A") & !is.na(c_diag)),"CLTI_dato"]<-pmin(cgps[with(cgps,(c_diag=="DI739C"|c_diag=="DI702A") & !is.na(c_diag)),"PAD_dato" ],
                                                                                         cgps[with(cgps,(c_diag=="DI739C"|c_diag=="DI702A") & !is.na(c_diag)),"CLTI_dato" ])

rm("cgps_imp", "tbl_imp_areg", "impute_arg","cgps_clti")



# Converting data - post imputation ---------------------------------------

#Ddropping individuals without statins and without NA LDL and NA apoB
cgps<-cgps[cgps$m5==0 & !is.na(cgps$m5) & !is.na(cgps$ldl), ]

cgps<-cgps[!is.na(cgps$apob_mgdl), ]

#Changing PAD definition based on operation codes
cgps$PAD<-ifelse(cgps$PAD==T | cgps$CLTI==T,T,F)
cgps[with(cgps,  CLTI_dato<PAD_dato) ,"c_diag"]<-cgps[with(cgps, CLTI_dato<PAD_dato)  ,"CLTI_op"]
cgps$PAD_dato<-pmin(cgps$PAD_dato,cgps$CLTI_dato)

# Removing 0s (for later log-conversion) and converting Particle count to mikromol/L (10^-6)
for (i in c("SL","ML","LL","I","XSVL","SVL","MVL","LVL","XLVL","XXLVL")){
  
  cgps[ ,paste0(i,"DLP")]<-cgps[ ,paste0(i,"DLP")]*1e9
  
  cgps[cgps[ ,paste0(i,"DLP")]==0  & !is.na(cgps[ ,paste0(i,"DLP")]) ,paste0(i,"DLP")]<-min(cgps[cgps[ ,paste0(i,"DLP")]>0,paste0(i,"DLP")],na.rm=T)
  cgps[cgps[ ,paste0(i,"DLC")]==0  & !is.na(cgps[ ,paste0(i,"DLC")]),paste0(i,"DLC")]<-min(cgps[cgps[ ,paste0(i,"DLC")]>0,paste0(i,"DLC")],na.rm=T)
  cgps[cgps[ ,paste0(i,"DLTG")]==0  & !is.na(cgps[ ,paste0(i,"DLTG")]),paste0(i,"DLTG")]<-min(cgps[cgps[ ,paste0(i,"DLTG")]>0,paste0(i,"DLTG")],na.rm=T)
  cgps[cgps[ ,paste0(i,"DLPL")]==0  & !is.na(cgps[ ,paste0(i,"DLPL")]),paste0(i,"DLPL")]<-min(cgps[cgps[ ,paste0(i,"DLPL")]>0,paste0(i,"DLPL")],na.rm=T)
}

cgps$VLDLPL<-(cgps$XXLVLDLPL+cgps$XLVLDLPL+cgps$LVLDLPL+cgps$MVLDLPL+cgps$SVLDLPL+cgps$XSVLDLPL)

#Creating total LDL and total VLDL particle number 
cgps$LDLP<-(cgps$SLDLP+cgps$MLDLP+cgps$LLDLP+cgps$IDLP)
cgps$VLDLP<-(cgps$XXLVLDLP+cgps$XLVLDLP+cgps$LVLDLP+cgps$MVLDLP+cgps$SVLDLP+cgps$XSVLDLP)
cgps$totalP<-(cgps$LDLP+cgps$VLDLP)

#Cholesterol content
cgps$LDLC<-(cgps$SLDLC+cgps$MLDLC+cgps$LLDLC+cgps$IDLC)
cgps$VLDLC<-(cgps$XXLVLDLC+cgps$XLVLDLC+cgps$LVLDLC+cgps$MVLDLC+cgps$SVLDLC+cgps$XSVLDLC)

cgps$totalC<-cgps$LDLC+cgps$VLDLC


#Proportion of VLDL in apoB
cgps$VLDL_prop<-cgps$VLDLP/cgps$totalP

cgps$VLDL_prop_cut4<-cut2(cgps$VLDL_prop,g=4)

cgps$LDL_prop<-cgps$LDLP/cgps$totalP
cgps$LDL_prop_cut4<-cut2(cgps$LDL_prop,g=4)

cgps$VLDLC_prop<-cgps$VLDLC/(cgps$VLDLC+cgps$IDLC+cgps$LDLC)

cgps$LDLC_total<-cgps$LDLC+cgps$IDLC

cgps$LDLTG_total<-cgps$LDLTG+cgps$IDLTG



#cOncerting statin treatment variable
cgps$m5<-as.numeric(cgps$m5)


#Correcting format, removing imputed data type
cgps$remncholcalc<-as.numeric(cgps$remncholcalc)
cgps$ldl<-as.numeric(cgps$ldl)
cgps$systolic<-as.numeric(cgps$systolic)
cgps$glucose<-as.numeric(cgps$glucose)
cgps$ppressure<-as.numeric(cgps$systolic)-as.numeric(cgps$diastolic)
cgps$cum_smoking<-as.numeric(cgps$cum_smoking)
cgps$m2a<-as.numeric(cgps$m2a)
cgps$dmall_pre<-as.numeric(cgps$dmall_pre)

cgps$followupPAD<-cgps$stopagePAD-cgps$startage

cgps$birthdate10<-cut2(cgps$birthdate, g=10)

cgps$remncholcalc_log<-log2(cgps$remncholcalc)
cgps$ldl_log<-log2(cgps$ldl)

#CLTI competing tisk variable
cgps$CLTI_comprisk<-as.integer(cgps$CLTI)
cgps$CLTI_comprisk[cgps$dod==TRUE]<-1
cgps$CLTI_comprisk[cgps$st131218==70]<-1
cgps$CLTI_comprisk[cgps$st131218==80]<-1
cgps$CLTI_comprisk[cgps$CLTI==TRUE]<-2
cgps$CLTI_comprisk<-factor(cgps$CLTI_comprisk, labels = c("Censored", "Death", "CLTI"))


#AMI competing tisk variable
cgps$AMI_comprisk<-as.integer(cgps$ami)
cgps$AMI_comprisk[cgps$dod==TRUE]<-1
cgps$AMI_comprisk[cgps$st131218==70]<-1
cgps$AMI_comprisk[cgps$st131218==80]<-1
cgps$AMI_comprisk[cgps$ami==TRUE]<-2
cgps$AMI_comprisk<-factor(cgps$AMI_comprisk, labels = c("Censored", "Death", "AMI"))


#IS competing tisk variable
cgps$IS_comprisk<-as.integer(cgps$is_uval)
cgps$IS_comprisk[cgps$dod==TRUE]<-1
cgps$IS_comprisk[cgps$st131218==70]<-1
cgps$IS_comprisk[cgps$st131218==80]<-1
cgps$IS_comprisk[cgps$is_uval==TRUE]<-2
cgps$IS_comprisk<-factor(cgps$IS_comprisk, labels = c("Censored", "Death", "IS"))


#ASCVD competing tisk variable
cgps$ASCVD_comprisk<-as.integer(cgps$is_uval)
cgps$ASCVD_comprisk[cgps$dod==TRUE]<-1
cgps$ASCVD_comprisk[cgps$st131218==70]<-1
cgps$ASCVD_comprisk[cgps$st131218==80]<-1
cgps$ASCVD_comprisk[cgps$is_uval==TRUE]<-2
cgps$ASCVD_comprisk[cgps$ami==TRUE]<-2
cgps$ASCVD_comprisk[cgps$PAD==TRUE]<-2
cgps$ASCVD_comprisk<-factor(cgps$ASCVD_comprisk, labels = c("Censored", "Death", "ASCVD"))


#ASCVD date
cgps$ASCVD_dato<-pmin(cgps$PAD_dato, cgps$amidato,cgps$isdato_uval)

cgps$stopageAMI<-as.numeric(cgps$amidato-cgps$fodedato)/365.25
cgps$stopageASCVD<-as.numeric(cgps$ASCVD_dato-cgps$fodedato)/365.25

#New names
cgps$IS<-cgps$is_uval
cgps$AMI<-cgps$ami

cgps$ASCVD<-as.logical(cgps$is_uval)
cgps$ASCVD[cgps$ami==TRUE]<-T
cgps$ASCVD[cgps$PAD==TRUE]<-T
cgps$IHD<-cgps$ihd

#CLTI

cgps$stopageCLTI<-as.numeric(cgps$CLTI_dato-cgps$fodedato)/365.25


cgps$followupCLTI<-cgps$stopageCLTI-cgps$startage
cgps$followupIS<-cgps$stopageIS-cgps$startage
cgps$followupAMI<-cgps$stopageAMI-cgps$startage
cgps$followupIHD<-cgps$stopageIHD-cgps$startage
cgps$followupASCVD<-cgps$stopageASCVD-cgps$startage


#Including PAD_pre in ascvd_pre
cgps$ascvd_pre<- with(cgps,ifelse(followupASCVD<0, TRUE, FALSE))

cgps$ascvd_pre[cgps$ihd_pre==1]<-T

# Preparing mediation analysis
cgps$VLDLP_log<-log2(cgps$VLDLP)
cgps$VLDLC_log<-log2(cgps$VLDLC)

cgps$VLDLTG_log<-log2(cgps$VLDLTG)
cgps$trig_log<-log2(cgps$trig)
cgps$remncholcalc_log<-log2(cgps$remncholcalc)

cgps$hscrp_log<-log2(cgps$hscrp)
cgps$Gp_log<-log2(cgps$Gp)

cgps$apob_mgdl_log<-log2(cgps$apob_mgdl)
cgps$LDLC_sqrt<-sqrt(cgps$LDLC)
cgps$fibrinog_log<-log2(cgps$fibrinog)


cgps$cum_smoking_log<-log2(cgps$cum_smoking)
cgps$systolic_log<-log2(cgps$systolic)
cgps$diastolic_log<-log2(cgps$diastolic)

cgps$cum_smoking_log[cgps$cum_smoking_log<0]<-(-2)


#Preparation for mediation analysis

cgps$VLDLP_log<-log2(cgps$VLDLP)

cgps$VLDLC_log<-log2(cgps$VLDLC)

cgps$VLDLTG_log<-log2(cgps$VLDLTG)

cgps$remncholcalc_log<-log2(cgps$remncholcalc)

cgps$apob_cut2<-cut2(cgps$apob_mgdl,100)



# Converting data- obush2 for RDR -----------------------------------------

#For regression dilution bias
cgps2_metabolomics<-read.csv("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstrom/obush2_reg_dil.csv")


for(q in c("","OBUSH1_")){
cgps2_metabolomics[,paste0(q,"LDLP")]<-(cgps2_metabolomics[paste0(q,"SLDLP")]+cgps2_metabolomics[paste0(q,"MLDLP")]+cgps2_metabolomics[paste0(q,"LLDLP")]+cgps2_metabolomics[paste0(q,"IDLP")])
cgps2_metabolomics[,paste0(q,"VLDLP")]<-(cgps2_metabolomics[paste0(q,"XXLVLDLP")]+cgps2_metabolomics[paste0(q,"XLVLDLP")]+cgps2_metabolomics[paste0(q,"LVLDLP")]+cgps2_metabolomics[paste0(q,"MVLDLP")]+cgps2_metabolomics[paste0(q,"SVLDLP")]+cgps2_metabolomics[paste0(q,"XSVLDLP")])
}


cgps2_metabolomics<- merge(cgps2_metabolomics,cgps[c("obushnr","m5","remncholcalc","ldl","sex","startage","birthdate",
                                                     "systolic","diastolic","status_smoking","cum_smoking")],by="obushnr", all=F)


lambda.VLDLP<-lm(log(VLDLP)~log(OBUSH1_VLDLP),data=cgps2_metabolomics,subset=m5==0 & statin2==0)$coefficients[2]
lambda.LDLP<-lm(LDLP~OBUSH1_LDLP,data=cgps2_metabolomics,subset=m5==0 & statin2==0)$coefficients[2]


lambda.apob_mgdl<-0.74
lambda.remncholcalc<-0.48
lambda.ldl<-0.60
lambda.hscrp_log<-0.58 # From article "C-reactive protein concentration and risk of coronary heart disease, stroke, and mortality: an individual participant meta-analysis."


# For lipoprotein subfractions
reg_dil_ratio_all<-read_dta("L:/LovbeskyttetMapper/Forskning KBA/Benjamin Nilsson Wadstrom/reg_dil_ratio_all.dta")

rm(cgps2_metabolomics)

# setting data distribution  ------------------------------------------------------------


#Describing variable distributions to rms package
dd<- datadist(cgps)
dd$limits$sex[2] <- "Female"
dd$limits$current_smoking[2] <- FALSE
dd$limits$m5[2] <- 0
dd$limits$dmall_pre[2] <- 1
dd$limits$ascvd_pre[2] <- FALSE
dd$limits$remncholcalc_cut05[2]<-"<0.5"

dd$limits$remncholcalc[2]<- quantile(cgps$remncholcalc,0.001,na.rm=T)
dd$limits$ldl[2]<-quantile(cgps$ldl,0.001,na.rm=T)
dd$limits$VLDLP[2]<- quantile(cgps$VLDLP,0.001,na.rm=T)
dd$limits$LDLP[2]<- quantile(cgps$LDLP,0.001,na.rm=T)
dd$limits$apob_mgdl[2]<- quantile(cgps$apob_mgdl,0.001,na.rm=T)
dd$limits$VLDL_prop[2]<- quantile(cgps$VLDL_prop,0.001,na.rm=T)

dd$limits$VLDLC[2]<- quantile(cgps$VLDLC,0.001,na.rm=T)
dd$limits$trig[2]<- quantile(cgps$trig,0.001,na.rm=T)
dd$limits$hscrp[2]<- quantile(cgps$hscrp,0.001,na.rm=T)



dd$limits$remncholcalc_log[2]<- log(0.1)
dd$limits$ldl_log[2]<- log(0.5)

options(datadist = "dd")
options(contrasts=c("contr.treatment", "contr.treatment"))

options(na.action = na.omit)



# Descriptive tables ------------------------------------------------------

#Preparation
library(tangram)

EHJ<-nejm
EHJ[['Cell']][["n"]]<-function (n, class = NULL, hdr = FALSE, possible = NULL, ...) 
{
  n <- formatC(n, format="d", big.mark=",")
  cell(n, class = c("cell_n", class), ...)
}
EHJ[['Cell']][['fraction']] <- function(numerator, denominator, format=3, ...)
{ paste0(formatC(numerator, format="d", big.mark=",")," (",render_f(round(100*numerator/denominator, 0), 0),'%)') }
EHJ[['Cell']][['iqr']] <-function (x, format = NA, na.rm = TRUE, names = FALSE, type = 8, 
                                   msd = FALSE, quant = c(0.25, 0.5, 0.75), ...) {
  if (length(quant)%%2 == 0) 
    stop("hmisc_iqr quant argument must be an odd length")
  m <- median(1:length(quant))
  y <- quantile(x, quant, na.rm=TRUE, names, type)
  
  if (median(x, na.rm=T)<10) {
    y<-round(y,1)
  } else {
    y<-round(y,0)
  }
  
  if (is.na(format)) 
    format <- format_guess(y)
  ql <- "-"
  if (sum(!is.na(y)) > 0) {
    ql <- sapply(y, function(x) as.character(x))
    ql <- paste0(ql[m], " (", paste0(ql[1:(m - 1)], 
                                     collapse = "-"), "-", paste0(ql[(m + 
                                                                        1):length(quant)], collapse = "-"), ")")
    if (msd) 
      ql <- paste0(ql, "  ", render_f(mean(x, na.rm = TRUE), 
                                      format), "?", render_f(sd(x, na.rm = TRUE), 
                                                             format))
  }
  cell(ql, ...)
}

EHJ[['Numerical']][['Categorical']] <-function (table, row, column, cell_style, pformat = NULL, msd = FALSE, 
                                                quant = c(0.25, 0.5, 0.75), overall = NULL, test = FALSE, 
                                                useNA = "no", ...) {
  overall_label <- if (is.character(overall)) 
    overall
  else "Overall"
  overall <- column$value != "1" && (isTRUE(overall) || 
                                       is.character(overall))
  datar <- row$data
  datac <- as.categorical(column$data)
  categories <- if (overall) 
    c(levels(datac), overall_label)
  else levels(datac)
  categories <- if (length(categories) == 1) 
    overall_label
  else categories
  format <- ifelse(is.na(row$format), format_guess(datar), 
                   row$format)
  useNA <- useNA == "always" || (sum(is.na(datar)) > 
                                   0 && useNA == "ifany")
  subN <- lapply(levels(datac), FUN = function(cat) {
    cell_style[["n"]](length(datac[datac == cat & !is.na(datac)]), 
                      subcol = cat, hdr = TRUE, possible = length(datac), 
                      ...)
  })
  if (overall) 
    subN[[length(subN) + 1]] <- cell_style[["n"]](sum(!is.na(column$data)), 
                                                  hdr = TRUE, subcol = "Overall", possible = length(column$data), 
                                                  ...)
  stat <- if (length(categories) == 1) {
    ""
  }
  else {
    tst <- suppressWarnings(spearman2(datac, datar, na.action = na.retain))
    cell_style[["fstat"]](f = render_f(tst["F"], 
                                       "%.2f"), df1 = tst["df1"], df2 = tst["df2"], 
                          p = cell_style[["p"]](tst["P"], pformat))
  }
  tbl <- table %>% row_header(derive_label(row))
  if (msd) 
    tbl <- tbl %>% row_header("   Mean?SD")
  if (useNA) 
    tbl <- tbl %>% row_header("   Missing (%)")
  tbl <- if (test) {
    col_header(tbl,  categories, "Test Statistic") %>% 
      col_header("", subN, "")
  }
  else {
    col_header(tbl,  categories) %>% col_header("", 
                                                subN)
  }
  
  tbl <- table_apply(tbl, categories, function(tbl, category) {
    x <- if (category == overall_label) 
      datar
    else datar[datac == category]
    tbl <- tbl %>% add_row(cell_style[["iqr"]](x, 
                                               row$format, subcol = category, na.rm = TRUE, msd = FALSE, 
                                               quant = quant))
    if (msd) {
      tbl <- add_row(tbl, cell(paste0(render_f(mean(x, 
                                                    na.rm = TRUE), row$format), "?", render_f(sd(x, 
                                                                                                 na.rm = TRUE), row$format))))
    }
    if (useNA) {
      tbl <- add_row(tbl, cell_style[["fraction"]](sum(is.na(x)), 
                                                   length(x), format = row$format, subcol = category))
    }
    tbl %>% new_col()
  })
  tbl <- home(tbl) %>% cursor_right(length(categories) + 1)
  if (test) 
    tbl <- add_col(tbl, stat)
  tbl
}

EHJ[['Categorical']][['Categorical']]<-function (table, row, column, cell_style, pformat = NULL, collapse_single = TRUE, 
                                                 overall = NULL, test = FALSE, row_percents = FALSE, useNA = "no", ...) 
{
  grid <- table(as.categorical(row$data), as.categorical(column$data), 
                useNA = useNA)
  if (is.na(colnames(grid)[ncol(grid)])) 
    grid <- grid[, 1:(ncol(grid) - 1)]
  validcol <- which(!apply(grid, 2, FUN = function(x) {
    all(x == 0)
  }))
  validrow <- which(!apply(grid, 1, FUN = function(x) {
    all(x == 0)
  }))
  stat <- if (length(validrow) < 2 || length(validcol) < 2) 
    NA
  else suppressWarnings(chisq.test(grid[validrow, validcol], 
                                   correct = FALSE))
  ncol <- dim(grid)[2]
  nrow <- dim(grid)[1]
  denominators <- if (row_percents) 
    matrix(rep(rowSums(grid), ncol), ncol = ncol, byrow = FALSE)
  else matrix(rep(colSums(grid), nrow), ncol = ncol, byrow = TRUE)
  rowlabels <- rownames(grid)
  subN <- lapply(colnames(grid), FUN = function(cat) cell_style[["n"]](sum(column$data == 
                                                                             cat, na.rm = TRUE), subcol = cat, possible = length(column$data), 
                                                                       hdr = TRUE, ...))
  if (!is.null(overall) && (!is.logical(overall) || overall)) {
    denominators <- cbind(denominators, rep(sum(grid), nrow))
    grid <- cbind(grid, rowSums(grid))
    colnames(grid)[ncol + 1] <- if (is.character(overall)) 
      overall
    else "Overall"
    subN[[ncol + 1]] <- cell_style[["n"]](sum(!is.na(column$data)), 
                                          possible = length(column$data), subcol = "Overall", 
                                          hdr = TRUE, ...)
    ncol <- ncol + 1
  }
  if (collapse_single && dim(grid)[1] <= 2) {
    name <- row$name()
    try({
      l2 <- attr(row$data, "label")
      if (!is.null(l2)) {
        name <- l2
      }
    })
    pos <- dim(grid)[1]
    x <- matrix(grid[pos, ], nrow = 1)
    colnames(x) <- colnames(grid)
    rownames(x) <- paste(name, ":", rownames(grid)[pos])
    grid <- x
    denominators <- matrix(denominators[pos, ], nrow = 1)
    nrow <- 1
  }
  else {
    if (is.na(rownames(grid)[nrow(grid)])) {
      rownames(grid)[nrow(grid)] <- "Missing (%)"
    }
    rownames(grid) <- lapply(rownames(grid), FUN = function(x) paste("  ", 
                                                                     x))
  }
  if (inherits(test, "function")) {
    test_result <- test(row, column, cell_style, ...)
    test <- TRUE
  }
  else if (test) {
    test_result <- if (any(is.na(stat))) 
      cell("NA")
    else cell_style[["chi2"]](render_f(stat$statistic, 
                                       2), stat$parameter, cell_style[["p"]](stat$p.value, 
                                                                             pformat))
  }
  if (test) {
    table <- col_header(table,  colnames(grid), 
                        "Test Statistic")
    table <- col_header(table, subN, "")
  }
  else {
    table <- col_header(table,  colnames(grid))
    table <- col_header(table,  subN)
  }
  if (nrow > 1) 
    table <- row_header(table, derive_label(row, ...))
  for (nm in rownames(grid)) table <- row_header(table, nm)
  for (j in 1:ncol) {
    if (nrow > 1) 
      table <- add_row(table, "")
    format <- if (is.na(row$format) || is.null(row$format)) 
      format_guess(as.vector(grid/denominators))
    else row$format
    for (i in 1:nrow) {
      table <- if (denominators[i, j] == 0) 
        add_row(table, "")
      else add_row(table, cell_style[["fraction"]](grid[i, 
                                                        j], denominators[i, j], format = format, subcol = colnames(grid)[i], 
                                                   subrow = rownames(grid)[j]))
    }
    table <- new_col(table)
  }
  if (test) {
    table <- add_row(table, test_result)
    if (nrow > 1) 
      table <- add_row(table, rep("", nrow))
  }
  table
}



cgps$ldlmgdl<-cgps$ldl*38.67
cgps$remncholcalcmgdl<-cgps$remncholcalc*38.67
cgps$trigmgdl<-cgps$trig*88.57
cgps$hdlmgdl<-cgps$hdl*38.67



# All individuals

all <- tangram("1~ sex + startage+apob_mgdl+ remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+ VLDLP+LDLP+hdl+hdlmgdl+ systolic+diastolic+factor(current_smoking)+ cum_smoking+ bmi+factor(dmall_pre)+ factor(ascvd_pre)+trig+trigmgdl+hdl+hdlmgdl",
               data = cgps, 
               test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)


csv(all, file="Temp.csv")


all <- tangram("1~ sex + startage+apob_mgdl+ remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+ VLDLP+LDLP+hdl+hdlmgdl+ systolic+diastolic+factor(current_smoking)+ cum_smoking+ bmi+factor(dmall_pre)+ factor(ascvd_pre)+trig+trigmgdl+hdl+hdlmgd",
               data = subset(cgps,!is.na(VLDLP) ), 
               test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)


csv(all, file="Temp.csv")


#Table S4
endpoint_descb_all <- tangram("1 ~c_diag",
                              data = subset(cgps,PAD==T & followupPAD>0), 
                              test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)

csv(endpoint_descb_all, file="Endpoint characteristics.csv")

#Table S5
endpoint_descb_all <- tangram("1 ~CLTI_op",
                              data = subset(cgps,CLTI==T & followupCLTI>0), 
                              test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)

csv(endpoint_descb_all, file="Endpoint characteristics.csv")



PAD_desc <- tangram("PAD ~ sex + startage+apob_mgdl+ remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+ VLDLP+LDLP+hdl+hdlmgdl+ systolic+diastolic+factor(current_smoking)+ cum_smoking+ bmi+factor(dmall_pre)+ factor(ascvd_pre)+trig+trigmgdl",
                    data = subset(cgps,followupPAD>0),
                    test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)

PAD_desc <- tangram("1 ~ sex + startage+apob_mgdl+ remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+ VLDLP+LDLP+hdl+hdlmgdl+ systolic+diastolic+factor(current_smoking)+ cum_smoking+ bmi+factor(dmall_pre)+ factor(ascvd_pre)+trig+trigmgdl",
                    data = subset(cgps,followupPAD>0),
                    test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)
csv(PAD_desc, file="Temp.csv")


CLTI_desc <- tangram("CLTI ~ sex + startage+apob_mgdl+ remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+ VLDLP+LDLP+hdl+hdlmgdl+ systolic+diastolic+factor(current_smoking)+ cum_smoking+ bmi+factor(dmall_pre)+ factor(ascvd_pre)+trig+trigmgdl",
                    data = subset(cgps,followupCLTI>0),
                    test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)

CLTI_desc <- tangram("1 ~ sex + startage+apob_mgdl+ remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+ VLDLP+LDLP+hdl+hdlmgdl+ systolic+diastolic+factor(current_smoking)+ cum_smoking+ bmi+factor(dmall_pre)+ factor(ascvd_pre)+trig+trigmgdl",
                    data = subset(cgps,followupCLTI>0),
                    test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)
csv(CLTI_desc, file="Temp.csv")


AMI_desc <- tangram("AMI ~ sex + startage+apob_mgdl+ remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+ VLDLP+LDLP+hdl+hdlmgdl+ systolic+diastolic+factor(current_smoking)+ cum_smoking+ bmi+factor(dmall_pre)+ factor(ascvd_pre)+trig+trigmgdl",
                    data = subset(cgps,followupAMI>0),
                    test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)

AMI_desc <- tangram("1 ~ sex + startage+apob_mgdl+ remncholcalc+ remncholcalcmgdl+ ldl+ ldlmgdl+ VLDLP+LDLP+hdl+hdlmgdl+ systolic+diastolic+factor(current_smoking)+ cum_smoking+ bmi+factor(dmall_pre)+ factor(ascvd_pre)+trig+trigmgdl",
                    data = subset(cgps,followupAMI>0),
                    test=FALSE, id="override", transforms = EHJ, msd=FALSE, digits = 1)


csv(AMI_desc, file="Temp.csv")



# Figure 1: bar plot particles fractions -------




subfractions<-c("SLDLP","MLDLP","LLDLP","IDLP","XSVLDLP","SVLDLP","MVLDLP","LVLDLP","XLVLDLP","XXLVLDLP")

ldl_remnant<-c("ldl","remncholcalc")

apob_mgdl<-"apob_mgdl"


fig2<-list()

for(k in c("apob_mgdl","ldl_remnant","subfractions")){
  df<-data.frame()
  for (p in rev(get(k))){
    
    dat<-subset(cgps,!is.na(XXLVLDLP))
    df[paste(k,p),c("median","particle","x")]<-c(median(dat[,p],na.rm=T),p,"All")
    
  }
  df$particle<-factor(df$particle,levels=df$particle,labels=df$particle)
  df$median<-as.numeric(df$median)
  fig2[k]<- list(ggplot(data=df, aes(x=x, y=median, fill=particle)) + 
                   geom_bar(position="stack", stat="identity",color="black")+scale_fill_manual(values=  colorRampPalette(colors = c("darkcyan","white", "darkorchid4"))(length(get(k))))+
                   labs(y="Number of particles, nmol/L")+
                   scale_y_continuous(expand = c(0,0)) + 
                   theme(panel.background =element_blank(),plot.title = element_text(size=12),axis.line=element_line(size=1),axis.ticks.length = unit(.2, "cm"),axis.ticks.y=element_line(size=1)))
  
  
}

#Changing function plot_grid
fargs  <- formals(plot_grid)
fargs$ncol           <- 3
formals(plot_grid) <- fargs


figure2<-do.call(plot_grid,fig2)
figure2



toffice(figure2,format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE, left = 0.15,top = 0.15,
        append = FALSE, width = 10, height = 6,devsize = FALSE, units = "in")



# Figure 2: Correlation matrix and scatterplots particle number and cholesterol content -------

#Correlation matrix

library(reshape2)
corr_list<-list()

corrmatrix_P<-round(cor(cgps[c("SLDLP","MLDLP","LLDLP","IDLP","XSVLDLP","SVLDLP","MVLDLP","LVLDLP","XLVLDLP","XXLVLDLP","ldl","remncholcalc","apob_mgdl")], use = "complete",method="spearman")^2,2)


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)}

corrmatrix_P<-get_upper_tri(corrmatrix_P)

corrmatrix_P <- melt(corrmatrix_P, na.rm = TRUE)

Figure_1<-ggplot(corrmatrix_P, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low ="white", high = "red",  
                       limit = c(0,1), space = "Lab", 
                       name="Spearman\nCorrelation^2") +
  theme_nothing()+ 
  coord_fixed() + 
  geom_text(aes(Var2, Var1, label = str_remove(value, "^0+")), color = "black", size = 4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.text.y = element_text(size=12),
        legend.justification = c(1, 0),
        legend.position = c(0.6, 0.7),
        legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))



toffice(figure = Figure_1, format = "pptx",title= "Copenhagen General Population Study",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = TRUE, width = 8, height = 8, devsize = FALSE, units = "in")




# Figure 4: RCS plots -----------------------------------------------


fig3<-list()
for (k in c("PAD","CLTI","AMI"))
  for (i in c("apob_mgdl","remncholcalc","ldl","VLDLP","LDLP")){
    model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~rcs(",i,",3)+strat(sex)+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate")) ,data=cgps,x=TRUE, y=TRUE)  
    pred<-Predict(model,name=i, ref.zero = TRUE, fun=exp)
    
    pred[c("yhat","lower","upper")]<-exp(log(pred[c("yhat","lower","upper")])/get(paste0("lambda.",i)))
    
    
    col<-if (i=="remncholcalc"|i=="VLDLP"){"darkcyan"} else if (i=="apob_mgdl") {"black"} else {"darkorchid4"}
    
    gg<-ggplot(pred,colfill = col)+
                              geom_line(size=1, colour=col) +
                              geom_hline(aes(yintercept = 1), linetype = 2)+ coord_cartesian( xlim=c(0,signif(quantile(cgps[i],0.995,na.rm=T),1)),ylim= c(0.5, 10))+
                              labs(x=i , y=paste("Hazard ratio for",k,"(95% CI)"))+scale_y_continuous(trans = "log",breaks = c(0.5,1,2,4,6,8,10))+
                              theme(axis.line.x.bottom =element_line(),axis.line.y.left =element_line() ,axis.line.x.top=element_blank() ,plot.margin=margin(0,3,0,3),plot.subtitle = element_blank(), plot.title = element_blank(),axis.text.y = element_text(size=10), 
                                    panel.grid=element_blank(),strip.background = element_rect(color = "white", fill = "white"), panel.background = element_rect(fill = "white"))
    
    gg_dens<-ggplot(data=cgps) + 
      geom_density(aes_string(x=i),color=col, fill=col, alpha=0.4)+
      theme_nothing()+scale_x_continuous(limits=c(0,signif(quantile(cgps[i],0.995,na.rm=T),1)))
    
    
    fig3[paste0(i,k)]<- list(plot_grid(if(k=="PAD"){gg_dens} else {ggplot()+theme_nothing()},
                                       gg,ncol=1,nrow=2,
                                       rel_heights = if(k=="PAD"){c(0.3,1)}else {c(0.1,1)}))
  }


#Changing function plot_grid
fargs  <- formals(plot_grid)
fargs$ncol           <- 5
formals(plot_grid) <- fargs

figure3<-do.call(plot_grid,fig3)


toffice(figure = figure3, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = TRUE, width = 11, height = 7, devsize = FALSE, units = "in")

# Figure 5: By mmol-increase in remnant and LDL ------------------------------------------------------


fig4<-data.frame()
for (k in c("PAD","CLTI","AMI"))
  for (i in c( "apob_mgdl","remncholcalc","ldl","VLDLP","LDLP")){
    
    fig4[paste(i,k),"var"]<-paste0(i,k)
    p_years<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~1")) ,data=subset(cgps,!is.na(get(i))),scale=1)
    
    fig4[paste(i,k),"N"]<-p_years$n
    fig4[paste(i,k),"N_events"]<-p_years$event
    fig4[paste(i,k),"events_pyears"]<-round(p_years$event/p_years$pyear*1000,1)
    
    
    for (g in c("","+hdl")){
    model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,"+",
                                 if(i=="remncholcalc") {"ldl"} else if(i=="ldl") {"remncholcalc"} else if(i=="VLDLP") {"LDLP"} else if(i=="LDLP") {"VLDLP"} ,
                                 "+strat(sex)+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate",g)) ,data=cgps,x=TRUE, y=TRUE)  
    
    
    pred<- if (i=="remncholcalc") {
      Predict(model,remncholcalc=c(dd$limits$remncholcalc[2],dd$limits$remncholcalc[2]+1), ref.zero = TRUE, fun=exp)
    } else if (i=="ldl"){
      Predict(model, ldl=c(dd$limits$ldl[2],dd$limits$ldl[2]+1),ref.zero = TRUE, fun=exp)
    } else if (i=="VLDLP") {
      Predict(model, VLDLP=c(dd$limits$VLDLP[2],dd$limits$VLDLP[2]+100),ref.zero = TRUE, fun=exp)
    } else if (i=="LDLP"){
      Predict(model, LDLP=c(dd$limits$LDLP[2],dd$limits$LDLP[2]+100),ref.zero = TRUE, fun=exp)
    } else if (i=="apob_mgdl"){
      Predict(model, apob_mgdl=c(dd$limits$apob_mgdl[2],dd$limits$apob_mgdl[2]+50),ref.zero = TRUE, fun=exp)
    }
    
    pred[c("yhat","lower","upper")]<-exp(log(pred[c("yhat","lower","upper")])/get(paste0("lambda.",i)))
   
  
    fig4[paste(i,k),paste0(c("yhat","lower","upper"),g)]<-pred[2,c("yhat","lower","upper")] 
    fig4[paste(i,k),paste0("pval",g)]<-pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2
    fig4[paste(i,k),paste0("coef",g)]<-model$coef[1]
    fig4[paste(i,k),paste0("se",g)]<-sqrt(model$var[1,1])
    fig4[paste(i,k),paste0("variable",g)]<-i
    }
    }

fig4

fig4$var<-factor(fig4$var,levels=fig4$var)

#z-test

for (i in c( "apob_mgdl","remncholcalc","ldl","VLDLP","LDLP"))
  for (k in c("PAD","CLTI")){

dat<-subset(fig4,variable==i)  
pdiff<-pnorm(abs(dat[grepl(k,dat$var),"coef"]-dat[grepl("AMI",dat$var),"coef"])/
                   sqrt(dat[grepl(k,dat$var),"se"]^2+dat[grepl("AMI",dat$var),"se"]^2),lower.tail = F)*2

print(paste(pdiff,i,k))
}



#Creating table
fig4_table <- tibble(cbind(fig4$var,fig4[grepl("yhat",colnames(fig4))|grepl("pval",colnames(fig4))])) %>% flextable()
fig4_table<- compose(fig4_table,j=paste0("yhat"), value=as_paragraph(round(fig4[,paste0("yhat")], digits = 1), " (", round(fig4[,paste0("lower")], digits = 1),"-", round(fig4[,paste0("upper")], digits = 1), ")"))
fig4_table<- compose(fig4_table,j=paste0("pval"), value=as_paragraph(signif.p(fig4[,paste0("pval")]),as_sup(signif.p(fig4[,paste0("pval")],exp=T))))
fig4_table<- compose(fig4_table,j=paste0("yhat+hdl"), value=as_paragraph(round(fig4[,paste0("yhat+hdl")], digits = 1), " (", round(fig4[,paste0("lower+hdl")], digits = 1),"-", round(fig4[,paste0("upper+hdl")], digits = 1), ")"))
fig4_table<- compose(fig4_table,j=paste0("pval+hdl"), value=as_paragraph(signif.p(fig4[,paste0("pval+hdl")]),as_sup(signif.p(fig4[,paste0("pval+hdl")],exp=T))))

fig4_table<-fontsize(fig4_table, size = 11, part = "body")
fig4_table<-align(autofit(fig4_table, part="all"),align="left", part = "all")
fig4_table


#figure

group.colors <- rep(c("darkblue","darkcyan","darkorchid4","darkcyan","darkorchid4"),3)

figure4<-ggplot(fig4, aes(y=var,x=yhat, xmin =lower, xmax = upper))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme(axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_point(size=3,colour=group.colors)+
  geom_errorbar(width=0.5, cex=0.7,color=group.colors, position=position_dodge(width=0.1))+ 
  scale_colour_grey(start = 0.75, end = 0)+scale_x_continuous(breaks=seq(1,3,1))+
  scale_y_discrete(limits=rev)  +coord_cartesian( xlim = c(0.75,3))

figure4b<-ggplot(data=fig4, aes(y=var,x=`yhat+hdl`, xmin =`lower+hdl`, xmax = `upper+hdl`))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme(axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_point(size=3,colour=group.colors)+
  geom_errorbar(width=0.5, cex=0.7,color=group.colors, position=position_dodge(width=0.1))+ 
  scale_colour_grey(start = 0.75, end = 0)+scale_x_continuous(breaks=seq(1,3,1))+
  scale_y_discrete(limits=rev)  +coord_cartesian( xlim = c(0.75,3))





toffice(figure = plot_grid(figure4,figure4b,ncol=2), format = "pptx", title = "Copenhagen General Population Study",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = FALSE, width = 2.5, height = 4.5,devsize = FALSE, units = "in")






# Figure 6: mediation analysis apob-----------------------------------


mediation_results<-data.frame()
for (k in c("PAD","CLTI","AMI"))
  for ( i in c("remncholcalc_log","ldl","VLDLP_log","LDLP")){
    
    dat<-cgps[!is.na(cgps[,i]),]
    
    obj<-cmest(data=subset(subset(dat, if (k=="PAD"){followupPAD>0} 
                                  else if (k=="AMI"){followupAMI>0}
                                  else if (k=="IS"){followupIS>0}
                                  else if (k=="ASCVD"){followupASCVD>0}
                                  else {sex!="hej"}), if (i=="hscrp_log"){hscrp>0} else {sex!="hej"}),
               model="rb", estimation = "paramfunc", inference = "delta", outcome=paste0("followup",k), event=k, 
               exposure = "apob_cut2", mediator=c(i), EMint=F, 
               basec= if (i=="VLDLP_log"| i=="LDLP"){c("startage","sex", "status_smoking", "cum_smoking", "birthdate", "systolic","diastolic")} else {
                 c("startage","sex", "status_smoking", "cum_smoking", "birthdate", "systolic","diastolic")},
               yreg="coxph", mreg=list("linear"),mval=list(1) , yval=TRUE, nboot=1, multimp = F)
    
    mediation_results[paste(i,k),c("mediation","lower","upper", "P.val")]<-as.data.frame(summary(obj)[c("effect.pe","effect.ci.low","effect.ci.high","effect.pval")])["pm", ]
    
    
    print(paste(i,k,"done!"))
    
  }

mediation_results$variable<-factor(rownames(mediation_results),levels=rownames(mediation_results))
mediation_results

#Converting P-values of negative mediation
mediation_results[mediation_results[,"mediation"]<0,"P.val"]<-1-mediation_results[mediation_results[,"mediation"]<0,"P.val"]

#Minumum mediation 0, maximum mediation 1
for (i in c("mediation", "lower", "upper")){
  mediation_results[mediation_results[,i]<0,i]<-0
  mediation_results[mediation_results[,i]>1,i]<-1
}
mediation_results


#Creating tables for figures
fig5_table<- tibble(mediation_results$variable,mediation=paste(signif(mediation_results$mediation,2)),P.val=paste(signif(mediation_results$P.val,2))) %>% flextable()
fig5_table<- compose(fig5_table,j="mediation", value=as_paragraph(sub("^.0^","",round(mediation_results$mediation*100))," (", 
                                                                  sub("^.0^","",round(mediation_results$lower*100)),"-", 
                                                                  sub("^.0^","",round(mediation_results$upper*100)),")"))
fig5_table<-align(autofit(fig5_table, part="all"),align="left", part = "all")
fig5_table


group.colors <- rep(c("darkcyan","darkorchid4","darkcyan","darkorchid4"),3)

figure6cgps<-ggplot(data=mediation_results, aes(y = variable ,  x = mediation*100,  xmin=lower*100,xmax=upper*100, fill = variable ))+
  labs(y="Variable", x="Explained risk, %")+
  theme_nothing()+geom_col(width = 0.6, fill=group.colors, size=1)+geom_vline(xintercept=0, size=0.5, linetype=1)+
  theme(axis.ticks.length.x = unit(0.2, "cm"),axis.title = element_text(size=10), axis.text.x = element_text(vjust=0.1),axis.text.y = element_text(hjust=0.4, vjust=0.3), axis.title.y = element_text(vjust = 3, angle=90,size=10), axis.ticks=element_line(size=1), axis.text = element_text(size=10), axis.line.x.bottom =element_line(size=1),axis.line.y.left =element_blank() ,axis.line.x.top=element_blank() ,plot.margin=margin(0,3,0,3),plot.subtitle = element_blank(), plot.title = element_text(), 
        panel.grid=element_blank())+scale_x_continuous(expand=c(0,0))+scale_y_discrete(limits=rev)+coord_cartesian(xlim=c(0,100))



toffice(figure = figure6cgps, format = "pptx", title = "Copenhagen General Population Study",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = FALSE, width = 4, height = 3,devsize = FALSE, units = "in")




# Figure SX: Venn diagram of endpoints ------------------------------------


library(ggvenn)


x <- list(
  "Peripheral artery disease" = subset(cgps,stopagePAD>startage &PAD==T)$obushnr, 
  "Chronic limb-threatening ischemia" = subset(cgps,stopageCLTI>startage &CLTI==T)$obushnr, 
  "Myocardial infarction" = subset(cgps,stopageAMI>startage &AMI==T)$obushnr
)


p<-ggvenn (x,  stroke_size = 0.5, set_name_size = 5)



toffice(figure = p, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = TRUE, width = 7, height = 7, devsize = FALSE, units = "in")



# Figure SX: LOESS particle concentration vs apoB --------------------------


fig2<-list()
for (i in c("remncholcalc","VLDLP","ldl","LDLP")){
  
  col<-if (i=="remncholcalc"|i=="VLDLP"){"darkcyan"} else if (i=="apob_mgdl") {"black"} else {"darkorchid4"}
  
  
  gg<-ggplot(data=cgps, aes(x=apob_mgdl))+theme_nothing()+
    theme(axis.line = element_line(), axis.title.x = element_text(),axis.title.y = element_text(angle = 90),axis.ticks = element_line(),axis.text = element_text(),axis.ticks.length = unit(.2, "cm"))+
    #geom_point(aes_string(y=i),color="grey30")+
    geom_smooth(aes_string(y=i),color=col,fill=col)+
    xlab("Apolipoprotein B, mg/dL")+ylab(i)+
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(ylim = c(0,quantile(cgps[,i],0.99, na.rm=T)),xlim=c(0,200))
  
  gg_dens<-ggplot(data=cgps) + 
    geom_density(aes_string(x=i),color=col, fill=col, alpha=0.4)+
    theme_nothing()+ coord_flip()+ scale_x_continuous(limits=c(0,quantile(cgps[,i],0.99, na.rm=T)))
  
  
  fig2[i]<- list(plot_grid(gg,gg_dens,ncol=2,rel_widths = c(1,0.3)))
  
}

#Changing function plot_grid
fargs  <- formals(plot_grid)
fargs$ncol           <- 2
formals(plot_grid) <- fargs

figure2<-do.call(plot_grid,fig2)

gg_dens_apob<-ggplot(data=cgps) + 
  geom_density(aes_string(x="apob_mgdl"),color="black", fill="black", alpha=0.4)+
  theme_nothing()+ scale_x_continuous(limits=c(0,200))

figure2<-plot_grid(plot_grid(gg_dens_apob,gg_dens_apob,ncol=2),figure2,nrow=2, rel_heights = c(0.15,1),ncol=NULL)


toffice(figure = figure2, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = TRUE, width = 7, height = 7, devsize = FALSE, units = "in")


# Figure SX: per IQR incresae ------------------------------------------------------


fig4<-data.frame()
for (k in c("PAD","CLTI","AMI"))
  for (i in c( "apob_mgdl","remncholcalc","ldl","VLDLP","LDLP")){
    
    p_years<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~1")) ,data=subset(cgps,!is.na(get(i))),scale=1)
    
    fig4[paste(i,k),"N"]<-p_years$n
    fig4[paste(i,k),"N_events"]<-p_years$event
    fig4[paste(i,k),"events_pyears"]<-p_years$event/p_years$pyear*1000
    
    
    for (g in c("","+hdl")){
      model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,"+",
                                   if(i=="remncholcalc") {"ldl"} else if(i=="ldl") {"remncholcalc"} else if(i=="VLDLP") {"LDLP"} else if(i=="LDLP") {"VLDLP"} ,
                                   "+strat(sex)+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate",g)) ,data=cgps,x=TRUE, y=TRUE)  
    
    
    pred<- if (i=="remncholcalc") {
      as.data.frame(summary(model,remncholcalc,est.all=F,antilog=F))
    } else if (i=="ldl"){
      as.data.frame(summary(model,ldl,est.all=F,antilog=F))
    } else if (i=="VLDLP") {
      as.data.frame(summary(model,VLDLP,est.all=F,antilog=F))
    } else if (i=="LDLP"){
      as.data.frame(summary(model,LDLP,est.all=F,antilog=F))
    } else if (i=="apob_mgdl"){
      as.data.frame(summary(model,apob_mgdl,est.all=F,antilog=F))
    }
    
    pred[c("yhat","lower","upper")]<-exp(pred[i,c("Effect","Lower 0.95", "Upper 0.95")]/get(paste0("lambda.",i)))
    
    
    
    fig4[paste(i,k),paste0(c("yhat","lower","upper"),g)]<-pred[1,c("yhat","lower","upper")] 
    fig4[paste(i,k),paste0("pval",g)]<-pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2
    fig4[paste(i,k),paste0("coef",g)]<-model$coef[1]
    fig4[paste(i,k),paste0("se",g)]<-sqrt(model$var[1,1])
    fig4[paste(i,k),paste0("variable",g)]<-i
    }

  }

fig4

fig4$var<-factor(rownames(fig4),levels=rownames(fig4))

#Formatting P-values for table with blank spaces
#Table
fig4_table <- tibble(cbind(fig4$var,fig4[grepl("yhat",colnames(fig4))|grepl("pval",colnames(fig4))])) %>% flextable()
fig4_table<- compose(fig4_table,j=paste0("yhat"), value=as_paragraph(round(fig4[,paste0("yhat")], digits = 1), " (", round(fig4[,paste0("lower")], digits = 1),"-", round(fig4[,paste0("upper")], digits = 1), ")"))
fig4_table<- compose(fig4_table,j=paste0("pval"), value=as_paragraph(signif.p(fig4[,paste0("pval")]),as_sup(signif.p(fig4[,paste0("pval")],exp=T))))
fig4_table<- compose(fig4_table,j=paste0("yhat+hdl"), value=as_paragraph(round(fig4[,paste0("yhat+hdl")], digits = 1), " (", round(fig4[,paste0("lower+hdl")], digits = 1),"-", round(fig4[,paste0("upper+hdl")], digits = 1), ")"))
fig4_table<- compose(fig4_table,j=paste0("pval+hdl"), value=as_paragraph(signif.p(fig4[,paste0("pval+hdl")]),as_sup(signif.p(fig4[,paste0("pval+hdl")],exp=T))))

fig4_table<-fontsize(fig4_table, size = 11, part = "body")
fig4_table<-align(autofit(fig4_table, part="all"),align="left", part = "all")
fig4_table

#figure


group.colors <- rep(c("darkblue","darkcyan","darkorchid4","darkcyan","darkorchid4"),3)

figure4<-ggplot(fig4, aes(y=var,x=yhat, xmin =lower, xmax = upper))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme(axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_point(size=3,colour=group.colors)+
  geom_errorbar(width=0.5, cex=0.7,color=group.colors, position=position_dodge(width=0.1))+ 
  scale_colour_grey(start = 0.75, end = 0)+scale_x_continuous(breaks=seq(1,2,0.5))+
  scale_y_discrete(limits=rev)  +coord_cartesian( xlim = c(0.5,2))

figure4b<-ggplot(data=fig4, aes(y=var,x=`yhat+hdl`, xmin =`lower+hdl`, xmax = `upper+hdl`))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme(axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_point(size=3,colour=group.colors)+
  geom_errorbar(width=0.5, cex=0.7,color=group.colors, position=position_dodge(width=0.1))+ 
  scale_colour_grey(start = 0.75, end = 0)+scale_x_continuous(breaks=seq(1,2,0.5))+
  scale_y_discrete(limits=rev)  +coord_cartesian( xlim = c(0.5,2))





toffice(figure = plot_grid(figure4,figure4b,ncol=2), format = "pptx", title = "Copenhagen General Population Study",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = FALSE, width = 2.5, height = 4.5,devsize = FALSE, units = "in")






# Figure Sx: Particle subfractions  ----------------------------------



fig4<-data.frame()
for (k in c("PAD","CLTI","AMI"))
  for (i in rev(c("SL","ML","LL","I","XSVL","SVL","MVL","LVL","XLVL","XXLVL"))){
    
    
    
    for (q in c("model1","model2","model3")) {
      
      adjust<-paste0("+strat(sex)+birthdate",
                     if (q=="model2" |q=="model3"){"+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)"} else {""},
                     if (q=="model3"){"+bmi+ascvd_pre+dmall_pre"} else {""})
      
      model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,"DLP",adjust)) ,data=cgps,x=TRUE, y=TRUE)  
      
      pred<-as.data.frame(summary(model))
      
      pred[c("Effect","Lower 0.95","Upper 0.95")]<-exp(pred[i,c("Effect","Lower 0.95", "Upper 0.95")]/subset(reg_dil_ratio_all,parm==paste0(i,"DLP"))$rdb)
      
      fig4[paste0(i,k),paste0(c("yhat","lower","upper"),q)]<-pred[2,c("Effect","Lower 0.95","Upper 0.95")]
      fig4[paste0(i,k),paste0("pval",q)]<-pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2
      
      
      
      
    }
    
    p_years<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~1")) ,data=subset(cgps,!is.na(get(paste0(i,"DLP")))),scale=1)
    
    fig4[paste0(i,k),"N"]<-p_years$n
    fig4[paste0(i,k),"N_events"]<-p_years$event
    fig4[paste0(i,k),"events_pyears"]<-p_years$event/p_years$pyear*1000
    
    print(paste(i,k,"done!"))
  }

fig4

fig4$var<-factor(rownames(fig4),levels=rownames(fig4))


#Table
fig4_table <- tibble(Endpoint=fig4$var,
                     individuals=fig4$N,
                     nevents= fig4$N_events,
                     personyears=round(fig4$events_pyears,1),
                     HR1=fig4[,paste0("yhat","model1")],P1=paste(signif(fig4[,paste0("pval","model1")],2)), HR2=fig4[,paste0("yhat","model2")], P2=paste(signif(fig4[,paste0("pval","model2")],2)),HR3=fig4[,paste0("yhat","model3")],P3=paste(signif(fig4[,paste0("pval","model3")],2))  ) %>% flextable()
fig4_table<- compose(fig4_table,j="HR1", value=as_paragraph(round(fig4[,paste0("yhat","model1")], digits = 1), " (", round(fig4[,paste0("lower","model1")], digits = 1),"-", round(fig4[,paste0("upper","model1")], digits = 1), ")"))
fig4_table<- compose(fig4_table,j="HR2", value=as_paragraph(round(fig4[,paste0("yhat","model2")], digits = 1), " (", round(fig4[,paste0("lower","model2")], digits = 1),"-", round(fig4[,paste0("upper","model2")], digits = 1), ")"))
fig4_table<- compose(fig4_table,j="HR3", value=as_paragraph(round(fig4[,paste0("yhat","model3")], digits = 1), " (", round(fig4[,paste0("lower","model3")], digits = 1),"-", round(fig4[,paste0("upper","model3")], digits = 1), ")"))

align(autofit(fig4_table, part="all"),align="left", part = "all")


#figure

group.colors<-rep(colorRampPalette(colors = c("darkcyan","darkgrey", "darkorchid4"))(10),3)

figure4a<-ggplot(fig4, aes(y=var,x=yhatmodel1, xmin =lowermodel1, xmax = uppermodel1))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme( axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_errorbar(width=0.5, cex=0.5,color="black", position=position_dodge(width=0.1))+ 
  geom_point(size=2,colour=group.colors)+
  scale_colour_grey(start = 0.75, end = 0)+
  scale_y_discrete(limits=rev)+scale_x_continuous(lim=c(0.5,2),breaks = c(0.5,0.75,1,1.5,2))



figure4b<-ggplot(fig4, aes(y=var,x=yhatmodel2, xmin =lowermodel2, xmax = uppermodel2))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme( axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_errorbar(width=0.5, cex=0.5,color="black", position=position_dodge(width=0.1))+ 
  geom_point(size=2,colour=group.colors)+
  scale_colour_grey(start = 0.75, end = 0)+
  scale_y_discrete(limits=rev) +scale_x_continuous(lim=c(0.5,2),breaks = c(0.5,0.75,1,1.5,2))

figure4c<-ggplot(fig4, aes(y=var,x=yhatmodel3, xmin =lowermodel3, xmax = uppermodel3))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme( axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_errorbar(width=0.5, cex=0.5,color="black", position=position_dodge(width=0.1))+ 
  geom_point(size=2,colour=group.colors)+
  scale_colour_grey(start = 0.75, end = 0)+
  scale_y_discrete(limits=rev)+scale_x_continuous(lim=c(0.5,2),breaks = c(0.5,0.75,1,1.5,2))


figure4<-plot_grid(figure4a,figure4b,figure4c,ncol=3)


toffice(figure = figure4, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = TRUE, width = 4, height = 6.5, devsize = FALSE, units = "in")



# Figure Sx: sensitivity analysis different models ----------------------------------


fig4<-data.frame()
for (k in c("PAD","AMI"))
    for (q in c("model1","model2","model3","model4","model5","model6","model7")){
      
dat<-cgps[c("apob_mgdl","remncholcalc","ldl","VLDLP","LDLP","sex","birthdate","obushnr","PAD","AMI","stopagePAD","stopageAMI","PAD_comprisk","AMI_comprisk","followupPAD","followupAMI","startage","systolic","diastolic","cum_smoking","status_smoking","chol", "XSVLDLP","trig","hdl","P1methyl2pyrrolidone","bmi","ascvd_pre","dmall_pre")]
      
        
      if (q=="model4"){
        dat$VLDLP<-dat$VLDLP-dat$XSVLDLP
        dat$LDLP<-dat$LDLP+dat$XSVLDLP
      } else if (q=="model7")  {
        dat$non_hdl_sampson<-(dat$chol-dat$hdl)*38.67
        dat$chol_sampson<-dat$chol*38.67
        dat$trig_sampson<-dat$trig*88.57
        dat$hdl_sampson<-dat$hdl*38.67
        
        dat$ldl_sampson<-dat$chol_sampson/0.948 - dat$hdl_sampson/0.971 - (dat$trig_sampson/8.56 + dat$trig_sampson*dat$non_hdl_sampson/2140 - dat$trig_sampson^2/16100) - 9.44
        dat$ldl<-dat$ldl_sampson/38.67
        
        dat$remncholcalc<-dat$chol-dat$ldl-dat$hdl
        
      } else {}
      
      if (q=="model5"){
        comprisk_dat<-finegray(as.formula(paste0("Surv(startage,stopage",k,",",k,"_comprisk)~.")), etype=k, data=subset(dat,!is.na(get(paste0(i))) & get(paste0("followup",k))>0),id=obushnr, count=1)  
        model<-cph(as.formula(paste0("Surv(fgstart,fgstop,fgstatus)~",i,adjust)) ,data=comprisk_dat,x=TRUE, y=TRUE,weights=fgwt)  
      }else {}
      
      for (i in c("apob_mgdl","remncholcalc","ldl","VLDLP","LDLP")){
      
      
      adjust<-paste0("+", if(i=="remncholcalc") {"ldl"} else if(i=="ldl") {"remncholcalc"} else if(i=="VLDLP") {"LDLP"} else if(i=="LDLP") {"VLDLP"},
                     "+strat(sex)+birthdate+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)",
                     if (q=="model2"){"+bmi+dmall_pre+ascvd_pre"} else {""},
                     if (q=="model6") {"+rcs(startage,3)"} else {})
      
      
      if (q=="model5"){
      model<-cph(as.formula(paste0("Surv(fgstart,fgstop,fgstatus)~",i,adjust)) ,data=comprisk_dat,x=TRUE, y=TRUE,weights=fgwt)  
      
      } else if (q=="model6"){
      model<-cph(as.formula(paste0("Surv(followup",k,",",k,")~",i,adjust)) ,data=subset(dat,get(paste0("followup",k))>0),x=TRUE, y=TRUE) 
      } else {
      model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,adjust)) ,data=dat,x=TRUE, y=TRUE)  
      }
      
      
      pred<- if (i=="remncholcalc") {
        Predict(model,remncholcalc=c(dd$limits$remncholcalc[2],dd$limits$remncholcalc[2]+1), ref.zero = TRUE, fun=exp)
      } else if (i=="ldl"){
        Predict(model, ldl=c(dd$limits$ldl[2],dd$limits$ldl[2]+1),ref.zero = TRUE, fun=exp)
      } else if (i=="VLDLP") {
        Predict(model, VLDLP=c(dd$limits$VLDLP[2],dd$limits$VLDLP[2]+100),ref.zero = TRUE, fun=exp)
      } else if (i=="LDLP"){
        Predict(model, LDLP=c(dd$limits$LDLP[2],dd$limits$LDLP[2]+100),ref.zero = TRUE, fun=exp)
      } else if (i=="apob_mgdl"){
        Predict(model, apob_mgdl=c(dd$limits$apob_mgdl[2],dd$limits$apob_mgdl[2]+50),ref.zero = TRUE, fun=exp)
      }
      
      if (q!="model3"){
      pred[c("yhat","lower","upper")]<-exp(log(pred[c("yhat","lower","upper")])/get(paste0("lambda.",i)))
      }else {}
     
      
      fig4[paste(i,k),paste0(c("yhat","lower","upper"),q)]<-pred[2,c("yhat","lower","upper")] 
      fig4[paste(i,k),paste0("pval",q)]<-pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2
    
    print(paste(i,k,q,"done!"))
  }}

fig4

fig4$var<-factor(rownames(fig4),levels=rownames(fig4))


#Table

fig4_table <- tibble(cbind(fig4$var,fig4[grepl("yhat",colnames(fig4))|grepl("pval",colnames(fig4))])) %>% flextable()

for (i in c("model1","model2","model3","model4","model5","model6","model7")){
fig4_table<- compose(fig4_table,j=paste0("yhat",i), value=as_paragraph(round(fig4[,paste0("yhat",i)], digits = 1), " (", round(fig4[,paste0("lower",i)], digits = 1),"-", round(fig4[,paste0("upper",i)], digits = 1), ")"))
fig4_table<- compose(fig4_table,j=paste0("pval",i), value=as_paragraph(signif.p(fig4[,paste0("pval",i)]),as_sup(signif.p(fig4[,paste0("pval",i)],exp=T))))

}
fig4_table<-fontsize(fig4_table, size = 10, part = "body")
fig4_table<-align(autofit(fig4_table, part="all"),align="left", part = "all")
fig4_table


#figure

group.colors<-rep(c("darkblue","darkcyan", "darkorchid4","darkcyan","darkorchid4"),2)

fig4_plot<-list()
for (i in c("model1","model2","model3","model4","model5","model6","model7")){
  fig4_plot[i]<-list(ggplot(fig4, aes_string(y="var",x=paste0("yhat",i), xmin =paste0("lower",i), xmax = paste0("upper",i)))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme( axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_errorbar(width=0.5, cex=0.5,color="black", position=position_dodge(width=0.1))+ 
  geom_point(size=2,colour=group.colors)+
  scale_colour_grey(start = 0.75, end = 0)+
  scale_y_discrete(limits=rev)+scale_x_continuous(lim=c(0.5,3),breaks = c(1,2,3)))
}


fargs  <- formals(plot_grid)
fargs$ncol           <- 4
formals(plot_grid) <- fargs


figure4<-do.call(plot_grid,fig4_plot)


toffice(figure = figure4, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = TRUE, width =3.5, height = 5, devsize = FALSE, units = "in")




# Figure Sx: sensitivity analysis different adjustment (bmi, diabetes, ascvd) ----------------------------------


fig4<-data.frame()
for (k in c("PAD","AMI"))
  for (q in c("model1","model2","model3","model4")){
    
    dat<-cgps[c("apob_mgdl","remncholcalc","ldl","VLDLP","LDLP","sex","birthdate","obushnr","PAD","AMI","stopagePAD","stopageAMI","PAD_comprisk","AMI_comprisk","followupPAD","followupAMI","startage","systolic","diastolic","cum_smoking","status_smoking","chol", "XSVLDLP","trig","hdl","P1methyl2pyrrolidone","bmi","ascvd_pre","dmall_pre")]
    
    
    
    for (i in c("apob_mgdl","remncholcalc","ldl","VLDLP","LDLP")){
      
      
      adjust<-paste0("+", if(i=="remncholcalc") {"ldl"} else if(i=="ldl") {"remncholcalc"} else if(i=="VLDLP" ) {"LDLP"} else if(i=="LDLP") {"VLDLP"},
                     "+strat(sex)+birthdate+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)",
                     if (q=="model2"){"+bmi"} else if (q=="model3") {"+dmall_pre"} else if (q=="model4") {"+ascvd_pre"} else {""})
      
      
      model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,adjust)) ,data=dat,x=TRUE, y=TRUE)  
      
      
      pred<- if (i=="remncholcalc") {
        Predict(model,remncholcalc=c(dd$limits$remncholcalc[2],dd$limits$remncholcalc[2]+1), ref.zero = TRUE, fun=exp)
      } else if (i=="ldl"){
        Predict(model, ldl=c(dd$limits$ldl[2],dd$limits$ldl[2]+1),ref.zero = TRUE, fun=exp)
      } else if (i=="VLDLP") {
        Predict(model, VLDLP=c(dd$limits$VLDLP[2],dd$limits$VLDLP[2]+100),ref.zero = TRUE, fun=exp)
      } else if (i=="LDLP"){
        Predict(model, LDLP=c(dd$limits$LDLP[2],dd$limits$LDLP[2]+100),ref.zero = TRUE, fun=exp)
      } else if (i=="apob_mgdl"){
        Predict(model, apob_mgdl=c(dd$limits$apob_mgdl[2],dd$limits$apob_mgdl[2]+50),ref.zero = TRUE, fun=exp)
      }
      
      pred[c("yhat","lower","upper")]<-exp(log(pred[c("yhat","lower","upper")])/get(paste0("lambda.",i)))
      
      fig4[paste(i,k),paste0(c("yhat","lower","upper"),q)]<-pred[2,c("yhat","lower","upper")] 
      fig4[paste(i,k),paste0("pval",q)]<-pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2
      
      print(paste(i,k,q,"done!"))
    }}

fig4

fig4$var<-factor(rownames(fig4),levels=rownames(fig4))


#Table

fig4_table <- tibble(cbind(fig4$var,fig4[grepl("yhat",colnames(fig4))|grepl("pval",colnames(fig4))])) %>% flextable()

for (i in c("model1","model2","model3","model4")){
  fig4_table<- compose(fig4_table,j=paste0("yhat",i), value=as_paragraph(round(fig4[,paste0("yhat",i)], digits = 1), " (", round(fig4[,paste0("lower",i)], digits = 1),"-", round(fig4[,paste0("upper",i)], digits = 1), ")"))
  fig4_table<- compose(fig4_table,j=paste0("pval",i), value=as_paragraph(signif.p(fig4[,paste0("pval",i)]),as_sup(signif.p(fig4[,paste0("pval",i)],exp=T))))
  
}
fig4_table<-fontsize(fig4_table, size = 10, part = "body")
fig4_table<-align(autofit(fig4_table, part="all"),align="left", part = "all")
fig4_table


#figure

group.colors<-rep(c("darkblue","darkcyan", "darkorchid4","darkcyan","darkorchid4"),2)

fig4_plot<-list()
for (i in c("model1","model2","model3","model4")){
  fig4_plot[i]<-list(ggplot(fig4, aes_string(y="var",x=paste0("yhat",i), xmin =paste0("lower",i), xmax = paste0("upper",i)))+
                       labs(x="Hazard ratio (95% CI)", y=NULL)+ 
                       theme_nothing()+
                       theme( axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
                       geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_errorbar(width=0.5, cex=0.5,color="black", position=position_dodge(width=0.1))+ 
                       geom_point(size=2,colour=group.colors)+
                       scale_colour_grey(start = 0.75, end = 0)+
                       scale_y_discrete(limits=rev)+scale_x_continuous(lim=c(0.5,3),breaks = c(1,2,3)))
}


fargs  <- formals(plot_grid)
fargs$ncol           <- 4
formals(plot_grid) <- fargs


figure4<-do.call(plot_grid,fig4_plot)


toffice(figure = figure4, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = TRUE, width =3.5, height = 2.5, devsize = FALSE, units = "in")




# Figure Sx: sensitivity analysis different adjustment (TG, HDL...) ----------------------------------


fig4<-data.frame()
for (k in c("PAD","AMI"))
  for (q in c("model1","model2","model4")){
    
    dat<-cgps[c("apob_mgdl","remncholcalc","ldl","VLDLP","LDLP","sex","birthdate","obushnr","PAD","AMI","stopagePAD","stopageAMI","PAD_comprisk","AMI_comprisk","followupPAD","followupAMI","startage","systolic","diastolic","cum_smoking","status_smoking","chol", "XSVLDLP","trig","hdl","P1methyl2pyrrolidone","bmi","ascvd_pre","dmall_pre")]
    

    
    for (i in c("apob_mgdl","remncholcalc","ldl","VLDLP","LDLP")){
      
      
      adjust<-paste0("+", if(i=="remncholcalc"& q!="model4") {"ldl"} else if(i=="ldl"& q!="model4") {"remncholcalc"} else if(i=="VLDLP" & q!="model4") {"LDLP"} else if(i=="LDLP"& q!="model4") {"VLDLP"},
                     "+strat(sex)+birthdate+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)",
                     if (q=="model2"){"+log(trig)"} else if (q=="model3") {"+hdl"} else {""})
      
      
        model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,adjust)) ,data=dat,x=TRUE, y=TRUE)  
      
      
      pred<- if (i=="remncholcalc") {
        Predict(model,remncholcalc=c(dd$limits$remncholcalc[2],dd$limits$remncholcalc[2]+1), ref.zero = TRUE, fun=exp)
      } else if (i=="ldl"){
        Predict(model, ldl=c(dd$limits$ldl[2],dd$limits$ldl[2]+1),ref.zero = TRUE, fun=exp)
      } else if (i=="VLDLP") {
        Predict(model, VLDLP=c(dd$limits$VLDLP[2],dd$limits$VLDLP[2]+100),ref.zero = TRUE, fun=exp)
      } else if (i=="LDLP"){
        Predict(model, LDLP=c(dd$limits$LDLP[2],dd$limits$LDLP[2]+100),ref.zero = TRUE, fun=exp)
      } else if (i=="apob_mgdl"){
        Predict(model, apob_mgdl=c(dd$limits$apob_mgdl[2],dd$limits$apob_mgdl[2]+50),ref.zero = TRUE, fun=exp)
      }
      
      pred[c("yhat","lower","upper")]<-exp(log(pred[c("yhat","lower","upper")])/get(paste0("lambda.",i)))
      
      fig4[paste(i,k),paste0(c("yhat","lower","upper"),q)]<-pred[2,c("yhat","lower","upper")] 
      fig4[paste(i,k),paste0("pval",q)]<-pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2
      
      print(paste(i,k,q,"done!"))
    }}

fig4

fig4$var<-factor(rownames(fig4),levels=rownames(fig4))


#Table

fig4_table <- tibble(cbind(fig4$var,fig4[grepl("yhat",colnames(fig4))|grepl("pval",colnames(fig4))])) %>% flextable()

for (i in c("model1","model2","model4")){
  fig4_table<- compose(fig4_table,j=paste0("yhat",i), value=as_paragraph(round(fig4[,paste0("yhat",i)], digits = 1), " (", round(fig4[,paste0("lower",i)], digits = 1),"-", round(fig4[,paste0("upper",i)], digits = 1), ")"))
  fig4_table<- compose(fig4_table,j=paste0("pval",i), value=as_paragraph(signif.p(fig4[,paste0("pval",i)]),as_sup(signif.p(fig4[,paste0("pval",i)],exp=T))))
  
}
fig4_table<-fontsize(fig4_table, size = 10, part = "body")
fig4_table<-align(autofit(fig4_table, part="all"),align="left", part = "all")
fig4_table



#figure

group.colors<-rep(c("darkblue","darkcyan", "darkorchid4","darkcyan","darkorchid4"),2)

fig4_plot<-list()
for (i in c("model1","model2","model3","model4")){
  fig4_plot[i]<-list(ggplot(fig4, aes_string(y="var",x=paste0("yhat",i), xmin =paste0("lower",i), xmax = paste0("upper",i)))+
                       labs(x="Hazard ratio (95% CI)", y=NULL)+ 
                       theme_nothing()+
                       theme( axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
                       geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_errorbar(width=0.5, cex=0.5,color="black", position=position_dodge(width=0.1))+ 
                       geom_point(size=2,colour=group.colors)+
                       scale_colour_grey(start = 0.75, end = 0)+
                       scale_y_discrete(limits=rev)+scale_x_continuous(lim=c(0.5,3),breaks = c(1,2,3)))
}


fargs  <- formals(plot_grid)
fargs$ncol           <- 4
formals(plot_grid) <- fargs


figure4<-do.call(plot_grid,fig4_plot)


toffice(figure = figure4, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = TRUE, width =3.5, height = 2.5, devsize = FALSE, units = "in")




# Figure SX: sensitivity analysis no ASCVD ------------------------------------------------------


fig4<-data.frame()
dat<-subset(cgps,ascvd_pre==F)
for (k in c("PAD","AMI"))
  for (i in c( "apob_mgdl","remncholcalc","ldl","VLDLP","LDLP")){
    
    fig4[paste(i,k),"var"]<-paste0(i,k)
    
    p_years<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~1")) ,data=subset(dat,!is.na(get(i))),scale=1)
    
    fig4[paste(i,k),"N"]<-p_years$n
    fig4[paste(i,k),"N_events"]<-p_years$event
    fig4[paste(i,k),"events_pyears"]<-round(p_years$event/p_years$pyear*1000,1)
    
     model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,"+",
                                 if(i=="remncholcalc") {"ldl"} else if(i=="ldl") {"remncholcalc"} else if(i=="VLDLP") {"LDLP"} else if(i=="LDLP") {"VLDLP"} ,
                                 "+strat(sex)+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate")) ,data=dat,x=TRUE, y=TRUE)  
    
    
    pred<- if (i=="remncholcalc") {
      Predict(model,remncholcalc=c(dd$limits$remncholcalc[2],dd$limits$remncholcalc[2]+1), ref.zero = TRUE, fun=exp)
    } else if (i=="ldl"){
      Predict(model, ldl=c(dd$limits$ldl[2],dd$limits$ldl[2]+1),ref.zero = TRUE, fun=exp)
    } else if (i=="VLDLP") {
      Predict(model, VLDLP=c(dd$limits$VLDLP[2],dd$limits$VLDLP[2]+100),ref.zero = TRUE, fun=exp)
    } else if (i=="LDLP"){
      Predict(model, LDLP=c(dd$limits$LDLP[2],dd$limits$LDLP[2]+100),ref.zero = TRUE, fun=exp)
    } else if (i=="apob_mgdl"){
      Predict(model, apob_mgdl=c(dd$limits$apob_mgdl[2],dd$limits$apob_mgdl[2]+50),ref.zero = TRUE, fun=exp)
    }
    
    pred[c("yhat","lower","upper")]<-exp(log(pred[c("yhat","lower","upper")])/get(paste0("lambda.",i)))
    
    fig4[paste(i,k),c("yhat","lower","upper")]<-pred[2,c("yhat","lower","upper")] 
    fig4[paste(i,k),"pval"]<-pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2
    
    

  }

rm(dat)


fig4$var<-factor(fig4$var,levels=fig4$var)


#Creating table
fig4_table <- tibble(fig4) %>% flextable(col_keys = colnames(fig4)[grepl("lower",colnames(fig4)) ==F & grepl("upper",colnames(fig4))==F])
fig4_table<- compose(fig4_table,j=paste0("yhat"), value=as_paragraph(round(fig4[,paste0("yhat")], digits = 1), " (", round(fig4[,paste0("lower")], digits = 1),"-", round(fig4[,paste0("upper")], digits = 1), ")"))
fig4_table<- compose(fig4_table,j=paste0("pval"), value=as_paragraph(signif.p(fig4[,paste0("pval")]),as_sup(signif.p(fig4[,paste0("pval")],exp=T))))

fig4_table<-fontsize(fig4_table, size = 11, part = "body")
fig4_table<-align(autofit(fig4_table, part="all"),align="left", part = "all")
fig4_table

save_as_pptx(fig4_table, path =  "Temp.pptx")




#figure

group.colors <- rep(c("darkblue","darkcyan","darkorchid4","darkcyan","darkorchid4"),2)

figure4<-ggplot(fig4, aes(y=var,x=yhat, xmin =lower, xmax = upper))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme(axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_point(size=3,colour=group.colors)+
  geom_errorbar(width=0.5, cex=0.7,color=group.colors, position=position_dodge(width=0.1))+ 
  scale_colour_grey(start = 0.75, end = 0)+
  scale_y_discrete(limits=rev)  +coord_cartesian( xlim = c(0.75,3))



toffice(figure = figure4, format = "pptx", title = "Copenhagen General Population Study",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = TRUE, width = 2.5, height = 3.3,devsize = FALSE, units = "in")







# Figure SX: sensitivity analysis excluding follow-up time ------------------------------------------------------

fig4<-data.frame()
for (q in c(0:6)) 
  for (k in c("PAD","AMI")){
    dat<-subset(cgps,get(paste0("followup",k))>q)
    for (i in c( "apob_mgdl","remncholcalc","ldl","VLDLP","LDLP")){
      fig4[paste(i,k),"var"]<-paste0(i,k)
      p_years<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~1")) ,data=subset(dat,!is.na(get(i))),scale=1)
      
      fig4[paste(i,k),paste0("N",q)]<-p_years$n
      fig4[paste(i,k),paste0("N_events",q)]<-p_years$event
      fig4[paste(i,k),paste0("events_pyears",q)]<-round(p_years$event/p_years$pyear*1000,1)
      
      
      model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,"+",
                                   if(i=="remncholcalc") {"ldl"} else if(i=="ldl") {"remncholcalc"} else if(i=="VLDLP") {"LDLP"} else if(i=="LDLP") {"VLDLP"} ,
                                   "+strat(sex)+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate")) ,data=dat,x=TRUE, y=TRUE)  
      
      
      pred<- if (i=="remncholcalc") {
        Predict(model,remncholcalc=c(dd$limits$remncholcalc[2],dd$limits$remncholcalc[2]+1), ref.zero = TRUE, fun=exp)
      } else if (i=="ldl"){
        Predict(model, ldl=c(dd$limits$ldl[2],dd$limits$ldl[2]+1),ref.zero = TRUE, fun=exp)
      } else if (i=="VLDLP") {
        Predict(model, VLDLP=c(dd$limits$VLDLP[2],dd$limits$VLDLP[2]+100),ref.zero = TRUE, fun=exp)
      } else if (i=="LDLP"){
        Predict(model, LDLP=c(dd$limits$LDLP[2],dd$limits$LDLP[2]+100),ref.zero = TRUE, fun=exp)
      } else if (i=="apob_mgdl"){
        Predict(model, apob_mgdl=c(dd$limits$apob_mgdl[2],dd$limits$apob_mgdl[2]+50),ref.zero = TRUE, fun=exp)
      }
      
      pred[c("yhat","lower","upper")]<-exp(log(pred[c("yhat","lower","upper")])/get(paste0("lambda.",i)))
      
      
      fig4[paste(i,k),paste0(c("yhat","lower","upper"),q)]<-pred[2,c("yhat","lower","upper")] 
      fig4[paste(i,k),paste0("pval",q)]<-pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2
      
      
    }
  }


fig4$var<-factor(fig4$var,levels=fig4$var)

#Creating table
fig4_table <- tibble(fig4) %>% flextable(col_keys = colnames(fig4)[grepl("lower",colnames(fig4)) ==F & grepl("upper",colnames(fig4))==F])
for (i in 0:6){
  fig4_table<- compose(fig4_table,j=paste0("yhat",i), value=as_paragraph(round(fig4[,paste0("yhat",i)], digits = 1), " (", round(fig4[,paste0("lower",i)], digits = 1),"-", round(fig4[,paste0("upper",i)], digits = 1), ")"))
  fig4_table<- compose(fig4_table,j=paste0("pval",i), value=as_paragraph(signif.p(fig4[,paste0("pval",i)]),as_sup(signif.p(fig4[,paste0("pval",i)],exp=T))))
  
}
fig4_table<-fontsize(fig4_table, size = 10, part = "body")
fig4_table<-align(autofit(fig4_table, part="all"),align="left", part = "all")
fig4_table



#figure

group.colors<-rep(c("darkblue","darkcyan", "darkorchid4","darkcyan","darkorchid4"),2)

fig4_plot<-list()
for (i in 0:6){
  fig4_plot[paste(i)]<-list(ggplot(fig4, aes_string(y="var",x=paste0("yhat",i), xmin =paste0("lower",i), xmax = paste0("upper",i)))+
                              labs(x="Hazard ratio (95% CI)", y=NULL)+ 
                              theme_nothing()+
                              theme( axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
                              geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_errorbar(width=0.5, cex=0.5,color="black", position=position_dodge(width=0.1))+ 
                              geom_point(size=2,colour=group.colors)+
                              scale_colour_grey(start = 0.75, end = 0)+
                              scale_y_discrete(limits=rev)+scale_x_continuous(lim=c(0.5,4),breaks = c(1,2,3,4)))
}


fargs  <- formals(plot_grid)
fargs$ncol           <- 4
formals(plot_grid) <- fargs


figure4<-do.call(plot_grid,fig4_plot)


toffice(figure = figure4, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = TRUE, width =3.5, height = 5, devsize = FALSE, units = "in")







# Figure SX: sensitivity analysis only NMR dataset ------------------------------------------------------


fig4<-data.frame()
dat<-subset(cgps,!is.na(VLDLP))
for (k in c("PAD","AMI"))
for (i in c( "apob_mgdl","remncholcalc","ldl","VLDLP","LDLP")){
  fig4[paste(i,k),"var"]<-paste0(i,k)
  p_years<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~1")) ,data=subset(dat,!is.na(get(i))),scale=1)
  
  fig4[paste(i,k),"N"]<-p_years$n
  fig4[paste(i,k),"N_events"]<-p_years$event
  fig4[paste(i,k),"events_pyears"]<-round(p_years$event/p_years$pyear*1000,1)
  
  model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,"+",
                               if(i=="remncholcalc") {"ldl"} else if(i=="ldl") {"remncholcalc"} else if(i=="VLDLP") {"LDLP"} else if(i=="LDLP") {"VLDLP"} ,
                               "+strat(sex)+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate")) ,data=dat,x=TRUE, y=TRUE)  
  
  
  pred<- if (i=="remncholcalc") {
    Predict(model,remncholcalc=c(dd$limits$remncholcalc[2],dd$limits$remncholcalc[2]+1), ref.zero = TRUE, fun=exp)
  } else if (i=="ldl"){
    Predict(model, ldl=c(dd$limits$ldl[2],dd$limits$ldl[2]+1),ref.zero = TRUE, fun=exp)
  } else if (i=="VLDLP") {
    Predict(model, VLDLP=c(dd$limits$VLDLP[2],dd$limits$VLDLP[2]+100),ref.zero = TRUE, fun=exp)
  } else if (i=="LDLP"){
    Predict(model, LDLP=c(dd$limits$LDLP[2],dd$limits$LDLP[2]+100),ref.zero = TRUE, fun=exp)
  } else if (i=="apob_mgdl"){
    Predict(model, apob_mgdl=c(dd$limits$apob_mgdl[2],dd$limits$apob_mgdl[2]+50),ref.zero = TRUE, fun=exp)
  }
  
  pred[c("yhat","lower","upper")]<-exp(log(pred[c("yhat","lower","upper")])/get(paste0("lambda.",i)))
  

  fig4[paste(i,k),c("yhat","lower","upper")]<-pred[2,c("yhat","lower","upper")] 
  fig4[paste(i,k),"pval"]<-pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2
  

}

rm(dat)
fig4

fig4$var<-factor(fig4$var,levels=fig4$var)



#Creating table
fig4_table <- tibble(fig4) %>% flextable(col_keys = colnames(fig4)[grepl("lower",colnames(fig4)) ==F & grepl("upper",colnames(fig4))==F])
fig4_table<- compose(fig4_table,j=paste0("yhat"), value=as_paragraph(round(fig4[,paste0("yhat")], digits = 1), " (", round(fig4[,paste0("lower")], digits = 1),"-", round(fig4[,paste0("upper")], digits = 1), ")"))
fig4_table<- compose(fig4_table,j=paste0("pval"), value=as_paragraph(signif.p(fig4[,paste0("pval")]),as_sup(signif.p(fig4[,paste0("pval")],exp=T))))

fig4_table<-fontsize(fig4_table, size = 11, part = "body")
fig4_table<-align(autofit(fig4_table, part="all"),align="left", part = "all")
fig4_table

save_as_pptx(fig4_table, path =  "Temp.pptx")



#figure

group.colors <- rep(c("darkblue","darkcyan","darkorchid4","darkcyan","darkorchid4"),2)

figure4<-ggplot(fig4, aes(y=var,x=yhat, xmin =lower, xmax = upper))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme(axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_point(size=3,colour=group.colors)+
  geom_errorbar(width=0.5, cex=0.7,color=group.colors, position=position_dodge(width=0.1))+ 
  scale_colour_grey(start = 0.75, end = 0)+
  scale_y_discrete(limits=rev)  +coord_cartesian( xlim = c(0.75,3))



toffice(figure = figure4, format = "pptx", title = "Copenhagen General Population Study",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = T, width = 2.5, height = 3.3,devsize = FALSE, units = "in")







# Figure SX: sensitivity analysis no ASCVD ------------------------------------------------------


cgps$PAD_AMI_int<-interaction(cgps$PAD,cgps$AMI)

fig4<-data.frame()
dat<-subset(cgps,PAD_AMI_int!="TRUE.TRUE")
for (k in c("PAD","AMI"))
  for (i in c( "apob_mgdl","remncholcalc","ldl","VLDLP","LDLP")){
    
    fig4[paste(i,k),"var"]<-paste0(i,k)
    
    p_years<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~1")) ,data=subset(dat,!is.na(get(i))),scale=1)
    
    fig4[paste(i,k),"N"]<-p_years$n
    fig4[paste(i,k),"N_events"]<-p_years$event
    fig4[paste(i,k),"events_pyears"]<-round(p_years$event/p_years$pyear*1000,1)
    
    model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,"+",
                                 if(i=="remncholcalc") {"ldl"} else if(i=="ldl") {"remncholcalc"} else if(i=="VLDLP") {"LDLP"} else if(i=="LDLP") {"VLDLP"} ,
                                 "+strat(sex)+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate")) ,data=dat,x=TRUE, y=TRUE)  
    
    
    pred<- if (i=="remncholcalc") {
      Predict(model,remncholcalc=c(dd$limits$remncholcalc[2],dd$limits$remncholcalc[2]+1), ref.zero = TRUE, fun=exp)
    } else if (i=="ldl"){
      Predict(model, ldl=c(dd$limits$ldl[2],dd$limits$ldl[2]+1),ref.zero = TRUE, fun=exp)
    } else if (i=="VLDLP") {
      Predict(model, VLDLP=c(dd$limits$VLDLP[2],dd$limits$VLDLP[2]+100),ref.zero = TRUE, fun=exp)
    } else if (i=="LDLP"){
      Predict(model, LDLP=c(dd$limits$LDLP[2],dd$limits$LDLP[2]+100),ref.zero = TRUE, fun=exp)
    } else if (i=="apob_mgdl"){
      Predict(model, apob_mgdl=c(dd$limits$apob_mgdl[2],dd$limits$apob_mgdl[2]+50),ref.zero = TRUE, fun=exp)
    }
    
    pred[c("yhat","lower","upper")]<-exp(log(pred[c("yhat","lower","upper")])/get(paste0("lambda.",i)))
    
    fig4[paste(i,k),c("yhat","lower","upper")]<-pred[2,c("yhat","lower","upper")] 
    fig4[paste(i,k),"pval"]<-pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2
    
    
    
  }

rm(dat)


fig4$var<-factor(fig4$var,levels=fig4$var)


#Creating table
fig4_table <- tibble(fig4) %>% flextable(col_keys = colnames(fig4)[grepl("lower",colnames(fig4)) ==F & grepl("upper",colnames(fig4))==F])
fig4_table<- compose(fig4_table,j=paste0("yhat"), value=as_paragraph(round(fig4[,paste0("yhat")], digits = 1), " (", round(fig4[,paste0("lower")], digits = 1),"-", round(fig4[,paste0("upper")], digits = 1), ")"))
fig4_table<- compose(fig4_table,j=paste0("pval"), value=as_paragraph(signif.p(fig4[,paste0("pval")]),as_sup(signif.p(fig4[,paste0("pval")],exp=T))))

fig4_table<-fontsize(fig4_table, size = 11, part = "body")
fig4_table<-align(autofit(fig4_table, part="all"),align="left", part = "all")
fig4_table


#figure

group.colors <- rep(c("darkblue","darkcyan","darkorchid4","darkcyan","darkorchid4"),2)

figure4<-ggplot(fig4, aes(y=var,x=yhat, xmin =lower, xmax = upper))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme(axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_point(size=3,colour=group.colors)+
  geom_errorbar(width=0.5, cex=0.7,color=group.colors, position=position_dodge(width=0.1))+ 
  scale_colour_grey(start = 0.75, end = 0)+
  scale_y_discrete(limits=rev)  +coord_cartesian( xlim = c(0.75,3))



toffice(figure = figure4, format = "pptx", title = "Copenhagen General Population Study",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = TRUE, width = 2.5, height = 3.3,devsize = FALSE, units = "in")







# Figure SX: sensitivity analysis men and women ------------------------------------------------------


fig4<-data.frame()
for (q in levels(cgps$sex)) {
dat<-subset(cgps,sex==q)
for (k in c("PAD","AMI"))
  for (i in c( "apob_mgdl","remncholcalc","ldl","VLDLP","LDLP")){
    fig4[paste(i,k),"var"]<-paste0(i,k)
    p_years<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~1")) ,data=subset(dat,!is.na(get(i))),scale=1)
    
    fig4[paste(i,k),paste0("N",q)]<-p_years$n
    fig4[paste(i,k),paste0("N_events",q)]<-p_years$event
    fig4[paste(i,k),paste0("events_pyears",q)]<-round(p_years$event/p_years$pyear*1000,1)
    
    
    model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,"+",
                                 if(i=="remncholcalc") {"ldl"} else if(i=="ldl") {"remncholcalc"} else if(i=="VLDLP") {"LDLP"} else if(i=="LDLP") {"VLDLP"} ,
                                 "+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate")) ,data=dat,x=TRUE, y=TRUE)  
    
    
    pred<- if (i=="remncholcalc") {
      Predict(model,remncholcalc=c(dd$limits$remncholcalc[2],dd$limits$remncholcalc[2]+1), ref.zero = TRUE, fun=exp)
    } else if (i=="ldl"){
      Predict(model, ldl=c(dd$limits$ldl[2],dd$limits$ldl[2]+1),ref.zero = TRUE, fun=exp)
    } else if (i=="VLDLP") {
      Predict(model, VLDLP=c(dd$limits$VLDLP[2],dd$limits$VLDLP[2]+100),ref.zero = TRUE, fun=exp)
    } else if (i=="LDLP"){
      Predict(model, LDLP=c(dd$limits$LDLP[2],dd$limits$LDLP[2]+100),ref.zero = TRUE, fun=exp)
    } else if (i=="apob_mgdl"){
      Predict(model, apob_mgdl=c(dd$limits$apob_mgdl[2],dd$limits$apob_mgdl[2]+50),ref.zero = TRUE, fun=exp)
    }
    
    pred[c("yhat","lower","upper")]<-exp(log(pred[c("yhat","lower","upper")])/get(paste0("lambda.",i)))
    
    
    fig4[paste(i,k),paste0(c("yhat","lower","upper"),q)]<-pred[2,c("yhat","lower","upper")] 
    fig4[paste(i,k),paste0("pval",q)]<-pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2
    
    
  }
}
rm(dat)


fig4$var<-factor(fig4$var,levels=fig4$var)

#Creating table
fig4_table <- tibble(fig4) %>% flextable(col_keys = colnames(fig4)[grepl("lower",colnames(fig4)) ==F & grepl("upper",colnames(fig4))==F])
for (i in levels(cgps$sex)){
  fig4_table<- compose(fig4_table,j=paste0("yhat",i), value=as_paragraph(round(fig4[,paste0("yhat",i)], digits = 1), " (", round(fig4[,paste0("lower",i)], digits = 1),"-", round(fig4[,paste0("upper",i)], digits = 1), ")"))
  fig4_table<- compose(fig4_table,j=paste0("pval",i), value=as_paragraph(signif.p(fig4[,paste0("pval",i)]),as_sup(signif.p(fig4[,paste0("pval",i)],exp=T))))
  
}
fig4_table<-fontsize(fig4_table, size = 11, part = "body")
fig4_table<-align(autofit(fig4_table, part="all"),align="left", part = "all")
fig4_table


#figure

group.colors<-rep(c("darkblue","darkcyan", "darkorchid4","darkcyan","darkorchid4"),2)

fig4_plot<-list()
for (i in levels(cgps$sex)){
  fig4_plot[i]<-list(ggplot(fig4, aes_string(y="var",x=paste0("yhat",i), xmin =paste0("lower",i), xmax = paste0("upper",i)))+
                       labs(x="Hazard ratio (95% CI)", y=NULL)+ 
                       theme_nothing()+
                       theme( axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
                       geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_errorbar(width=0.5, cex=0.5,color="black", position=position_dodge(width=0.1))+ 
                       geom_point(size=2,colour=group.colors)+
                       scale_colour_grey(start = 0.75, end = 0)+
                       scale_y_discrete(limits=rev)+scale_x_continuous(lim=c(0.5,3),breaks = c(1,2,3)))
}


fargs  <- formals(plot_grid)
fargs$ncol           <- 2
formals(plot_grid) <- fargs


figure4<-do.call(plot_grid,fig4_plot)

toffice(figure = figure4, format = "pptx", title = "Copenhagen General Population Study",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = T, width = 3, height = 3.3,devsize = FALSE, units = "in")




fig4<-data.frame()
# P values for interaction
for (k in c("PAD","AMI"))
  for (i in c( "apob_mgdl","remncholcalc","ldl","VLDLP","LDLP")){
  
    dat<-cgps
    
    model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,"*sex+",
                                 if(i=="remncholcalc") {"ldl"} else if(i=="ldl") {"remncholcalc"} else if(i=="VLDLP") {"LDLP"} else if(i=="LDLP") {"VLDLP"} ,
                                 "+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate")) ,data=dat,x=TRUE, y=TRUE)  
    
    
s<-length(model$coefficients)

    
    fig4[paste(i,k),paste0("pval interaction")]<-pnorm(abs(model$coef[s])/sqrt(model$var[s,s]),lower.tail=F)*2
    
    
  }

fig4


# Figure SX: sensitivity analysis by PAD ICD-10 code ------------------------------------------------------


fig4<-data.frame()
for (q in c("DI702","DI702A","DI739","DI739A","DI739C","K"))
  for (k in c("PAD")){
    
    dat<-subset(cgps,get(paste0("followup",k))>0)
    
    q<-if(q=="K"){grep("K",unique(cgps$c_diag),value=T)} else {q}
    
    dat[dat$PAD==T & dat$c_diag %in% q==F,]$stopagePAD<- dat[dat$PAD==T & dat$c_diag %in% q==F,]$stopagemort
    
    dat[dat$c_diag %in% q==F |is.na(dat$c_diag),]$PAD<-F

    q<-q[1]
      
    if (sum(dat$PAD==T)<5) {next}
    
    for (i in c( "apob_mgdl","remncholcalc","ldl","VLDLP","LDLP")){
      fig4[paste(i,k),"var"]<-paste0(i,k)
      p_years<-pyears(as.formula(paste0("Surv(startage,stopage",k,",",k,")~1")) ,data=subset(dat,!is.na(get(i))),scale=1)
      
      fig4[paste(i,k),paste0("N",q)]<-p_years$n
      fig4[paste(i,k),paste0("N_events",q)]<-p_years$event
      fig4[paste(i,k),paste0("events_pyears",q)]<-round(p_years$event/p_years$pyear*1000,1)
      
      
      model<-cph(as.formula(paste0("Surv(startage,stopage",k,",",k,")~",i,"+",
                                   if(i=="remncholcalc") {"ldl"} else if(i=="ldl") {"remncholcalc"} else if(i=="VLDLP") {"LDLP"} else if(i=="LDLP") {"VLDLP"} ,
                                   "+strat(sex)+systolic+diastolic+strat(status_smoking)+rcs(cum_smoking, 3)+birthdate"
                                   )) ,data=dat,x=TRUE, y=TRUE)  
      
      
      pred<- if (i=="remncholcalc") {
        Predict(model,remncholcalc=c(dd$limits$remncholcalc[2],dd$limits$remncholcalc[2]+1), ref.zero = TRUE, fun=exp)
      } else if (i=="ldl"){
        Predict(model, ldl=c(dd$limits$ldl[2],dd$limits$ldl[2]+1),ref.zero = TRUE, fun=exp)
      } else if (i=="VLDLP") {
        Predict(model, VLDLP=c(dd$limits$VLDLP[2],dd$limits$VLDLP[2]+100),ref.zero = TRUE, fun=exp)
      } else if (i=="LDLP"){
        Predict(model, LDLP=c(dd$limits$LDLP[2],dd$limits$LDLP[2]+100),ref.zero = TRUE, fun=exp)
      } else if (i=="apob_mgdl"){
        Predict(model, apob_mgdl=c(dd$limits$apob_mgdl[2],dd$limits$apob_mgdl[2]+50),ref.zero = TRUE, fun=exp)
      }
      
      pred[c("yhat","lower","upper")]<-exp(log(pred[c("yhat","lower","upper")])/get(paste0("lambda.",i)))
      
      
      fig4[paste(i,k),paste0(c("yhat","lower","upper"),q)]<-pred[2,c("yhat","lower","upper")] 
      fig4[paste(i,k),paste0("pval",q)]<-pnorm(abs(model$coef[1])/sqrt(model$var[1,1]),lower.tail=F)*2
      
      
    }
  }


fig4$var<-factor(fig4$var,levels=fig4$var)


#Creating table
fig4_table <- tibble(fig4) %>% flextable(col_keys = colnames(fig4)[grepl("lower",colnames(fig4)) ==F & grepl("upper",colnames(fig4))==F])
for (i in c("DI702","DI702A","DI739","DI739A","DI739C","KNFQ19")){
  fig4_table<- compose(fig4_table,j=paste0("yhat",i), value=as_paragraph(round(fig4[,paste0("yhat",i)], digits = 1), " (", round(fig4[,paste0("lower",i)], digits = 1),"-", round(fig4[,paste0("upper",i)], digits = 1), ")"))
  fig4_table<- compose(fig4_table,j=paste0("pval",i), value=as_paragraph(signif.p(fig4[,paste0("pval",i)]),as_sup(signif.p(fig4[,paste0("pval",i)],exp=T))))
  
}
fig4_table<-fontsize(fig4_table, size = 10, part = "body")
fig4_table<-align(autofit(fig4_table, part="all"),align="left", part = "all")
fig4_table

save_as_pptx(fig4_table, path =  "Temp.pptx")


#figure

group.colors<-rep(c("darkblue","darkcyan", "darkorchid4","darkcyan","darkorchid4"),1)

fig4_plot<-list()
for (i in c("DI702","DI702A","DI739","DI739A","DI739C","KNFQ19")){
  fig4_plot[paste(i)]<-list(ggplot(fig4, aes_string(y="var",x=paste0("yhat",i), xmin =paste0("lower",i), xmax = paste0("upper",i)))+
                              labs(x="Hazard ratio (95% CI)", y=NULL)+ 
                              theme_nothing()+
                              theme( axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
                              geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_errorbar(width=0.5, cex=0.5,color="black", position=position_dodge(width=0.1))+ 
                              geom_point(size=2,colour=group.colors)+
                              scale_colour_grey(start = 0.75, end = 0)+
                              scale_y_discrete(limits=rev)+scale_x_continuous(lim=c(0.5,4),breaks = c(1,2,3,4)))
}


fargs  <- formals(plot_grid)
fargs$ncol           <- 4
formals(plot_grid) <- fargs


figure4<-do.call(plot_grid,fig4_plot)


toffice(figure = figure4, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = TRUE, width =3.5, height = 5, devsize = FALSE, units = "in")







# Figure SX: indirect and direct effects ------------------------------------------------------



fig4<-data.frame()
for (k in c("PAD","CLTI","AMI"))
  for ( i in c("remncholcalc_log","ldl","VLDLP_log","LDLP")){
    
    dat<-cgps[!is.na(cgps[,i]),]
    
    
    obj<-cmest(data=subset(subset(dat, if (k=="PAD"){followupPAD>0} 
                                  else if (k=="AMI"){followupAMI>0}
                                  else if (k=="IS"){followupIS>0}
                                  else if (k=="ASCVD"){followupASCVD>0}
                                  else {sex!="hej"}), if (i=="hscrp_log"){hscrp>0} else {sex!="hej"}),
               model="rb", estimation = "paramfunc", inference = "delta", outcome=paste0("followup",k), event=k, 
               exposure = "apob_cut2", mediator=c(i), EMint=F, 
               basec= if (i=="VLDLP_log"| i=="LDLP"){c("startage","sex", "status_smoking", "cum_smoking", "birthdate", "systolic","diastolic")} else {
                 c("startage","sex", "status_smoking", "cum_smoking", "birthdate", "systolic","diastolic")},
               yreg="coxph", mreg=list("linear"),mval=list(1) , yval=TRUE, nboot=1, multimp = F)
    
    fig4[paste(i,k,c("indirect","direct","total")),c("yhat","lower","upper", "P.val")]<-as.data.frame(summary(obj)[c("effect.pe","effect.ci.low","effect.ci.high","effect.pval")])[c("Rpnie","Rpnde","Rte"), ]
    
    print(paste(i,k,"done!"))
    
  }

fig4$variable<-factor(rownames(fig4),levels=rownames(fig4))

#---Creating table
fig4_table <- tibble(fig4) %>% flextable(col_keys = colnames(fig4)[grepl("lower",colnames(fig4)) ==F & grepl("upper",colnames(fig4))==F])
for (i in ""){
  fig4_table<- compose(fig4_table,j=paste0("yhat",i), value=as_paragraph(round(fig4[,paste0("yhat",i)], digits = 1), " (", round(fig4[,paste0("lower",i)], digits = 1),"-", round(fig4[,paste0("upper",i)], digits = 1), ")"))
  fig4_table<- compose(fig4_table,j=paste0("P.val",i), value=as_paragraph(signif.p(fig4[,paste0("P.val",i)]),as_sup(signif.p(fig4[,paste0("P.val",i)],exp=T))))
  
}
fig4_table<-fontsize(fig4_table, size = 10, part = "body")
fig4_table<-align(autofit(fig4_table, part="all"),align="left", part = "all")
fig4_table


group.colors <- rep(c("darkcyan","darkblue","black","darkorchid4","darkblue","black"),6)

#figure

figure4<-ggplot(fig4, aes(y=variable,x=yhat, xmin =lower, xmax = upper))+
  labs(x="Hazard ratio (95% CI)", y=NULL)+ 
  theme_nothing()+
  theme(axis.line.x=element_line(size=1),axis.ticks.x=element_line(size=1),axis.ticks.length.x = unit(0.2, "cm"), axis.text.x =element_text(vjust=0.4), axis.title.x = element_text() , plot.margin=margin(0,3,0,3),plot.title =element_text(face="bold",size=10,hjust = 0.5),strip.text.x = element_text(size = 10))+
  geom_vline(aes(xintercept = 1), linetype = 2, size=0.35) + geom_point(size=2,colour=group.colors)+
  geom_errorbar(width=0.5, cex=0.75,color=group.colors, position=position_dodge(width=0.1))+ 
  scale_y_discrete(limits=rev)  +coord_cartesian( xlim = c(0.7,2))



toffice(figure = figure4, format = "pptx", title = "Copenhagen General Population Study",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = T, width = 2.5, height = 5,devsize = FALSE, units = "in")




# Figure SX: mediation analysis only NMR dataset-----------------------------------


dat<-subset(cgps,!is.na(VLDLP))
mediation_results<-data.frame()
for (k in c("PAD","CLTI","AMI"))
  for ( i in c("remncholcalc_log","ldl","VLDLP_log","LDLP")){
    
    
    obj<-cmest(data=subset(subset(dat, if (k=="PAD"){followupPAD>0} 
                                  else if (k=="AMI"){followupAMI>0}
                                  else if (k=="IS"){followupIS>0}
                                  else if (k=="ASCVD"){followupASCVD>0}
                                  else {sex!="hej"}), if (i=="hscrp_log"){hscrp>0} else {sex!="hej"}),
               model="rb", estimation = "paramfunc", inference = "delta", outcome=paste0("followup",k), event=k, 
               exposure = "apob_cut2", mediator=c(i), EMint=F, 
               basec= if (i=="VLDLP_log"| i=="LDLP"){c("startage","sex", "status_smoking", "cum_smoking", "birthdate", "systolic","diastolic")} else {
                 c("startage","sex", "status_smoking", "cum_smoking", "birthdate", "systolic","diastolic")},
               yreg="coxph", mreg=list("linear"),mval=list(1) , yval=TRUE, nboot=1, multimp = F)
    
    mediation_results[paste(i,k),c("mediation","lower","upper", "P.val")]<-as.data.frame(summary(obj)[c("effect.pe","effect.ci.low","effect.ci.high","effect.pval")])["pm", ]
    
    
    print(paste(i,k,"done!"))
    
  }
rm(dat)
mediation_results$variable<-factor(rownames(mediation_results),levels=rownames(mediation_results))
mediation_results

#Converting P-values of negative mediation
mediation_results[mediation_results[,"mediation"]<0,"P.val"]<-1-mediation_results[mediation_results[,"mediation"]<0,"P.val"]

#Minumum mediation 0, maximum mediation 1
for (i in c("mediation", "lower", "upper")){
  mediation_results[mediation_results[,i]<0,i]<-0
  mediation_results[mediation_results[,i]>1,i]<-1
}
mediation_results


#Creating tables for figures
fig5_table<- tibble(mediation_results$variable,mediation=paste(signif(mediation_results$mediation,2)),P.val=paste(signif(mediation_results$P.val,2))) %>% flextable()
fig5_table<- compose(fig5_table,j="mediation", value=as_paragraph(sub("^.0^","",round(mediation_results$mediation*100))," (", 
                                                                  sub("^.0^","",round(mediation_results$lower*100)),"-", 
                                                                  sub("^.0^","",round(mediation_results$upper*100)),")"))
fig5_table<-align(autofit(fig5_table, part="all"),align="left", part = "all")
fig5_table


group.colors <- rep(c("darkcyan","darkorchid4","darkcyan","darkorchid4"),3)

figure6cgps<-ggplot(data=mediation_results, aes(y = variable ,  x = mediation*100,  xmin=lower*100,xmax=upper*100, fill = variable ))+
  labs(y="Variable", x="Explained risk, %")+
  theme_nothing()+geom_col(width = 0.6, fill=group.colors, size=1)+geom_vline(xintercept=0, size=0.5, linetype=1)+
  theme(axis.ticks.length.x = unit(0.2, "cm"),axis.title = element_text(size=10), axis.text.x = element_text(vjust=0.1),axis.text.y = element_text(hjust=0.4, vjust=0.3), axis.title.y = element_text(vjust = 3, angle=90,size=10), axis.ticks=element_line(size=1), axis.text = element_text(size=10), axis.line.x.bottom =element_line(size=1),axis.line.y.left =element_blank() ,axis.line.x.top=element_blank() ,plot.margin=margin(0,3,0,3),plot.subtitle = element_blank(), plot.title = element_text(), 
        panel.grid=element_blank())+scale_x_continuous(expand=c(0,0))+scale_y_discrete(limits=rev)+coord_cartesian(xlim=c(0,100))



toffice(figure = figure6cgps, format = "pptx", title = "Copenhagen General Population Study",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,left = 0.15,top = 0.15,
        append = FALSE, width = 4, height = 3,devsize = FALSE, units = "in")




# Figure SX: aalen johansen --------------------------------------

fig4<-list()
for (k in c("PAD","CLTI","AMI"))
  for (i in c( "apob_mgdl","remncholcalc","ldl")){
    
    dat<-subset(cgps,startage>30)
  
    dat[,paste0(i,"_cut2")]<-cut2(dat[,i],if(i=="apob_mgdl"){c(100,150)} else if(i=="remncholcalc"){c(0.5,1)} else if (i=="ldl"){c(3,4)})
    
    fit<-survfit(as.formula(paste0("Surv(startage,","stopage",k,",",k,"_comprisk,type='mstate')~",i,"_cut2")),data = dat,id=obushnr)
    
    dat<-summary(fit)

    dat<-data.frame(group=factor(dat[["strata"]]),x=dat[["time"]],y=dat[["pstate"]][,3]*100)
    
   fig4[[paste0(k,i)]]<-ggplot(dat,aes(x=x,y=y,group=group,color=group))+
                                geom_step(size=1)+
                                scale_color_manual(values=c("seagreen","orange","red"))+
                                theme_classic()+theme(legend.position = "none")+labs(title=i,x="Age, years",y=paste0("Cumulative incidence of ",k,", %"))+
                                scale_y_continuous(limits=c(0,20),expand=c(0,0))+scale_x_continuous(limits=c(30,90),breaks=seq(30,90,10),expand=c(0,0))
                            
    
    
    
  }

#Changing function plot_grid
fargs  <- formals(plot_grid)
fargs$ncol           <- 3
formals(plot_grid) <- fargs

figure4<-do.call(plot_grid,fig4)
figure4

toffice(figure = figure4, format = "pptx",
        filename = "Temp.pptx", nr = 1, nc = 1, irow = 1, icol = 1,onsame = FALSE,
        append = TRUE, width = 11, height = 7, devsize = FALSE, units = "in")


