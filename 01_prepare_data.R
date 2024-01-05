library(tidyverse)
library(magrittr)
library(cocor)
library(Metrics)
library(multDM)
library(flextable)
library(officer)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggrepel)
library(gtools)
library(berryFunctions)

#######################
# MAKE FUNCTIONS
#######################

#Custom function for formatting small p and q values in SciNot
DigFormat<-function(d,cut,fixed,scientific){
  d<-as.numeric(d)
  if (is.na(d)){
    chard=NA
  }else if (abs(d) < cut){
    chard=sprintf(scientific,d)
  }else{
    chard=sprintf(fixed,d)
  }
  return(chard)
}

DigFormat_V<-Vectorize(DigFormat)

# cut= cutoff value after which numbers need scientific formatting
# fixed= format for when a number is above the cutoff value
# scientific = format for when a number is below the cutoff value
#
# Standard usage for NeRD Lab:
# DigFormat_V(VARNAME,cut=0.001,fixed="%.3f",scientific="%.2e")
# Usage for multiple variables in a data frame (requires tidyverse packages):
# across(c(VAR1,VAR2,VAR3), DigFormat_V,cut=0.001,fixed="%.3f",scientific="%.2e"))


#Function for running correlation tests and extracting needed values
runcor<-function(xdata,ydata,regions){
  rows<-length(regions)
  out.df<-as.data.frame(matrix(nrow=rows,ncol=6))
  colnames(out.df)<-c("brain_region","t","r","pval","CIlow","CIhigh")
  brain.reg<-paste(regions,sep=",")
  
  for (i in 1:length(brain.reg)){
    var=brain.reg[[i]]
    save<-cor.test(xdata[[var]],ydata[[var]])
    out.df[i,"t"]<-save$statistic
    out.df[i,"r"]<-save$estimate
    out.df[i,"pval"]<-save$p.value
    out.df[i,"CIlow"]<-save$conf.int[1]
    out.df[i,"CIhigh"]<-save$conf.int[2]
    out.df[i,"brain_region"]<-brain.reg[i]
  }
  return(out.df)
}

#Function for running correlation tests with demographic variables and extracting needed values
rundemcor<-function(xdata,ydata,demvar){
  regions=colnames(xdata[,-1])
  rows<-length(regions)
  out.df<-as.data.frame(matrix(nrow=rows,ncol=6))
  colnames(out.df)<-c("brain_region","t","r","pval","CIlow","CIhigh")
  brain.reg<-paste(regions,sep=",")
  
  for (i in 1:length(brain.reg)){
    var=brain.reg[[i]]
    save<-cor.test(xdata[[var]],ydata[[demvar]])
    out.df[i,"t"]<-save$statistic
    out.df[i,"r"]<-save$estimate
    out.df[i,"pval"]<-save$p.value
    out.df[i,"CIlow"]<-save$conf.int[1]
    out.df[i,"CIhigh"]<-save$conf.int[2]
    out.df[i,"brain_region"]<-brain.reg[i]
  }
  return(out.df)
}


#Function for calculating a Steiger's Z test
test2r.steigerz1 <- function (ry.x1, ry.x2, rx1.x2, n){ #included with thanks to Bruce Dudek
  fz1 <- 0.5 * log((1 + ry.x1)/(1 - ry.x1))           # https://rdrr.io/github/bcdudek/bcdstats/
  fz2 <- 0.5 * log((1 + ry.x2)/(1 - ry.x2))
  fz3 <- 0.5 * log((1 + rx1.x2)/(1 - rx1.x2))
  dif <- fz1-fz2
  av <- (ry.x1 + ry.x2)/2
  covnum1 <- rx1.x2 * (1-(av^2)-(av^2))
  covnum2 <- (.5*(av^2)) * (1-(av^2)-(av^2)-(rx1.x2^2))
  covdenom <- (1-(av^2))^2
  cov <- (covnum1-covnum2)/covdenom
  zteststat = (((n-3)^.5)*dif)/((2-(2*cov))^.5)
  p <- pnorm(abs(zteststat),0,1, lower.tail = FALSE)
  two_p <- p*2
  #if (twotailed)
  #    p <- 2 * p
  return(list(test = "test of difference between dependent correlations r(yx1) and r(yx2)",
              ry.x1=ry.x1, ry.x2=ry.x2,rx1.x2=rx1.x2,cov=cov,num1=covnum1,num2=covnum2,denom=covdenom,
              Difference_between_fishersz_correlations=dif,
              z.teststatistic = zteststat,onetail_p_value = p,twotail_p_value=two_p))
  
}


#Function for running Steiger's Z on a vectorized dataset and extracting needed values
runcocorvec<-function(regions,r.jk,r.jh,r.kh,n){
  rows=length(regions)
  out.df<-as.data.frame(matrix(nrow=rows,ncol=3))
  colnames(out.df)<-c("brain_region","steigers_Z","steigers_pval")
  brain.reg<-paste(regions,sep=",")
  
  for (i in 1:length(brain.reg)){
    var=i
    save<-test2r.steigerz1(ry.x1=r.jk[[var]],ry.x2=r.jh[[var]],rx1.x2=r.kh[[var]],n=n)
    out.df[i,"steigers_Z"]<-save$z.teststatistic
    out.df[i,"steigers_pval"]<-save$twotail_p_value
    out.df[i,"brain_region"]<-brain.reg[i]
  }  
  return(out.df)
}

#Function for identifying outliers
outlier.sd <- function(x) {
  return(x < mean(x)-(3*sd(x)) | x > mean(x)+(3*sd(x)))
}

#Function for making boxplots that compare individuals before & after synthsr processing to 3T
make_boxes<-function(df,x,y,ylab=y){
  df %>%
    group_by(get(x)) %>%
    mutate(outlier = ifelse(outlier.sd(get(y)), SubjID, NA)) %>%
    ggplot(aes(x=get(x),y=get(y))) +
    geom_point()+
    geom_line(aes(group=SubjID),linewidth=0.1,alpha=0.2)+
    geom_boxplot(alpha=0.4, width=0.5,aes(fill=get(x)))+
    theme_light()+
    ylab(ylab)+
    theme(axis.title.x=element_blank(),
          panel.grid = element_blank(),
          legend.position = 'none')
}

#Function for making a scaled difference score from a data frame
make_diff_dfs<-function(dfa,dfb,varlist){
  newdf<-tibble(SubjID=dfa$SubjID)
  for (i in 1:length(varlist)){
    varname<-varlist[[i]]
    newname<-paste0(varname)
    newdf<-newdf %>%
      mutate(newcol=as.vector(scale(abs(pull(dfa,varname)-pull(dfb,varname))))) %>%
      rename({{newname}}:=newcol)
  }
  return(newdf)
}


#######################
# READ IN DATA
#######################
#Get list of excluded participants and identify which aseg2table and aparc2table files to use
exclude<-as.character(readLines("exclusion_list.txt"))
date='20231016'

#Read in cortical thickness files
lh.CT<-read_delim(file=paste0(date,'_desikan.lh.CT.csv')) %>% select(-BrainSegVolNotVent)
rh.CT<-read_delim(file=paste0(date,'_desikan.rh.CT.csv')) %>% select(-BrainSegVolNotVent,-eTIV)
CT<-merge(lh.CT,rh.CT,by.x="lh.aparc.thickness",by.y="rh.aparc.thickness") %>% 
  select(-eTIV,-lh_MeanThickness_thickness, -rh_MeanThickness_thickness) %>%
  rename(SubjID=lh.aparc.thickness) %>%
  filter(!str_detect(SubjID,exclude)) 
ctregs<-names(CT)

#This is for WM surface area variables; we used pial surface area instead
# lh.SA<-read_delim(file=paste0(drive,date,'_polyneuro_desikan.lh.SA.csv')) %>% select(-BrainSegVolNotVent)
# rh.SA<-read_delim(file=paste0(drive,date,'_polyneuro_desikan.rh.SA.csv')) %>% select(-BrainSegVolNotVent,-eTIV)
# 
# SA<-merge(lh.SA,rh.SA,by.x="lh.aparc.area",by.y="rh.aparc.area") %>% 
#   select(-eTIV,-lh_WhiteSurfArea_area, -rh_WhiteSurfArea_area) %>%
#   rename(SubjID=lh.aparc.area) %>% 
#   filter(!str_detect(SubjID,exclude))
# saregs<-names(SA)

#Read in pial SA files
lh.pialSA<-read_delim(file=paste0(date,'_pial.lh.SA.csv')) %>% select(-BrainSegVolNotVent)
rh.pialSA<-read_delim(file=paste0(date,'_pial.rh.SA.csv')) %>% select(-BrainSegVolNotVent,-eTIV)

pialSA<-merge(lh.pialSA,rh.pialSA,by.x="lh.aparc.pial.area",by.y="rh.aparc.pial.area") %>% 
  select(-eTIV) %>%
  rename(SubjID=lh.aparc.pial.area) %>% 
  filter(!str_detect(SubjID,exclude)) %>%
  rename_with(~paste0(.,"pial"), -SubjID)
pialsaregs<-names(pialSA)

#Read in cortical volume files
lh.CV<-read_delim(file=paste0(date,'_desikan.lh.vol.csv')) %>% select(-BrainSegVolNotVent,-eTIV)
rh.CV<-read_delim(file=paste0(date,'_desikan.rh.vol.csv')) %>% select(-BrainSegVolNotVent,-eTIV)

CV<-merge(lh.CV,rh.CV,by.x="lh.aparc.volume",by.y="rh.aparc.volume") %>%
  rename(SubjID=lh.aparc.volume) %>% 
  filter(!str_detect(SubjID,exclude)) 
cvregs<-names(CV)

#Read in subcortical volume files and remove areas that aren't of interest
bl.sub<-read_delim(file=paste0(date,'_aseg.volume.csv'))
sub<-bl.sub%>% 
  select(-`5th-Ventricle`,-`Left-WM-hypointensities`,-`Right-WM-hypointensities`,
         -`Left-non-WM-hypointensities`,-`Right-non-WM-hypointensities`,
         -`WM-hypointensities`,-`non-WM-hypointensities`,-`BrainSegVol`,
         -`BrainSegVolNotVent`,-`lhCortexVol`,-`rhCortexVol`,-`CortexVol`,
         -`lhCerebralWhiteMatterVol`,-`rhCerebralWhiteMatterVol`,-`CerebralWhiteMatterVol`,
         -`SubCortGrayVol`,-`TotalGrayVol`,-`SupraTentorialVol`,-`SupraTentorialVolNotVent`,
         -`MaskVol`,-`BrainSegVol-to-eTIV`,-`MaskVol-to-eTIV`,
         -`lhSurfaceHoles`,-`rhSurfaceHoles`,-`SurfaceHoles`,-`EstimatedTotalIntraCranialVol`)
subregs<-names(sub)

sub<-sub %>% rename(SubjID="Measure:volume")  %>%
  filter(!str_detect(SubjID,exclude)) 


#Make a dataframe of all FreeSurfer regions of interest and a list of variables contained in it
RegionalVals<-full_join(CT, # SA) %>% full_join(
                        pialSA) %>% full_join(CV) %>% full_join(sub)
reglists<-list("CT"=ctregs[-1],#"SA"=saregs[-1],
               "pialSA"=pialsaregs[-1],"CV"=cvregs[-1],"sub"=subregs[-1])
allregs<-c(unlist(reglists))

#Make a list of global variables of interest
vars<-c("SubjID","Total_pial_area","Mean_thick","EstimatedTotalIntraCranialVol",
        "SubCortGrayVol","CortexVol","CerebralWhiteMatterVol","BrainSegVol")

#######################
# SELECT GLOBAL DATA
#######################

#Get only the global measures of interest into a separate dataframe
LHCT<-lh.CT %>% select(1,lh_MeanThickness_thickness) %>% rename("SubjID"=1)
RHCT<-rh.CT %>% select(1,rh_MeanThickness_thickness)%>% rename("SubjID"=1)

#No longer using WM surface area
# LHSA<-lh.SA %>% select(1,lh_WhiteSurfArea_area)%>% rename("SubjID"=1)
# RHSA<-rh.SA %>% select(1,rh_WhiteSurfArea_area)%>% rename("SubjID"=1)

#Total Pial surface area isn't part of the usually-extracted data so we had to pull it separately
RHPSA<-read_lines("extract_pialsurf_rh.txt") %>% str_replace(.,":",",") %>% 
  str_remove("/stats/rh.aparc.pial.stats")%>% as.tibble() %>% 
  separate(col=1,into=c("SubjID",NA,NA,NA,"rh.PialSurfArea",NA),sep=",") %>% mutate(rh.PialSurfArea=as.numeric(rh.PialSurfArea))
LHPSA<-read_lines("extract_pialsurf_lh.txt") %>% str_replace(.,":",",") %>% 
  str_remove("/stats/lh.aparc.pial.stats")%>% as.tibble() %>% 
  separate(col=1,into=c("SubjID",NA,NA,NA,"lh.PialSurfArea",NA),sep=",") %>% mutate(lh.PialSurfArea=as.numeric(lh.PialSurfArea))

SUBC<-bl.sub %>% select(1,EstimatedTotalIntraCranialVol,SubCortGrayVol,CortexVol,CerebralWhiteMatterVol,BrainSegVol) %>% 
  rename("SubjID"=1)

GlobalVals<-full_join(LHCT,RHCT) %>% 
  # full_join(LHSA) %>% full_join(RHSA) %>% 
  full_join(LHPSA) %>% full_join(RHPSA) %>%
  full_join(SUBC) %>%
  filter(!str_detect(SubjID,exclude)) %>% 
  mutate(#Total_area=lh_WhiteSurfArea_area+rh_WhiteSurfArea_area,
    Total_pial_area=lh.PialSurfArea+rh.PialSurfArea,
    Mean_thick=(lh_MeanThickness_thickness+rh_MeanThickness_thickness)/2)

remove(lh.CT,lh.SA,lh.CV,
       rh.CT,rh.SA,rh.CV,
       ctregs,saregs,cvregs,subregs,
       LHCT,LHSA,LHPSA,RHCT,RHSA,RHPSA,SUBC)
#######################
# GLOBAL PROCESSING
#######################
#Make empty tibble and identify the axial scan pattern
all_merged<-tibble(brain_region=character())
i="-axi[12]-1to1"

#Split the data into tibbles by scan type (e.g. 3T, SynthSR-processed, minimally-processed "raw") 
global.3Trun01<-GlobalVals %>%  filter(!str_detect(SubjID,exclude)) %>%
  filter(grepl('run-01',SubjID)) %>% separate(SubjID,into=c("sub","ses","acq"),sep="_") %>%
  mutate(SubjID=paste(sub,ses,sep="_")) %>% filter(rh_MeanThickness_thickness != 0) %>% 
  select(-sub,-ses)
global.synthsr<-GlobalVals %>% filter(!str_detect(SubjID,exclude)) %>%
  filter(grepl(paste0('synthsr',i,'$'),SubjID)) %>% separate(SubjID,into=c("sub","ses","acq"),sep="_") %>%
  mutate(SubjID=paste(sub,ses,sep="_")) %>% filter(rh_MeanThickness_thickness != 0) %>% 
  select(-sub,-ses) %>%
  group_by(SubjID) %>% slice_head(n=1) %>% ungroup() %>% mutate(acq=str_replace(acq,"axi[12]","axi"))
global.rawLF<-GlobalVals %>%  filter(!str_detect(SubjID,exclude)) %>%
  filter(grepl(paste0('rawLF',i,'$'),SubjID)) %>% separate(SubjID,into=c("sub","ses","acq"),sep="_") %>%
  mutate(SubjID=paste(sub,ses,sep="_")) %>% filter(rh_MeanThickness_thickness != 0) %>% 
  select(-sub,-ses)%>%
  group_by(SubjID) %>% slice_head(n=1) %>% ungroup()%>% mutate(acq=str_replace(acq,"axi[12]","axi"))
global.rawLF.all<-GlobalVals %>%  filter(!str_detect(SubjID,exclude)) %>%
  filter(grepl(paste0('rawLF','$'),SubjID)) %>% separate(SubjID,into=c("sub","ses","acq"),sep="_") %>%
  mutate(SubjID=paste(sub,ses,sep="_")) %>% filter(rh_MeanThickness_thickness != 0) %>% 
  select(-sub,-ses)%>%
  group_by(SubjID) %>% slice_head(n=1) %>% ungroup()
global.synthsr.all<-GlobalVals %>%  filter(!str_detect(SubjID,exclude)) %>%
  filter(grepl(paste0('synthsr','$'),SubjID)) %>% separate(SubjID,into=c("sub","ses","acq"),sep="_") %>%
  mutate(SubjID=paste(sub,ses,sep="_")) %>% filter(rh_MeanThickness_thickness != 0) %>% 
  select(-sub,-ses) %>%
  group_by(SubjID) %>% slice_head(n=1) %>% ungroup()

#Make sure all the participants are the same in each tibble so we don't get odd correlation results
global.3Trun01 %<>% filter(SubjID %in% global.synthsr$SubjID) %>%
  filter(SubjID %in% global.rawLF$SubjID)
global.synthsr %<>% filter(SubjID %in% global.3Trun01$SubjID) %>%
  filter(SubjID %in% global.rawLF$SubjID)
global.rawLF %<>% filter(SubjID %in% global.synthsr$SubjID) %>%
  filter(SubjID %in% global.3Trun01$SubjID)
global.rawLF.all%<>%filter(SubjID %in% global.rawLF$SubjID)%>%
  filter(SubjID %in% global.3Trun01$SubjID)
global.synthsr.all%<>%filter(SubjID %in% global.rawLF$SubjID)%>%
  filter(SubjID %in% global.3Trun01$SubjID)

#save the tibbles to a data file for easy access later
savelist=c("global.3Trun01","global.synthsr","global.rawLF","global.rawLF.all","global.synthsr.all","vars","GlobalVals")
save(list=savelist,file="global_data.Rdata")

#If all the filtering worked then save a list of participants for future reference
test<-nrow(global.3Trun01)==nrow(global.synthsr) & nrow(global.synthsr)==nrow(global.rawLF)
if (test) { 
  n=nrow(global.3Trun01)
  print(paste0("Total N = ",n))
  write(global.3Trun01$SubjID,file="sublist.txt")
}

##### RUN CORRELATIONS
#AIm 1 - axi1 raw vs axi1 synth
synth.3T.cor<-runcor(global.3Trun01,global.synthsr,regions=vars[-1]) %>%
  rename(synth3T_r=r,synth3T_pval=pval)

raw.3T.cor<-runcor(global.3Trun01,global.rawLF,regions=vars[-1]) %>%
  rename(raw3T_r=r,raw3T_pval=pval)

raw.synth.cor<-runcor(global.synthsr, global.rawLF,regions=vars[-1]) %>%
  rename(synthraw_r=r,synthraw_pval=pval)

#Aim 2 - axi1 raw vs all raw
rawall.3T.cor<-runcor(global.3Trun01,global.rawLF.all,regions=vars[-1]) %>%
  rename(rawall3T_r=r,rawall3T_pval=pval)

raw.rawall.cor<-runcor(global.rawLF,global.rawLF.all,regions=vars[-1]) %>%
  rename(rawallraw_r=r,rawallraw_pval=pval)

#Aim 3 - axi1 raw vs all synth
synthall.3T.cor<- runcor(global.3Trun01,global.synthsr.all,regions=vars[-1]) %>%
  rename(synthall3T_r=r,synthall3T_pval=pval)

raw.synthall.cor<-runcor(global.synthsr.all, global.rawLF,regions=vars[-1]) %>%
  rename(synthallraw_r=r,synthallraw_pval=pval)

#Aim 4 - axi1 synth vs all synth
synth.synthall.cor<-runcor(global.synthsr,global.synthsr.all,regions=vars[-1]) %>%
  rename(synthallsynth_r=r,synthallsynth_pval=pval)

#Merge all cortest results together
merged<-tibble(brain_region=vars[-1])

merged<-merged %>% left_join(select(synth.3T.cor,brain_region,synth3T_r,synth3T_pval)) %>%
  left_join(select(raw.3T.cor,brain_region,raw3T_r,raw3T_pval)) %>%
  left_join(select(raw.synth.cor,brain_region,synthraw_r,synthraw_pval)) %>% 
  left_join(select(rawall.3T.cor,brain_region,rawall3T_r,rawall3T_pval)) %>% 
  left_join(select(raw.rawall.cor,brain_region,rawallraw_r,rawallraw_pval)) %>% 
  left_join(select(synthall.3T.cor,brain_region,synthall3T_r,synthall3T_pval)) %>% 
  left_join(select(raw.synthall.cor,brain_region,synthallraw_r,synthallraw_pval)) %>%
  left_join(select(synth.synthall.cor,brain_region,synthallsynth_r,synthallsynth_pval))


##### RUN STEIGER TESTS
#Aim 1 tests
steigers_aim1<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$synth3T_r,r.jh=merged$raw3T_r,r.kh=merged$synthraw_r,n=n) %>%
  rename(aim1steigers_Z=steigers_Z,aim1steigers_pval=steigers_pval)

#Aim 2 tests
steigers_aim2<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$rawall3T_r,r.jh=merged$raw3T_r,r.kh=merged$rawallraw_r,n=n) %>%
  rename(aim2steigers_Z=steigers_Z,aim2steigers_pval=steigers_pval)

#Aim 3 tests
steigers_aim3<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$synthall3T_r,r.jh=merged$raw3T_r,r.kh=merged$synthallraw_r,n=n) %>%
  rename(aim3steigers_Z=steigers_Z,aim3steigers_pval=steigers_pval)

#Aim 4 tests
steigers_aim4<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$synthall3T_r,r.jh=merged$synth3T_r,r.kh=merged$synthallsynth_r,n=n) %>%
  rename(aim4steigers_Z=steigers_Z,aim4steigers_pval=steigers_pval)

#merge all together
merged<-left_join(merged,steigers_aim1) %>% 
  left_join(steigers_aim2) %>% 
  left_join(steigers_aim3) %>% 
  left_join(steigers_aim4) %>%
    select(-contains(c("synthraw","rawallraw","synthallraw")))


##### ADJUST FOR MULTIPLE COMPARISONS

#Calculate qvals for entire global dataset
qvals<-merged %>% select(brain_region,contains("pval")) %>% 
  pivot_longer(cols=contains("pval"),names_to=c("measure",NA), names_sep="_",values_to="pval") %>%
  mutate(qval=p.adjust(as.numeric(pval),"fdr")) %>%
  pivot_wider(names_from=measure,values_from=c(pval,qval),names_glue="{measure}_{.value}")

#Merge qvals into dataframe and order columns for easier reading
all_merged<-left_join(merged,qvals) %>% select(brain_region,
                                               raw3T_r,raw3T_pval,raw3T_qval,
                                               synth3T_r,synth3T_pval,synth3T_qval,
                                               aim1steigers_Z,aim1steigers_pval,aim1steigers_qval,
                                               rawall3T_r,rawall3T_pval,rawall3T_qval,
                                               aim2steigers_Z,aim2steigers_pval,aim2steigers_qval,
                                               synthall3T_r,synthall3T_pval,synthall3T_qval,
                                               aim3steigers_Z,aim3steigers_pval,aim3steigers_qval,
                                               aim4steigers_Z,aim4steigers_pval,aim4steigers_qval) %>%
  mutate(across(contains('val'), DigFormat_V,cut=0.001,fixed="%.3f",scientific="%.2e"),
        across(ends_with(c("_r","_Z")),DigFormat_V,cut=0.01,fixed="%.2f",scientific="%.1e"),
        brain_region=recode(brain_region,"Mean_thick"="Mean Cortical Thickness",
                            # "Total_area"="Total WM Surface Area", #Not using WM SA
                            "Total_pial_area"="Total Surface Area",
                             "EstimatedTotalIntraCranialVol"="Estimated Intracranial Volume",
                             "SubCortGrayVol"="Subcortical Gray Matter Volume",
                              "CortexVol" = "Cortical Volume",
                            "CerebralWhiteMatterVol" = "Cerebral White Matter Volume",
                            "BrainSegVol" = "Total Brain Volume"))

save(all_merged,file="global_all_merged.Rdat")

#######################
# GLOBAL TABLE GENERATION
#######################

#Fresh start point for generating tables
load("global_all_merged.Rdat")

#Prepare table without aim 4 (presented separately from other aims in manuscript)
flexprep<-all_merged %>% select(-contains("aim4"))%>%
  pivot_longer(cols=c(5:22),names_to=c("comp","meas"),names_sep="_",values_to="val") %>% 
  mutate(meas=ifelse(str_detect(comp,"steigers"),paste0("steigers",meas),meas),
         comp=case_when(comp=="synth3T" ~ "aim1",
                        comp=="rawall3T" ~ "aim2",
                        comp=="synthall3T" ~ "aim3",
                        TRUE ~ comp), 
         comp=str_remove(comp,"steigers")) %>% 
  pivot_wider(id_cols=c(1:5),names_from=meas,names_prefix="comparison_",values_from=val) %>%
  select(comp,everything())%>% arrange(comp)

#Make tibble with headers for each aim in global measures table because it is stacked
insertdf<-tibble(comp=c("aim1","aim1","aim2","aim2","aim3","aim3"),
                 brain_region=rep("Brain measure",times=6),
                 raw3T_r=rep(c("Standard Axial 64mT Correlations with 3T","r"),times=3),
                 raw3T_pval=rep(c("Standard Axial 64mT Correlations with 3T","p"),times=3),
                 raw3T_qval=rep(c("Standard Axial 64mT Correlations with 3T","q"),times=3),
                 comparison_r=c("SynthSR-Processed Axial 64mT Correlations with 3T","r",
                                "Standard Multi-Orientation 64mT Correlations with 3T","r",
                                "SynthSR-Processed Multi-Orientation 64mT Correlations with 3T","r"),
                 comparison_pval=c("SynthSR-Processed Axial 64mT Correlations with 3T","p",
                                "Standard Multi-Orientation 64mT Correlations with 3T","p",
                                "SynthSR-Processed Multi-Orientation 64mT Correlations with 3T","p"),
                 comparison_qval=c("SynthSR-Processed Axial 64mT Correlations with 3T","q",
                                "Standard Multi-Orientation 64mT Correlations with 3T","q",
                                "SynthSR-Processed Multi-Orientation 64mT Correlations with 3T","q"),
                 comparison_steigersZ=rep(c("Steiger's","Z"),times=3),
                 comparison_steigerspval=rep(c("Steiger's","p"),times=3),
                 comparison_steigersqval=rep(c("Steiger's","q"),times=3)
                 )

#Use insertRows function to add headers to subtables
flexprep<-flexprep %>% 
  insertRows(r=c(0,1,length(vars),length(vars),2*length(vars)-1,2*length(vars)-1),rcurrent=TRUE,new=insertdf)

#Use flextable to make presentable global vals table with formatted headers and save to html
fullflex<-flexprep %>% select(-1) %>% flextable () %>% 
  theme_booktabs() %>% delete_part("header") %>%
  vline(j=c(1,4,7),border=fp_border(color = "black",width=2)) %>%
  hline(i=c(2,9,11,18,20),border=fp_border(color = "black",width=2)) %>%
  merge_h(i=c(1,10,19),part="body") %>% merge_v(j=1,part="body") %>%
  border_outer(border=fp_border(color="black",width=2),part="body") %>%
  bold(i= ~ as.numeric(`raw3T_qval`)<0.05, j=2:4,part="body") %>%
  bold(i= ~ as.numeric(`comparison_qval`)<0.05, j=5:7,part="body") %>%
  bold(i= ~ as.numeric(`comparison_steigersqval`)<0.05, j=8:10,part="body") %>%
  bold(i=c(1,2,10,11,19,20),part="body") %>%
  italic(i=c(2,11,20),j=c(2:10)) %>%
  align(align="center",part="body") %>%
  fix_border_issues() %>% width(j=1,width=2.5,unit="in") %>%
  autofit()%>%
  add_header_lines(values="Table 2. Correspondence between low-field (64mT) and high-field (3T) MR images and comparisons across super-resolution approaches")
  
save_as_html(fullflex,path="../manuscript/Table2.html")


#Make separate flextable for Aim1
flexAim1<-select(all_merged,brain_region,
                 raw3T_r,raw3T_pval,raw3T_qval,
                 synth3T_r,synth3T_pval,synth3T_qval,
                 aim1steigers_Z,aim1steigers_pval,aim1steigers_qval) %>%
  flextable()  %>% theme_booktabs() %>% 
  vline(j=c(1,4,7),border=fp_border(color = "black",width=2)) %>% 
  border_outer(border=fp_border(color="black",width=2),part="all") %>%
  bold(i= ~ as.numeric(`raw3T_qval`)<0.05, j=2:4,part="body") %>%
  bold(i= ~ as.numeric(`synth3T_qval`)<0.05, j=5:7,part="body") %>%
  bold(i= ~ as.numeric(`aim1steigers_qval`)<0.05, j=8:10,part="body") %>%
  bold(part="header") %>%
  add_header_row(top=TRUE,values=c("Measurement",
                                   "Freesurfer 3T T1/T2 vs. Raw LF Axial Orientation T1/T2",
                                   "Freesurfer 3T T1/T2 vs. Fully Coregistered SynthSR T1/T2",
                                   "Steiger"),
                 colwidths = c(1,3,3,3)) %>% 
  set_header_labels(values=c("Measurement",
                             "r","p","q",
                             "r","p","q",
                             "z","p","q")) %>%
  italic(i=2,part="header") %>% merge_v(part="header") %>% merge_v(j=1,part="body") %>% 
  align(align="center",part="all") %>% 
  fix_border_issues() #%>% save_as_html(path="Aim1_Global_rawAxi_vs_synthAxi.html",title="Aim 1, Global Measures: Raw Axi T1/T2 vs SynthSR Axi T1/T2")

#Make separate flextable for Aim2
flexAim2<-select(all_merged,brain_region,
                 raw3T_r,raw3T_pval,raw3T_qval,
                 rawall3T_r,rawall3T_pval,rawall3T_qval,
                 aim2steigers_Z,aim2steigers_pval,aim2steigers_qval) %>%
  flextable()  %>% theme_booktabs() %>% 
  vline(j=c(1,4,7),border=fp_border(color = "black",width=2)) %>% 
  border_outer(border=fp_border(color="black",width=2),part="all") %>%
  bold(i= ~ as.numeric(`raw3T_qval`)<0.05, j=2:4,part="body") %>%
  bold(i= ~ as.numeric(`rawall3T_qval`)<0.05, j=5:7,part="body") %>%
  bold(i= ~ as.numeric(`aim2steigers_qval`)<0.05, j=8:10,part="body") %>%
  bold(part="header") %>%
  add_header_row(top=TRUE,values=c("Measurement",
                                   "Freesurfer 3T T1/T2 vs. Raw LF Axial Orientation T1/T2",
                                   "Freesurfer 3T T1/T2 vs. Raw LF All Scans T1/T2",
                                   "Steiger"),
                 colwidths = c(1,3,3,3)) %>% 
  set_header_labels(values=c("Measurement",
                             "r","p","q",
                             "r","p","q",
                             "z","p","q")) %>%
  italic(i=2,part="header") %>% merge_v(part="header") %>% merge_v(j=1,part="body") %>% 
  align(align="center",part="all") %>%
  fix_border_issues() # %>% save_as_html(path="Aim2_Global_rawAxi_vs_rawall.html",title="Aim 2, Global Measures: Raw Axi T1/T2 vs Raw All T1/T2")

#Make separate flextable for Aim3
flexAim3<-select(all_merged,brain_region,
                 raw3T_r,raw3T_pval,raw3T_qval,
                 synthall3T_r,synthall3T_pval,synthall3T_qval,
                 aim3steigers_Z,aim3steigers_pval,aim3steigers_qval) %>%
  flextable()  %>% theme_booktabs() %>% 
  vline(j=c(1,4,7),border=fp_border(color = "black",width=2)) %>% 
  border_outer(border=fp_border(color="black",width=2),part="all") %>%
  bold(i= ~ as.numeric(`raw3T_qval`)<0.05, j=2:4,part="body") %>%
  bold(i= ~ as.numeric(`synthall3T_qval`)<0.05, j=5:7,part="body") %>%
  bold(i= ~ as.numeric(`aim3steigers_qval`)<0.05, j=8:10,part="body") %>%
  bold(part="header") %>%
  add_header_row(top=TRUE,values=c("Measurement",
                                   "Freesurfer 3T T1/T2 vs. Raw LF Axial Orientation T1/T2",
                                   "Freesurfer 3T T1/T2 vs. SynthSR All Scans T1/T2",
                                   "Steiger"),
                 colwidths = c(1,3,3,3)) %>% 
  set_header_labels(values=c("Measurement",
                             "r","p","q",
                             "r","p","q",
                             "z","p","q")) %>%
  italic(i=2,part="header") %>% merge_v(part="header") %>% merge_v(j=1,part="body") %>% 
  align(align="center",part="all") %>% 
  fix_border_issues() # %>% save_as_html(path="Aim3_Global_rawAxi_vs_synthall.html",title="Aim 3, Global Measures: Raw Axi T1/T2 vs SynthSR All T1/T2")

#Make separate flextable for Aim4
flexAim4<-all_merged %>% select(brain_region,synth3T_r,synth3T_pval,synth3T_qval,
                                synthall3T_r,synthall3T_pval,synthall3T_qval,
                                aim4steigers_Z,aim4steigers_pval,aim4steigers_qval) %>%
  flextable()%>%theme_booktabs() %>%
  vline(j=c(1,4,7),border=fp_border(color = "black",width=2),part="all") %>%
  border_outer(border=fp_border(color="black",width=2),part="all") %>%
  bold(i= ~ as.numeric(`synth3T_qval`)<0.05, j=2:4,part="body") %>%
  bold(i= ~ as.numeric(`synthall3T_qval`)<0.05, j=5:7,part="body") %>%
  bold(i= ~ as.numeric(`aim4steigers_qval`)<0.05, j=8:10,part="body") %>%
  bold(part="header") %>%
  add_header_row(top=TRUE,values=c("Measurement",
                                   "Freesurfer 3T T1/T2 vs. SynthSR Axial Orientation T1/T2",
                                   "Freesurfer 3T T1/T2 vs. SynthSR All Scans T1/T2",
                                   "Steiger"),
                 colwidths = c(1,3,3,3)) %>% 
  set_header_labels(values=c("Measurement",
                             "r","p","q",
                             "r","p","q",
                             "z","p","q")) %>%
  italic(i=2,part="header") %>% merge_v(part="header") %>% merge_v(j=1,part="body") %>% 
  align(align="center",part="all") %>%
  fix_border_issues() 

save_as_html(flexAim4, path="Global_synthAxi_vs_synthall.html",title="Global Measures: SynthSR Axi T1/T2 vs SynthSR All T1/T2")

#######################
# Global FIGURE GENERATION
#######################

#Make readable labels for global vars
varlabs=c(NA,recode(all_merged$brain_region,"Mean_thick"="Mean Cortical Thickness",
                    # "Total_area"="Total WM Surface Area",
                    "Total_pial_area"="Total Surface Area",
                    "EstimatedTotalIntraCranialVol"="Estimated Intracranial Volume",
                    "SubCortGrayVol"="Subcortical Gray Matter Volume",
                    "CortexVol" = "Cortical Volume",
                    "CerebralWhiteMatterVol" = "Cerebral White Matter Volume",
                    "BrainSegVol" = "Total Brain Volume"))

####PLOT AIM 1
#Make a tibble of data needed for Aim1 
aim1df<-full_join(global.3Trun01,global.rawLF) %>% full_join(global.synthsr) %>% 
  mutate(acq=case_when(acq=='run-01'~'3T',
                       acq=='rawLF-axi-1to1' ~ "0.64mT",
                       acq=='synthsr-axi-1to1' ~ "SynthSR 0.64mT"),
         acq=factor(acq,levels=c('0.64mT','SynthSR','3T')))

#Loop through variables and make a comparison boxplot for each
plots<-list()
for (i in 2:length(vars)){
  var=vars[[i]]
  varlab=varlabs[[i]]
  plot<-make_boxes(aim1df,"acq",var,ylab=varlab)
  if(i!=3){
    plot=plot+scale_y_continuous(labels = function(x) format(x, scientific = TRUE,digits=2))
  }
  plots[[i-1]]<-plot
}

#Arrange boxplots using patchwork and save output
pattern<-"AABBCCDD#\n#EEFFGG##"
allAim1<-wrap_elements(plots[[1]])+wrap_elements(plots[[2]])+wrap_elements(plots[[3]])+wrap_elements(plots[[4]])+
  wrap_elements(plots[[5]])+wrap_elements(plots[[6]])+wrap_elements(plots[[7]])+plot_layout(design=pattern)
ggsave("Aim1_GlobalVals_Boxplots.pdf",plot=allAim1,width=10,height=5)

####PLOT AIM 2
#Make a tibble of data needed for Aim 2

aim2df<-full_join(global.3Trun01,global.rawLF) %>% full_join(global.rawLF.all) %>% 
  mutate(acq=case_when(acq=='run-01'~'3T',
                       acq=='rawLF-axi-1to1' ~ "0.64mT",
                       acq=='rawLF' ~ "Repeated 0.64mT"),
         acq=factor(acq,levels=c('0.64mT','Repeated 0.64mT','3T')))

#Loop through variables and make a comparison boxplot for each
plots<-list()
for (i in 2:length(vars)){
  var=vars[[i]]
  varlab=varlabs[[i]]
  plot<-make_boxes(aim2df,"acq",var,ylab=varlab)
  if(i!=3){
    plot=plot+scale_y_continuous(labels = function(x) format(x, scientific = TRUE,digits=2))
  }
  plots[[i-1]]<-plot
}

#Arrange boxplots using patchwork and save output
allAim2<-wrap_elements(plots[[1]])+wrap_elements(plots[[2]])+wrap_elements(plots[[3]])+wrap_elements(plots[[4]])+
  wrap_elements(plots[[5]])+wrap_elements(plots[[6]])+wrap_elements(plots[[7]])+plot_layout(design=pattern)
ggsave("Aim2_GlobalVals_Boxplots.pdf",plot=allAim2,width=10,height=5)

####PLOT AIM 3
#Make a tibble of data needed for Aim 3
aim3df<-full_join(global.3Trun01,global.rawLF) %>% full_join(global.synthsr.all) %>% 
  mutate(acq=case_when(acq=='run-01'~'3T',
                       acq=='rawLF-axi-1to1' ~ "0.64mT",
                       acq=='SynthSR' ~ "Repeated + SynthSR 0.64mT"),
         acq=factor(acq,levels=c('0.64mT','Repeated + SynthSR 0.64mT','3T')))

#Loop through variables and make a comparison boxplot for each
plots<-list()
for (i in 2:length(vars)){
  var=vars[[i]]
  varlab=varlabs[[i]]
  plot<-make_boxes(aim3df,"acq",var,ylab=varlab)
  if(i!=3){
    plot=plot+scale_y_continuous(labels = function(x) format(x, scientific = TRUE,digits=2))
  }
  plots[[i-1]]<-plot
}

#Arrange boxplots using patchwork and save output
allAim3<-wrap_elements(plots[[1]])+wrap_elements(plots[[2]])+wrap_elements(plots[[3]])+wrap_elements(plots[[4]])+
  wrap_elements(plots[[5]])+wrap_elements(plots[[6]])+wrap_elements(plots[[7]])+plot_layout(design=pattern)
ggsave("Aim3_GlobalVals_Boxplots.pdf",plot=allAim3,width=10,height=5)

#######################
# REGIONAL PROCESSING
#######################

#Make empty tibble and identify the axial scan pattern

all_merged<-tibble(brain_region=character())
i="-axi[12]-1to1"

#Split the data into tibbles by scan type (e.g. 3T, SynthSR-processed, minimally-processed "raw") 
regional.3Trun01<-RegionalVals %>%  filter(!str_detect(SubjID,exclude)) %>%
  filter(grepl('run-01',SubjID)) %>% separate(SubjID,into=c("sub","ses","acq"),sep="_") %>%
  mutate(SubjID=paste(sub,ses,sep="_")) %>% select(-sub,-ses)
regional.synthsr<-RegionalVals %>% filter(!str_detect(SubjID,exclude)) %>%
  filter(grepl(paste0('synthsr',i,'$'),SubjID)) %>% separate(SubjID,into=c("sub","ses","acq"),sep="_") %>%
  mutate(SubjID=paste(sub,ses,sep="_")) %>% select(-sub,-ses) %>%
  group_by(SubjID) %>% slice_head(n=1) %>% ungroup() %>% mutate(acq=str_replace(acq,"axi[12]","axi"))
regional.rawLF<-RegionalVals %>%  filter(!str_detect(SubjID,exclude)) %>%
  filter(grepl(paste0('rawLF',i,'$'),SubjID)) %>% separate(SubjID,into=c("sub","ses","acq"),sep="_") %>%
  mutate(SubjID=paste(sub,ses,sep="_")) %>% select(-sub,-ses)%>%
  group_by(SubjID) %>% slice_head(n=1) %>% ungroup()%>% mutate(acq=str_replace(acq,"axi[12]","axi"))
regional.rawLF.all<-RegionalVals %>%  filter(!str_detect(SubjID,exclude)) %>%
  filter(grepl(paste0('rawLF','$'),SubjID)) %>% separate(SubjID,into=c("sub","ses","acq"),sep="_") %>%
  mutate(SubjID=paste(sub,ses,sep="_")) %>% select(-sub,-ses)%>%
  group_by(SubjID) %>% slice_head(n=1) %>% ungroup()
regional.synthsr.all<-RegionalVals %>%  filter(!str_detect(SubjID,exclude)) %>%
  filter(grepl(paste0('synthsr','$'),SubjID)) %>% separate(SubjID,into=c("sub","ses","acq"),sep="_") %>%
  mutate(SubjID=paste(sub,ses,sep="_")) %>% select(-sub,-ses) %>%
  group_by(SubjID) %>% slice_head(n=1) %>% ungroup()

#Make sure all the participants are the same in each tibble so we don't get odd correlation results
regional.3Trun01 %<>% filter(SubjID %in% regional.synthsr$SubjID) %>%
  filter(SubjID %in% regional.rawLF$SubjID)
regional.synthsr %<>% filter(SubjID %in% regional.3Trun01$SubjID) %>%
  filter(SubjID %in% regional.rawLF$SubjID)
regional.rawLF %<>% filter(SubjID %in% regional.synthsr$SubjID) %>%
  filter(SubjID %in% regional.3Trun01$SubjID)
regional.rawLF.all%<>%filter(SubjID %in% regional.rawLF$SubjID)%>%
  filter(SubjID %in% regional.3Trun01$SubjID)
regional.synthsr.all%<>%filter(SubjID %in% regional.rawLF$SubjID)%>%
  filter(SubjID %in% regional.3Trun01$SubjID)

#save the tibbles to a data file for easy access later
savelist=c("regional.3Trun01","regional.synthsr","regional.rawLF","regional.rawLF.all","regional.synthsr.all","allregs")
save(list=savelist,file="regional_data.Rdata")

#If all the filtering worked then save a list of participants for future reference
test<-nrow(regional.3Trun01)==nrow(regional.synthsr) & nrow(regional.synthsr)==nrow(regional.rawLF)
if (test) { 
  n=nrow(regional.3Trun01)
  print(paste0("Total N = ",n))
}

##### RUN CORRELATIONS
#AIm 1 - axi1 raw vs axi1 synth
synth.3T.cor<-runcor(regional.3Trun01,regional.synthsr,regions=allregs) %>%
  rename(synth3T_r=r,synth3T_pval=pval)

raw.3T.cor<-runcor(regional.3Trun01,regional.rawLF,regions=allregs) %>%
  rename(raw3T_r=r,raw3T_pval=pval)

raw.synth.cor<-runcor(regional.synthsr, regional.rawLF,regions=allregs) %>%
  rename(synthraw_r=r,synthraw_pval=pval)

#Aim 2 - axi1 raw vs all raw
rawall.3T.cor<-runcor(regional.3Trun01,regional.rawLF.all,regions=allregs) %>%
  rename(rawall3T_r=r,rawall3T_pval=pval)

raw.rawall.cor<-runcor(regional.rawLF,regional.rawLF.all,regions=allregs) %>%
  rename(rawallraw_r=r,rawallraw_pval=pval)

#Aim 3 - axi1 raw vs all synth
synthall.3T.cor<- runcor(regional.3Trun01,regional.synthsr.all,regions=allregs) %>%
  rename(synthall3T_r=r,synthall3T_pval=pval)

raw.synthall.cor<-runcor(regional.synthsr.all, regional.rawLF,regions=allregs) %>%
  rename(synthallraw_r=r,synthallraw_pval=pval)

#Aim 4 - axi1 synth vs all synth
synth.synthall.cor<-runcor(regional.synthsr,regional.synthsr.all,regions=allregs) %>%
  rename(synthallsynth_r=r,synthallsynth_pval=pval)

#Merge all cortest results together
merged<-tibble(brain_region=allregs)

merged<-merged %>% left_join(select(synth.3T.cor,brain_region,synth3T_r,synth3T_pval)) %>%
  left_join(select(raw.3T.cor,brain_region,raw3T_r,raw3T_pval)) %>%
  left_join(select(raw.synth.cor,brain_region,synthraw_r,synthraw_pval)) %>% 
  left_join(select(rawall.3T.cor,brain_region,rawall3T_r,rawall3T_pval)) %>% 
  left_join(select(raw.rawall.cor,brain_region,rawallraw_r,rawallraw_pval)) %>% 
  left_join(select(synthall.3T.cor,brain_region,synthall3T_r,synthall3T_pval)) %>% 
  left_join(select(raw.synthall.cor,brain_region,synthallraw_r,synthallraw_pval)) %>%
  left_join(select(synth.synthall.cor,brain_region,synthallsynth_r,synthallsynth_pval))

##### RUN STEIGER TESTS
#Aim 1 tests
steigers_aim1<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$synth3T_r,r.jh=merged$raw3T_r,r.kh=merged$synthraw_r,n=n) %>%
  rename(aim1steigers_Z=steigers_Z,aim1steigers_pval=steigers_pval)

#Aim 2 tests
steigers_aim2<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$rawall3T_r,r.jh=merged$raw3T_r,r.kh=merged$rawallraw_r,n=n) %>%
  rename(aim2steigers_Z=steigers_Z,aim2steigers_pval=steigers_pval)

#Aim 3 tests
steigers_aim3<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$synthall3T_r,r.jh=merged$raw3T_r,r.kh=merged$synthallraw_r,n=n) %>%
  rename(aim3steigers_Z=steigers_Z,aim3steigers_pval=steigers_pval)

#Aim 4 tests
steigers_aim4<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$synthall3T_r,r.jh=merged$synth3T_r,r.kh=merged$synthallsynth_r,n=n) %>%
  rename(aim4steigers_Z=steigers_Z,aim4steigers_pval=steigers_pval)

#merge all together
merged<-left_join(merged,steigers_aim1) %>% 
  left_join(steigers_aim2) %>% 
  left_join(steigers_aim3) %>% 
  left_join(steigers_aim4) %>%   
  select(-contains(c("synthraw","rawallraw","synthallraw")))


##### ADJUST FOR MULTIPLE COMPARISONS

#Calculate qvals for entire regional dataset
qvals<-merged %>% select(brain_region,contains("pval")) %>% 
  pivot_longer(cols=contains("pval"),names_to=c("measure",NA), names_sep="_",values_to="pval") %>%
  mutate(qval=p.adjust(as.numeric(pval),"fdr")) %>%
  pivot_wider(names_from=measure,values_from=c(pval,qval),names_glue="{measure}_{.value}")

#Merge qvals into dataframe and order columns for easier reading
all_merged<-left_join(merged,qvals) %>% select(brain_region,
                                               raw3T_r,raw3T_pval,raw3T_qval,
                                               synth3T_r,synth3T_pval,synth3T_qval,
                                               aim1steigers_Z,aim1steigers_pval,aim1steigers_qval,
                                               rawall3T_r,rawall3T_pval,rawall3T_qval,
                                               aim2steigers_Z,aim2steigers_pval,aim2steigers_qval,
                                               synthall3T_r,synthall3T_pval,synthall3T_qval,
                                               aim3steigers_Z,aim3steigers_pval,aim3steigers_qval,
                                               aim4steigers_Z,aim4steigers_pval,aim4steigers_qval) %>%
  mutate(across(contains('val'), DigFormat_V,cut=0.001,fixed="%.3f",scientific="%.2e"),
         across(contains(c("_r$","_Z$")),\(x) round(x, digits=2))) 

save(all_merged,file="regional_all_merged.Rdat")

#######################
# REGIONAL TABLE GENERATION
#######################
#Too many variables for a combined table

#Make separate flextable for Aim1
flexAim1<-select(all_merged,brain_region,
                 raw3T_r,raw3T_pval,raw3T_qval,
                 synth3T_r,synth3T_pval,synth3T_qval,
                 aim1steigers_Z,aim1steigers_pval,aim1steigers_qval) %>%
  flextable()  %>% theme_booktabs() %>% 
  vline(j=c(1,4,7),border=fp_border(color = "black",width=2)) %>% 
  hline(i=seq(from=68,to=242,by=68),border=fp_border(color = "black",width=2)) %>%
  border_outer(border=fp_border(color="black",width=2),part="all") %>%
  bold(i= ~ as.numeric(`raw3T_qval`)<0.05, j=2:4,part="body") %>%
  bold(i= ~ as.numeric(`synth3T_qval`)<0.05, j=5:7,part="body") %>%
  bold(i= ~ as.numeric(`aim1steigers_qval`)<0.05, j=8:10,part="body") %>%
  bold(part="header") %>%
  add_header_row(top=TRUE,values=c("Measurement",
                                   "Freesurfer 3T T1/T2 vs. Raw LF Axial Orientation T1/T2",
                                   "Freesurfer 3T T1/T2 vs. Fully Coregistered SynthSR T1/T2",
                                   "Steiger"),
                 colwidths = c(1,3,3,3)) %>% 
  set_header_labels(values=c("Measurement",
                             "r","p","q",
                             "r","p","q",
                             "z","p","q")) %>%
  italic(i=2,part="header") %>% merge_v(part="header") %>% merge_v(j=1,part="body") %>% 
  align(align="center",part="all") %>%
  fix_border_issues() %>% save_as_html(path="Aim1_Regional_rawAxi_vs_synthAxi.html",title="Aim 1, Regional Measures: Raw Axi T1/T2 vs SynthSR Axi T1/T2")

#Make separate flextable for Aim2
flexAim2<-select(all_merged,brain_region,
                 raw3T_r,raw3T_pval,raw3T_qval,
                 rawall3T_r,rawall3T_pval,rawall3T_qval,
                 aim2steigers_Z,aim2steigers_pval,aim2steigers_qval) %>%
  flextable()  %>% theme_booktabs() %>% 
  vline(j=c(1,4,7),border=fp_border(color = "black",width=2)) %>% 
  hline(i=seq(from=68,to=242,by=68),border=fp_border(color = "black",width=2)) %>%
  border_outer(border=fp_border(color="black",width=2),part="all") %>%
  bold(i= ~ as.numeric(`raw3T_qval`)<0.05, j=2:4,part="body") %>%
  bold(i= ~ as.numeric(`rawall3T_qval`)<0.05, j=5:7,part="body") %>%
  bold(i= ~ as.numeric(`aim2steigers_qval`)<0.05, j=8:10,part="body") %>%
  bold(part="header") %>%
  add_header_row(top=TRUE,values=c("Measurement",
                                   "Freesurfer 3T T1/T2 vs. Raw LF Axial Orientation T1/T2",
                                   "Freesurfer 3T T1/T2 vs. Raw LF All Scans T1/T2",
                                   "Steiger"),
                 colwidths = c(1,3,3,3)) %>% 
  set_header_labels(values=c("Measurement",
                             "r","p","q",
                             "r","p","q",
                             "z","p","q")) %>%
  italic(i=2,part="header") %>% merge_v(part="header") %>% merge_v(j=1,part="body") %>% 
  align(align="center",part="all") %>% 
  fix_border_issues() %>% save_as_html(path="Aim2_Regional_rawAxi_vs_rawall.html",title="Aim 2, Regional Measures: Raw Axi T1/T2 vs Raw All T1/T2")

#Make separate flextable for Aim3
flexAim3<-select(all_merged,brain_region,
                 raw3T_r,raw3T_pval,raw3T_qval,
                 synthall3T_r,synthall3T_pval,synthall3T_qval,
                 aim3steigers_Z,aim3steigers_pval,aim3steigers_qval) %>%
  flextable()  %>% theme_booktabs() %>% 
  vline(j=c(1,4,7),border=fp_border(color = "black",width=2)) %>% 
  hline(i=seq(from=68,to=242,by=68),border=fp_border(color = "black",width=2)) %>%
  border_outer(border=fp_border(color="black",width=2),part="all") %>%
  bold(i= ~ as.numeric(`raw3T_qval`)<0.05, j=2:4,part="body") %>%
  bold(i= ~ as.numeric(`synthall3T_qval`)<0.05, j=5:7,part="body") %>%
  bold(i= ~ as.numeric(`aim3steigers_qval`)<0.05, j=8:10,part="body") %>%
  bold(part="header") %>%
  add_header_row(top=TRUE,values=c("Measurement",
                                   "Freesurfer 3T T1/T2 vs. Raw LF Axial Orientation T1/T2",
                                   "Freesurfer 3T T1/T2 vs. SynthSR All Scans T1/T2",
                                   "Steiger"),
                 colwidths = c(1,3,3,3)) %>% 
  set_header_labels(values=c("Measurement",
                             "r","p","q",
                             "r","p","q",
                             "z","p","q")) %>%
  italic(i=2,part="header") %>% merge_v(part="header") %>% merge_v(j=1,part="body") %>% 
  align(align="center",part="all") %>% 
  fix_border_issues() %>% save_as_html(path="Aim3_Regional_rawAxi_vs_synthall.html",title="Aim 3, Regional Measures: Raw Axi T1/T2 vs SynthSR All T1/T2")

#Make separate flextable for Aim4
flexAim4<-select(all_merged,brain_region,
                 synth3T_r,synth3T_pval,synth3T_qval,
                 synthall3T_r,synthall3T_pval,synthall3T_qval,
                 aim4steigers_Z,aim4steigers_pval,aim4steigers_qval) %>%
  flextable()  %>% theme_booktabs() %>% 
  vline(j=c(1,4,7),border=fp_border(color = "black",width=2)) %>% 
  hline(i=seq(from=68,to=242,by=68),border=fp_border(color = "black",width=2)) %>%
  border_outer(border=fp_border(color="black",width=2),part="all") %>%
  bold(i= ~ as.numeric(`synth3T_qval`)<0.05, j=2:4,part="body") %>%
  bold(i= ~ as.numeric(`synthall3T_qval`)<0.05, j=5:7,part="body") %>%
  bold(i= ~ as.numeric(`aim4steigers_qval`)<0.05, j=8:10,part="body") %>%
  bold(part="header") %>%
  add_header_row(top=TRUE,values=c("Measurement",
                                   "Freesurfer 3T T1/T2 vs. SynthSR Axial Orientation T1/T2",
                                   "Freesurfer 3T T1/T2 vs. SynthSR All Scans T1/T2",
                                   "Steiger"),
                 colwidths = c(1,3,3,3)) %>% 
  set_header_labels(values=c("Measurement",
                             "r","p","q",
                             "r","p","q",
                             "z","p","q")) %>%
  italic(i=2,part="header") %>% merge_v(part="header") %>% merge_v(j=1,part="body") %>% 
  align(align="center",part="all") %>% 
  fix_border_issues() %>% save_as_html(path="Regional_synthAxi_vs_synthAll.html",title="Regional Measures: SynthSR Axi T1/T2 vs SynthSR Repeated Scans T1/T2")



#############
#FOLLOW-UP TABLES
#############

####PULL MISCELLANEOUS INFORMATION
#Get mean and median of correlations for minimally processed single pair scans
get_reg_cors<-all_merged %>% separate(brain_region,into=c(NA,NA,"meas"),sep="_",remove=FALSE) %>%
  mutate(meas=ifelse(meas=="Anterior" | meas=="Posterior",NA,meas),
         meas=ifelse(is.na(meas),"sub",meas)) %>%
  group_by(meas) %>%
  summarize(rawaxirmean=mean(raw3T_r),
            rawaxirmedian=median(raw3T_r))

#Get average measures for 3T scan regions
Reg_summary<-summarise(regional.3Trun01,across(where(is.numeric),mean,names="mean_{.col}")) %>% 
  t() %>% 
  cbind(row.names(.)) %>%
  as.tibble(select(V2,V1)) %>%
  rename("brain_region"=V2,"avg_meas"=V1)

####CHECK CORRELATIONS WITH DEMOGRAPHICS AND MOTION

#Read in global data, motion info, and demographic info
load(file="global_data.Rdata")
motion<-read_csv("MotionAverages.csv") %>% filter(!str_detect(subj,"ses-02"))
demo<-read_delim("NIH_individual_level_data.csv") %>% mutate(subj=paste0("sub-",studyid,"_ses-01"))

#Get number of participants
n=nrow(global.3Trun01)

#Merge motion and demographic info together for easy access
motion_and_demo<-full_join(motion,demo) %>%
  filter(subj %in% global.3Trun01$SubjID)

#Get tibbles of the scaled difference score between 3T scans and lowfield scans
rawLF_3T_stdiff<-make_diff_dfs(dfa=global.3Trun01,dfb=global.rawLF,varlist=vars[-1])
rawLF.all_3T_stdiff<-make_diff_dfs(dfa=global.3Trun01,dfb=global.rawLF.all,varlist=vars[-1])
synthsr_3T_stdiff<-make_diff_dfs(dfa=global.3Trun01,dfb=global.synthsr,varlist=vars[-1])
synthsr.all_3T_stdiff<-make_diff_dfs(dfa=global.3Trun01,dfb=global.synthsr.all,varlist=vars[-1])

#Test relationship between motion during LF scans and the scaled 3T-LF difference score
rawLFdiff.motcor<-rundemcor(xdata=rawLF_3T_stdiff,ydata=motion_and_demo,demvar="fd") %>%
  rename(rawLFdiff.motcor_r=r,rawLFdiff.motcor_pval=pval)
rawLFalldiff.motcor<-rundemcor(xdata=rawLF.all_3T_stdiff,ydata=motion_and_demo,demvar="fd") %>%
  rename(rawLFalldiff.motcor_r=r,rawLFalldiff.motcor_pval=pval)
synthsrdiff.motcor<-rundemcor(xdata=synthsr_3T_stdiff,ydata=motion_and_demo,demvar="fd") %>%
  rename(synthsrdiff.motcor_r=r,synthsrdiff.motcor_pval=pval)
synthsralldiff.motcor<-rundemcor(xdata=synthsr.all_3T_stdiff,ydata=motion_and_demo,demvar="fd") %>%
  rename(synthsralldiff.motcor_r=r,synthsralldiff.motcor_pval=pval)

#Test relationship between participant age age and the scaled 3T-LF difference score
rawLFdiff.agecor<-rundemcor(xdata=rawLF_3T_stdiff,ydata=motion_and_demo,demvar="Age") %>%
  rename(rawLFdiff.agecor_r=r,rawLFdiff.agecor_pval=pval)
rawLFalldiff.agecor<-rundemcor(xdata=rawLF.all_3T_stdiff,ydata=motion_and_demo,demvar="Age") %>%
  rename(rawLFalldiff.agecor_r=r,rawLFalldiff.agecor_pval=pval)
synthsrdiff.agecor<-rundemcor(xdata=synthsr_3T_stdiff,ydata=motion_and_demo,demvar="Age") %>%
  rename(synthsrdiff.agecor_r=r,synthsrdiff.agecor_pval=pval)
synthsralldiff.agecor<-rundemcor(xdata=synthsr.all_3T_stdiff,ydata=motion_and_demo,demvar="Age") %>%
  rename(synthsralldiff.agecor_r=r,synthsralldiff.agecor_pval=pval)

#Test relationship between scaled 3T-LF difference score with different LF processing
rawLFdiff.synthsrdiff.cor<-runcor(xdata=rawLF_3T_stdiff,ydata=synthsr_3T_stdiff,regions=vars[-1]) %>%
  rename(rawLF_synthsr_r=r,rawLF_synthsr_pval=pval)
rawLFdiff.rawLFall.cor<-runcor(xdata=rawLF_3T_stdiff,ydata=rawLF.all_3T_stdiff,regions=vars[-1]) %>%
  rename(rawLF_rawLFall_r=r,rawLF_rawLFall_pval=pval)
rawLFdiff.synthsrall.cor<-runcor(xdata=rawLF_3T_stdiff,ydata=synthsr.all_3T_stdiff,regions=vars[-1]) %>%
  rename(rawLF_synthsrall_r=r,rawLF_synthsrall_pval=pval)

#Make empty table for reporting scaled difference score results
merged<-tibble(brain_region=vars[-1])

#Combine tibbles from preceding correlations
merged<-merged %>%
  left_join(select(rawLFdiff.motcor,brain_region,rawLFdiff.motcor_r,rawLFdiff.motcor_pval)) %>%
  left_join(select(synthsrdiff.motcor,brain_region,synthsrdiff.motcor_r,synthsrdiff.motcor_pval)) %>%
  left_join(select(rawLFdiff.synthsrdiff.cor,brain_region,rawLF_synthsr_r,rawLF_synthsr_pval))%>%
  
  left_join(select(rawLFalldiff.motcor,brain_region,rawLFalldiff.motcor_r,rawLFalldiff.motcor_pval)) %>%
  left_join(select(rawLFdiff.rawLFall.cor,brain_region,rawLF_rawLFall_r,rawLF_rawLFall_pval))%>%
  
  left_join(select(synthsralldiff.motcor,brain_region,synthsralldiff.motcor_r,synthsralldiff.motcor_pval)) %>%
  left_join(select(rawLFdiff.synthsrall.cor,brain_region,rawLF_synthsrall_r,rawLF_synthsrall_pval)) %>%
  
  left_join(select(rawLFdiff.agecor,brain_region,rawLFdiff.agecor_r,rawLFdiff.agecor_pval)) %>%
  left_join(select(synthsrdiff.agecor,brain_region,synthsrdiff.agecor_r,synthsrdiff.agecor_pval)) %>%
  left_join(select(rawLFalldiff.agecor,brain_region,rawLFalldiff.agecor_r,rawLFalldiff.agecor_pval)) %>%
  left_join(select(synthsralldiff.agecor,brain_region,synthsralldiff.agecor_r,synthsralldiff.agecor_pval)) 

#Did motion affect different LF processing differently?
steigers_synth.v.raw.motion<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$synthsrdiff.motcor_r,
                           r.jh=merged$rawLFdiff.motcor_r,
                           r.kh=merged$rawLF_synthsr_r,n=n) %>%
  rename(synthvraw.mot.steigers_Z=steigers_Z,synthvraw.mot.steigers_pval=steigers_pval)

steigers_rawall.v.raw.motion<-runcocorvec(regions=merged$brain_region,
                                          r.jk=merged$rawLFalldiff.motcor_r,
                                          r.jh=merged$rawLFdiff.motcor_r,
                                          r.kh=merged$rawLF_rawLFall_r,n=n)%>%
  rename(rawallvraw.mot.steigers_Z=steigers_Z,rawallvraw.mot.steigers_pval=steigers_pval)

steigers_synthall.v.raw.motion<-runcocorvec(regions=merged$brain_region,
                                          r.jk=merged$synthsralldiff.motcor_r,
                                          r.jh=merged$rawLFdiff.motcor_r,
                                          r.kh=merged$rawLF_synthsrall_r,n=n)%>%
  rename(synthallvraw.mot.steigers_Z=steigers_Z,synthallvraw.mot.steigers_pval=steigers_pval)

steigers_synth.v.raw.age<-runcocorvec(regions=merged$brain_region,
                                         r.jk=merged$synthsrdiff.agecor_r,
                                         r.jh=merged$rawLFdiff.agecor_r,
                                         r.kh=merged$rawLF_synthsr_r,n=n)%>%
  rename(synthvraw.age.steigers_Z=steigers_Z,synthvraw.age.steigers_pval=steigers_pval)

steigers_rawall.v.raw.age<-runcocorvec(regions=merged$brain_region,
                                          r.jk=merged$rawLFalldiff.agecor_r,
                                          r.jh=merged$rawLFdiff.agecor_r,
                                          r.kh=merged$rawLF_rawLFall_r,n=n)%>%
  rename(rawallvraw.age.steigers_Z=steigers_Z,rawallvraw.age.steigers_pval=steigers_pval)

steigers_synthall.v.raw.age<-runcocorvec(regions=merged$brain_region,
                                            r.jk=merged$synthsralldiff.agecor_r,
                                            r.jh=merged$rawLFdiff.agecor_r,
                                            r.kh=merged$rawLF_synthsrall_r,n=n)%>%
  rename(synthallvraw.age.steigers_Z=steigers_Z,synthallvraw.age.steigers_pval=steigers_pval)

#Put all the steigers tests into the aggregate dataframe
merged<-merged %>%
  left_join(select(steigers_synth.v.raw.motion,brain_region,ends_with("_Z"),ends_with("_pval"))) %>%
  left_join(select(steigers_rawall.v.raw.motion,brain_region,ends_with("_Z"),ends_with("_pval"))) %>%
  left_join(select(steigers_synthall.v.raw.motion,brain_region,ends_with("_Z"),ends_with("_pval"))) %>%
  left_join(select(steigers_synth.v.raw.age,brain_region,ends_with("_Z"),ends_with("_pval"))) %>%
  left_join(select(steigers_rawall.v.raw.age,brain_region,ends_with("_Z"),ends_with("_pval"))) %>%
  left_join(select(steigers_synthall.v.raw.age,brain_region,ends_with("_Z"),ends_with("_pval"))) %>%
  select(-contains(c("rawLF_synthsr","rawLF_rawLFall","rawLF_synthsrall")))

#Adjust for multiple comparisons for motion/demo processing
qvals<-merged %>% select(brain_region,contains("pval")) %>% 
  pivot_longer(cols=contains("pval"),names_to=c("measure",NA), names_sep="_",values_to="pval") %>%
  mutate(qval=p.adjust(as.numeric(pval),"fdr")) %>%
  pivot_wider(names_from=measure,values_from=c(pval,qval),names_glue="{measure}_{.value}")

#Merge corrected significance tests back into tibble
all_merged<-left_join(merged,qvals) %>% select(brain_region,
                                               rawLFdiff.motcor_r,rawLFdiff.motcor_pval,rawLFdiff.motcor_qval,
                                               synthsrdiff.motcor_r,synthsrdiff.motcor_pval,synthsrdiff.motcor_qval,
                                               synthvraw.mot.steigers_Z,synthvraw.mot.steigers_pval,synthvraw.mot.steigers_qval,
                                               
                                               rawLFalldiff.motcor_r,rawLFalldiff.motcor_pval,rawLFalldiff.motcor_qval,
                                               rawallvraw.mot.steigers_Z,rawallvraw.mot.steigers_pval,rawallvraw.mot.steigers_qval,
                                               
                                               synthsralldiff.motcor_r,synthsralldiff.motcor_pval,synthsralldiff.motcor_qval,
                                               synthallvraw.mot.steigers_Z,synthallvraw.mot.steigers_pval,synthallvraw.mot.steigers_qval,
                                               
                                               rawLFdiff.agecor_r,rawLFdiff.agecor_pval,rawLFdiff.agecor_qval,
                                               synthsrdiff.agecor_r,synthsrdiff.agecor_pval,synthsrdiff.agecor_qval,
                                               synthvraw.age.steigers_Z,synthvraw.age.steigers_pval,synthvraw.age.steigers_qval,
                                               
                                               rawLFalldiff.agecor_r,rawLFalldiff.agecor_pval,rawLFalldiff.agecor_qval,
                                               rawallvraw.age.steigers_Z,rawallvraw.age.steigers_pval,rawallvraw.age.steigers_qval,
                                               
                                               synthsralldiff.agecor_r,synthsralldiff.agecor_pval,synthsralldiff.agecor_qval,
                                               synthallvraw.age.steigers_Z,synthallvraw.age.steigers_pval,synthallvraw.age.steigers_qval) %>%
  mutate(across(contains('val'), DigFormat_V,cut=0.001,fixed="%.3f",scientific="%.2e"),
         across(ends_with(c("_r","_Z")),DigFormat_V,cut=0.01,fixed="%.2f",scientific="%.1e"),
         brain_region=recode(brain_region,"Mean_thick"="Mean Cortical Thickness",
                             # "Total_area"="Total WM Surface Area",
                             "Total_pial_area"="Total Surface Area",
                             "EstimatedTotalIntraCranialVol"="Estimated Intracranial Volume",
                             "SubCortGrayVol"="Subcortical Gray Matter Volume",
                             "CortexVol" = "Cortical Volume",
                             "CerebralWhiteMatterVol" = "Cerebral White Matter Volume",
                             "BrainSegVol" = "Total Brain Volume"))


#Make flextable to show motion relationship with 3T-LF correlations with different LF processing
table1<-select(all_merged,brain_region,contains("mot")) %>%
  flextable() %>% theme_booktabs() %>%
  set_header_labels(values=c('Brain measure',
                      'r','p','q',
                      'r','p','q',
                      'Z','p','q',
                      'r','p','q',
                      'Z','p','q',
                      'r','p','q',
                      'Z','p','q'))%>%
  add_header_row(values=c("Brain measure",
                          "Standard LF Axial-3T x Motion","SynthSR-3T x Motion",
                          "Steiger's Test: SynthSR-3T x Motion vs. Standard LF Axial-3T x Motion",
                          "Standard LF Multiple-3T x Motion","Steigers Test: Standard LF Multiple-3T x Motion vs. Standard LF Axial-3T x Motion",
                          "SynthSR Multiple-3T x Motion","Steiger's test: SynthSR multiple-3T x Motion vs. Standard LF Axial-3T x Motion"),
                 colwidths=c(1,3,3,3,3,3,3,3)) %>%
  vline(j=c(1,4,7,10,13,16,19),border=fp_border(color = "black",width=2)) %>% 
  merge_v(part="header") %>% italic(i=2,part="header") %>%
  bold(part="header") %>%
  border_outer(border=fp_border(color="black",width=2),part="all") %>%
  bold(i= ~ as.numeric(`rawLFdiff.motcor_qval`)<0.05, j=2:4,part="body") %>%
  bold(i= ~ as.numeric(`synthsrdiff.motcor_qval`)<0.05, j=5:7,part="body") %>%
  bold(i= ~ as.numeric(`synthvraw.mot.steigers_qval`)<0.05, j=8:10,part="body") %>%
  bold(i= ~ as.numeric(`rawLFalldiff.motcor_qval`)<0.05, j=11:13,part="body") %>%
  bold(i= ~ as.numeric(`rawallvraw.mot.steigers_qval`)<0.05, j=14:16,part="body") %>%
  bold(i= ~ as.numeric(`synthsralldiff.motcor_qval`)<0.05, j=17:19,part="body") %>%
  bold(i= ~ as.numeric(`synthallvraw.mot.steigers_qval`)<0.05, j=20:22,part="body") %>%
  fix_border_issues()
# save_as_html(table1,path="code/Diffscore_motion_cor_tests.html",title="Effects of repeated scans and SynthSR processing on the association between standardized 0.64mT-3T discrepancy and motion during scan")

#Make flextable to show age relationship with 3T-LF correlations with different LF processing

table2<-select(all_merged,brain_region,contains("age")) %>%
  flextable() %>% theme_booktabs() %>%
  set_header_labels(values=c('Brain measure',
                             'r','p','q',
                             'r','p','q',
                             'Z','p','q',
                             'r','p','q',
                             'Z','p','q',
                             'r','p','q',
                             'Z','p','q'))%>%
  add_header_row(values=c("Brain measure",
                          "Standard LF Axial-3T x Age","SynthSR-3T x Age",
                          "Steiger's Test: SynthSR-3T x Age vs. Standard LF Axial-3T x Age",
                          "Standard LF Multiple-3T x Age","Steigers Test: Standard LF Multiple-3T x Age vs. Standard LF Axial-3T x Age",
                          "SynthSR Multiple-3T x Age","Steiger's test: SynthSR multiple-3T x Age vs. Standard LF Axial-3T x Age"),
                 colwidths=c(1,3,3,3,3,3,3,3)) %>%
  vline(j=c(1,4,7,10,13,16,19),border=fp_border(color = "black",width=2)) %>% 
  merge_v(part="header") %>% italic(i=2,part="header") %>%
  bold(part="header") %>%
  # italic(i=2,part="header",j=c(-1,-2,-5,-8,-11,-14,-17,-20)) %>%
  border_outer(border=fp_border(color="black",width=2),part="all") %>%
  bold(i= ~ as.numeric(`rawLFdiff.agecor_qval`)<0.05, j=2:4,part="body") %>%
  bold(i= ~ as.numeric(`synthsrdiff.agecor_qval`)<0.05, j=5:7,part="body") %>%
  bold(i= ~ as.numeric(`synthvraw.age.steigers_qval`)<0.05, j=8:10,part="body") %>%
  bold(i= ~ as.numeric(`rawLFalldiff.agecor_qval`)<0.05, j=11:13,part="body") %>%
  bold(i= ~ as.numeric(`rawallvraw.age.steigers_qval`)<0.05, j=14:16,part="body") %>%
  bold(i= ~ as.numeric(`synthsralldiff.agecor_qval`)<0.05, j=17:19,part="body") %>%
  bold(i= ~ as.numeric(`synthallvraw.age.steigers_qval`)<0.05, j=20:22,part="body") %>%
  fix_border_issues()

save(list=c("table1","table2"),file="Diffscore_cor_tests.Rdata")

# save_as_html(table2,path="code/Diffscore_age_cor_tests.html",title="Effects of repeated scans and SynthSR processing on the association between standardized 0.64mT-3T discrepancy and age")

