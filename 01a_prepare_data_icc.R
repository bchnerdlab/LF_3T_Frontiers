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
library(psych)

##NOTE: This script is not fully commented.  It is a subset of the operations performed in 01_prepare_data.R
##      It solely replaces Pearson's correlations with ICCs.  Check that script's comments for more information.


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

runicc<-function(xdata,ydata,regions){
  rows<-length(regions)
  out.df<-as.data.frame(matrix(nrow=rows,ncol=6))
  colnames(out.df)<-c("brain_region","t","r","pval","CIlow","CIhigh")
  brain.reg<-paste(regions,sep=",")
  
  for (i in 1:length(brain.reg)){
    var=brain.reg[[i]]
    save<-ICC(data.frame(xdata[[var]],ydata[[var]]),lmer=FALSE)$results %>% filter(type=="ICC3")
    out.df[i,"t"]<-save$F[1]
    out.df[i,"r"]<-save$ICC[1]
    out.df[i,"pval"]<-save$p[1]
    out.df[i,"CIlow"]<-save$`lower bound`[1]
    out.df[i,"CIhigh"]<-save$`upper bound`[1]
    out.df[i,"brain_region"]<-brain.reg[i]
  }
  return(out.df)
}

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
#######################
# READ IN REGIONAL DATA
#######################

#Make list of excluded participants and identify which aseg2table and aparc2table files to use
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
# LHSA<-lh.SA %>% select(1,lh_WhiteSurfArea_area)%>% rename("SubjID"=1)
# RHSA<-rh.SA %>% select(1,rh_WhiteSurfArea_area)%>% rename("SubjID"=1)
RHPSA<-read_lines("../data/extract_pialsurf_rh.txt") %>% str_replace(.,":",",") %>% 
  str_remove("/stats/rh.aparc.pial.stats")%>% as.tibble() %>% 
  separate(col=1,into=c("SubjID",NA,NA,NA,"rh.PialSurfArea",NA),sep=",") %>% mutate(rh.PialSurfArea=as.numeric(rh.PialSurfArea))
LHPSA<-read_lines("../data/extract_pialsurf_lh.txt") %>% str_replace(.,":",",") %>% 
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

#################
# GLOBAL PROCESSING
#################
all_merged<-tibble(brain_region=character())
i="-axi[12]-1to1"

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
savelist=c("global.3Trun01","global.synthsr","global.rawLF","global.rawLF.all","global.synthsr.all","vars","GlobalVals")
test<-nrow(global.3Trun01)==nrow(global.synthsr) & nrow(global.synthsr)==nrow(global.rawLF)
# save(list=savelist,file="global_data.Rdata")

if (test) { 
  n=nrow(global.3Trun01)
  print(paste0("Total N = ",n))
  # write(global.3Trun01$SubjID,file="G:/Shared drives/NeRD_lab/research_projects/LF_3T/data/sublist.txt")
}

#AIm 1 - axi1 raw vs axi1 synth
synth.3T.cor<-runicc(global.3Trun01,global.synthsr,regions=vars[-1]) %>%
  rename(synth3T_r=r,synth3T_pval=pval)

raw.3T.cor<-runicc(global.3Trun01,global.rawLF,regions=vars[-1]) %>%
  rename(raw3T_r=r,raw3T_pval=pval)

raw.synth.cor<-runicc(global.synthsr, global.rawLF,regions=vars[-1]) %>%
  rename(synthraw_r=r,synthraw_pval=pval)

#Aim 2 - axi1 raw vs all raw
rawall.3T.cor<-runicc(global.3Trun01,global.rawLF.all,regions=vars[-1]) %>%
  rename(rawall3T_r=r,rawall3T_pval=pval)

raw.rawall.cor<-runicc(global.rawLF,global.rawLF.all,regions=vars[-1]) %>%
  rename(rawallraw_r=r,rawallraw_pval=pval)

#Aim 3 - axi1 raw vs all synth
synthall.3T.cor<- runicc(global.3Trun01,global.synthsr.all,regions=vars[-1]) %>%
  rename(synthall3T_r=r,synthall3T_pval=pval)

raw.synthall.cor<-runicc(global.synthsr.all, global.rawLF,regions=vars[-1]) %>%
  rename(synthallraw_r=r,synthallraw_pval=pval)

#Aim 4 - axi1 synth vs all synth
synth.synthall.cor<-runicc(global.synthsr,global.synthsr.all,regions=vars[-1]) %>%
  rename(synthallsynth_r=r,synthallsynth_pval=pval)

#Merge all together
merged<-tibble(brain_region=vars[-1])

merged<-merged %>% left_join(select(synth.3T.cor,brain_region,synth3T_r,synth3T_pval)) %>%
  left_join(select(raw.3T.cor,brain_region,raw3T_r,raw3T_pval)) %>%
  left_join(select(raw.synth.cor,brain_region,synthraw_r,synthraw_pval)) %>% 
  left_join(select(rawall.3T.cor,brain_region,rawall3T_r,rawall3T_pval)) %>% 
  left_join(select(raw.rawall.cor,brain_region,rawallraw_r,rawallraw_pval)) %>% 
  left_join(select(synthall.3T.cor,brain_region,synthall3T_r,synthall3T_pval)) %>% 
  left_join(select(raw.synthall.cor,brain_region,synthallraw_r,synthallraw_pval)) %>%
  left_join(select(synth.synthall.cor,brain_region,synthallsynth_r,synthallsynth_pval))

#Aim 1 tests
steigers_aim1<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$synth3T_r,r.jh=merged$raw3T_r,r.kh=merged$synthraw_r,n=n) %>%
  rename(aim1steigers_Z=steigers_Z,aim1steigers_pval=steigers_pval)

# mape_synth_axi1<-runmape(actual=global.3Trun01, predicted=global.synthsr, regions=vars[-1]) %>%
#   rename(synthaxi13T_mape=mape)
# mape_raw_axi1<-runmape(actual=global.3Trun01, predicted=global.rawLF, regions=vars[-1]) %>%
#   rename(rawaxi13T_mape=mape)

#Aim 2 tests
steigers_aim2<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$rawall3T_r,r.jh=merged$raw3T_r,r.kh=merged$rawallraw_r,n=n) %>%
  rename(aim2steigers_Z=steigers_Z,aim2steigers_pval=steigers_pval)
# mape_raw_all<-runmape(actual=global.3Trun01, predicted=global.rawLF.all, regions=vars[-1]) %>%
#   rename(rawall3T_mape=mape)

#Aim 3 tests
steigers_aim3<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$synthall3T_r,r.jh=merged$raw3T_r,r.kh=merged$synthallraw_r,n=n) %>%
  rename(aim3steigers_Z=steigers_Z,aim3steigers_pval=steigers_pval)
# mape_synth_all<-runmape(actual=global.3Trun01, predicted=global.synthsr.all, regions=vars[-1]) %>%
#   rename(synthall3T_mape=mape)

steigers_aim4<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$synthall3T_r,r.jh=merged$synth3T_r,r.kh=merged$synthallsynth_r,n=n) %>%
  rename(aim4steigers_Z=steigers_Z,aim4steigers_pval=steigers_pval)




merged<-left_join(merged,steigers_aim1) %>% 
  # left_join(mape_synth_axi1) %>% 
  # left_join(mape_raw_axi1) %>% 
  left_join(steigers_aim2) %>% 
  # left_join(mape_raw_all) %>% 
  left_join(steigers_aim3) %>% 
  # left_join(mape_synth_all) %>% 
  left_join(steigers_aim4) %>%
    select(-contains(c("synthraw","rawallraw","synthallraw")))




qvals<-merged %>% select(brain_region,contains("pval")) %>% 
  pivot_longer(cols=contains("pval"),names_to=c("measure",NA), names_sep="_",values_to="pval") %>%
  mutate(qval=p.adjust(as.numeric(pval),"fdr")) %>%
  pivot_wider(names_from=measure,values_from=c(pval,qval),names_glue="{measure}_{.value}")

all_merged<-left_join(merged,qvals) %>% select(brain_region,
                                               raw3T_r,raw3T_pval,raw3T_qval,#rawaxi13T_mape,
                                               synth3T_r,synth3T_pval,synth3T_qval,#synthaxi13T_mape,
                                               aim1steigers_Z,aim1steigers_pval,aim1steigers_qval,
                                               rawall3T_r,rawall3T_pval,rawall3T_qval,#rawall3T_mape,
                                               aim2steigers_Z,aim2steigers_pval,aim2steigers_qval,
                                               synthall3T_r,synthall3T_pval,synthall3T_qval,#synthall3T_mape,
                                               aim3steigers_Z,aim3steigers_pval,aim3steigers_qval,
                                               aim4steigers_Z,aim4steigers_pval,aim4steigers_qval) %>%
  mutate(across(contains('val'), DigFormat_V,cut=0.001,fixed="%.3f",scientific="%.2e"),
        across(ends_with(c("_r","_Z")),DigFormat_V,cut=0.01,fixed="%.2f",scientific="%.1e"),
        #across(contains("mape"),DigFormat_V,cut=.001,fixed="%.2f",scientific="%.2f"),
        brain_region=recode(brain_region,"Mean_thick"="Mean Cortical Thickness",
                             # "Total_area"="Total WM Surface Area",
                            "Total_pial_area"="Total Surface Area",
                             "EstimatedTotalIntraCranialVol"="Estimated Intracranial Volume",
                             "SubCortGrayVol"="Subcortical Gray Matter Volume",
                              "CortexVol" = "Cortical Volume",
                            "CerebralWhiteMatterVol" = "Cerebral White Matter Volume",
                            "BrainSegVol" = "Total Brain Volume"))

save(all_merged,file="icc_global_all_merged.Rdat")

load("icc_global_all_merged.Rdat")
vars<-c("SubjID","Total_pial_area","Mean_thick","EstimatedTotalIntraCranialVol",
        "SubCortGrayVol","CortexVol","CerebralWhiteMatterVol","BrainSegVol")


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

insertdf<-tibble(comp=c("aim1","aim1","aim2","aim2","aim3","aim3"),
                 brain_region=rep("Brain measure",times=6),
                 raw3T_r=rep(c("Standard Axial 64mT ICCs with 3T","ICC"),times=3),
                 raw3T_pval=rep(c("Standard Axial 64mT ICCs with 3T","p"),times=3),
                 raw3T_qval=rep(c("Standard Axial 64mT ICCs with 3T","q"),times=3),
                 comparison_r=c("SynthSR-Processed Axial 64mT ICCs with 3T","ICC",
                                "Standard Multi-Orientation 64mT ICCs with 3T","ICC",
                                "SynthSR-Processed Multi-Orientation 64mT ICCs with 3T","ICC"),
                 comparison_pval=c("SynthSR-Processed Axial 64mT ICCs with 3T","p",
                                "Standard Multi-Orientation 64mT ICCs with 3T","p",
                                "SynthSR-Processed Multi-Orientation 64mT ICCs with 3T","p"),
                 comparison_qval=c("SynthSR-Processed Axial 64mT ICCs with 3T","q",
                                "Standard Multi-Orientation 64mT ICCs with 3T","q",
                                "SynthSR-Processed Multi-Orientation 64mT ICCs with 3T","q"),
                 comparison_steigersZ=rep(c("Steiger's","Z"),times=3),
                 comparison_steigerspval=rep(c("Steiger's","p"),times=3),
                 comparison_steigersqval=rep(c("Steiger's","q"),times=3)
                 
                 )


  
  
test2<-flexprep %>% 
  insertRows(r=c(0,1,length(vars),length(vars),2*length(vars)-1,2*length(vars)-1),rcurrent=TRUE,new=insertdf)

fullflex<-test2 %>% select(-1) %>% flextable () %>% 
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
  add_header_lines(values="Table 3. Intra-class correlations (ICCs) between low-field (64mT) and high-field (3T) MR images and comparisons across super-resolution approaches")
  
save_as_html(fullflex,path="table3.html")


#################
# REGIONAL PROCESSING
#################
all_merged<-tibble(brain_region=character())
i="-axi[12]-1to1"

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
savelist=c("regional.3Trun01","regional.synthsr","regional.rawLF","regional.rawLF.all","regional.synthsr.all","allregs")
save(list=savelist,file="regional_data.Rdata")

test<-nrow(regional.3Trun01)==nrow(regional.synthsr) & nrow(regional.synthsr)==nrow(regional.rawLF)

if (test) { 
  n=nrow(regional.3Trun01)
  print(paste0("Total N = ",n))
}

#AIm 1 - axi1 raw vs axi1 synth
synth.3T.cor<-runicc(regional.3Trun01,regional.synthsr,regions=allregs) %>%
  rename(synth3T_r=r,synth3T_pval=pval)

raw.3T.cor<-runicc(regional.3Trun01,regional.rawLF,regions=allregs) %>%
  rename(raw3T_r=r,raw3T_pval=pval)

raw.synth.cor<-runicc(regional.synthsr, regional.rawLF,regions=allregs) %>%
  rename(synthraw_r=r,synthraw_pval=pval)

#Aim 2 - axi1 raw vs all raw
rawall.3T.cor<-runicc(regional.3Trun01,regional.rawLF.all,regions=allregs) %>%
  rename(rawall3T_r=r,rawall3T_pval=pval)

raw.rawall.cor<-runicc(regional.rawLF,regional.rawLF.all,regions=allregs) %>%
  rename(rawallraw_r=r,rawallraw_pval=pval)

#Aim 3 - axi1 raw vs all synth
synthall.3T.cor<- runicc(regional.3Trun01,regional.synthsr.all,regions=allregs) %>%
  rename(synthall3T_r=r,synthall3T_pval=pval)

raw.synthall.cor<-runicc(regional.synthsr.all, regional.rawLF,regions=allregs) %>%
  rename(synthallraw_r=r,synthallraw_pval=pval)

#Aim 4 - axi1 synth vs all synth
synth.synthall.cor<-runicc(regional.synthsr,regional.synthsr.all,regions=allregs) %>%
  rename(synthallsynth_r=r,synthallsynth_pval=pval)
#Merge all together
merged<-tibble(brain_region=allregs)

merged<-merged %>% left_join(select(synth.3T.cor,brain_region,synth3T_r,synth3T_pval)) %>%
  left_join(select(raw.3T.cor,brain_region,raw3T_r,raw3T_pval)) %>%
  left_join(select(raw.synth.cor,brain_region,synthraw_r,synthraw_pval)) %>% 
  left_join(select(rawall.3T.cor,brain_region,rawall3T_r,rawall3T_pval)) %>% 
  left_join(select(raw.rawall.cor,brain_region,rawallraw_r,rawallraw_pval)) %>% 
  left_join(select(synthall.3T.cor,brain_region,synthall3T_r,synthall3T_pval)) %>% 
  left_join(select(raw.synthall.cor,brain_region,synthallraw_r,synthallraw_pval)) %>%
  left_join(select(synth.synthall.cor,brain_region,synthallsynth_r,synthallsynth_pval))

#Aim 1 tests
steigers_aim1<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$synth3T_r,r.jh=merged$raw3T_r,r.kh=merged$synthraw_r,n=n) %>%
  rename(aim1steigers_Z=steigers_Z,aim1steigers_pval=steigers_pval)
# mape_synth_axi1<-runmape(actual=regional.3Trun01, predicted=regional.synthsr, regions=vars[-1]) %>%
#   rename(synthaxi13T_mape=mape)
# mape_raw_axi1<-runmape(actual=regional.3Trun01, predicted=regional.rawLF, regions=vars[-1]) %>%
#   rename(rawaxi13T_mape=mape)

#Aim 2 tests
steigers_aim2<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$rawall3T_r,r.jh=merged$raw3T_r,r.kh=merged$rawallraw_r,n=n) %>%
  rename(aim2steigers_Z=steigers_Z,aim2steigers_pval=steigers_pval)
# mape_raw_all<-runmape(actual=regional.3Trun01, predicted=regional.rawLF.all, regions=vars[-1]) %>%
#   rename(rawall3T_mape=mape)

#Aim 3 tests
steigers_aim3<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$synthall3T_r,r.jh=merged$raw3T_r,r.kh=merged$synthallraw_r,n=n) %>%
  rename(aim3steigers_Z=steigers_Z,aim3steigers_pval=steigers_pval)
# mape_synth_all<-runmape(actual=regional.3Trun01, predicted=regional.synthsr.all, regions=vars[-1]) %>%
#   rename(synthall3T_mape=mape)

steigers_aim4<-runcocorvec(regions=merged$brain_region,
                           r.jk=merged$synthall3T_r,r.jh=merged$synth3T_r,r.kh=merged$synthallsynth_r,n=n) %>%
  rename(aim4steigers_Z=steigers_Z,aim4steigers_pval=steigers_pval)


merged<-left_join(merged,steigers_aim1) %>% 
  # left_join(mape_synth_axi1) %>% 
  # left_join(mape_raw_axi1) %>% 
  left_join(steigers_aim2) %>% 
  # left_join(mape_raw_all) %>% 
  left_join(steigers_aim3) %>% 
  # left_join(mape_synth_all) %>% 
  left_join(steigers_aim4) %>%   
  select(-contains(c("synthraw","rawallraw","synthallraw")))




qvals<-merged %>% select(brain_region,contains("pval")) %>% 
  pivot_longer(cols=contains("pval"),names_to=c("measure",NA), names_sep="_",values_to="pval") %>%
  mutate(qval=p.adjust(as.numeric(pval),"fdr")) %>%
  pivot_wider(names_from=measure,values_from=c(pval,qval),names_glue="{measure}_{.value}")

all_merged<-left_join(merged,qvals) %>% select(brain_region,
                                               raw3T_r,raw3T_pval,raw3T_qval,#rawaxi13T_mape,
                                               synth3T_r,synth3T_pval,synth3T_qval,#synthaxi13T_mape,
                                               aim1steigers_Z,aim1steigers_pval,aim1steigers_qval,
                                               rawall3T_r,rawall3T_pval,rawall3T_qval,#rawall3T_mape,
                                               aim2steigers_Z,aim2steigers_pval,aim2steigers_qval,
                                               synthall3T_r,synthall3T_pval,synthall3T_qval,#synthall3T_mape,
                                               aim3steigers_Z,aim3steigers_pval,aim3steigers_qval,
                                               aim4steigers_Z,aim4steigers_pval,aim4steigers_qval) %>%
  mutate(across(contains('val'), DigFormat_V,cut=0.001,fixed="%.3f",scientific="%.2e"),
         across(contains(c("_r$","_Z$")),\(x) round(x, digits=2))) #,
         #across(contains("mape"),DigFormat_V,cut=.001,fixed="%.2f",scientific="%.2f"),
         # brain_region=recode(brain_region,"Mean_thick"="Mean Cortical Thickness",
         #                     "Total_area"="Total Surface Area",
         #                     "EstimatedTotalIntraCranialVol"="Estimated Intracranial Volume",
         #                     "SubCortGrayVol"="Subcortical Gray Matter Volume"))

save(all_merged,file="icc_regional_all_merged.Rdat")

load("icc_regional_all_merged.Rdat")
SA_dk_ref<-read_delim(file='../code/ref/SA_ggseg_ref_desikan.txt',delim="\t") ##ggseg reference
CT_dk_ref<-read_delim(file='../code/ref/CT_ggseg_ref_desikan.txt',delim="\t") ##ggseg reference
CV_dk_ref<-read_delim(file='../code/ref/vol_ggseg_ref_desikan.txt',delim="\t") ##ggseg reference
sub_dk_ref<-read_delim(file='../code/ref/sub_ggseg_ref_aseg_new.txt',delim="\t") ##ggseg reference

###############################################
#MAKE NEW FLEXTABLE STABLES
###############################################
all_ref<-full_join(select(SA_dk_ref,find,hemi,replacewith),select(CT_dk_ref,find,hemi,replacewith)) %>%
  full_join(select(CV_dk_ref,find,hemi,replacewith)) %>%
  full_join(select(sub_dk_ref,find,hemi,replacewith)) %>%
  mutate(hemi=ifelse(is.na(hemi),"left",hemi),
         replacewith=paste(hemi,replacewith,sep=" "),
         measure=case_when(str_detect(find,"_areapial") ~ "area",
                           str_detect(find,"_thickness") ~ "thickness",
                           str_detect(find,"_volume") ~ "volume",
                           TRUE ~ "subcortical volume")) %>% 
  select(-hemi) %>% rename("brain_region"=find) 


flexAim1<-left_join(all_merged, all_ref) %>% unique() %>%
  mutate(brain_region=ifelse(!is.na(replacewith),replacewith,tolower(str_replace_all(brain_region,"[-_]"," "))),
         brain_region=str_remove(brain_region,"midline "),
         measure=ifelse(is.na(measure),"subcortical volume",measure),
         brain_region=paste(brain_region,measure,sep=" "),
         brain_region=str_replace(brain_region,"cc ","corpus callosum "),
         brain_region=str_replace(brain_region,"DC ","diencephalon "),
         brain_region=str_replace(brain_region,"bankssts","banks of superior temporal sulcus")) %>%
  mutate(across(ends_with("_r"),DigFormat_V,cut=0.01,fixed="%.2f",scientific="%.1e"),
         across(ends_with("_Z"),DigFormat_V,cut=0.01,fixed="%.2f",scientific="%.1e"),
         across(ends_with("val"),DigFormat_V,cut=0.001,fixed="%.3f",scientific="%.2e")) %>%
  select(brain_region,
         raw3T_r,raw3T_pval,raw3T_qval,#rawaxi13T_mape,
         synth3T_r,synth3T_pval,synth3T_qval,#synthaxi13T_mape,
         aim1steigers_Z,aim1steigers_pval,aim1steigers_qval) %>%
  flextable()  %>% theme_booktabs() %>% fontsize(size=9,part="all") %>%
  autofit() %>%
  vline(j=c(1,4,7),border=fp_border(color = "black",width=2)) %>% 
  hline(i=seq(from=68,to=242,by=68),border=fp_border(color = "black",width=2)) %>%
  border_outer(border=fp_border(color="black",width=2),part="all") %>%
  bold(i= ~ as.numeric(`raw3T_qval`)<0.05, j=2:4,part="body") %>%
  bold(i= ~ as.numeric(`synth3T_qval`)<0.05, j=5:7,part="body") %>%
  bold(i= ~ as.numeric(`aim1steigers_qval`)<0.05, j=8:10,part="body") %>%
  bold(part="header") %>%
  add_header_row(top=TRUE,values=c("Measurement",
                                   "Standard Axial 64mT Correlations with 3T",
                                   "SynthSR−Processed Axial 64mT Correlations with 3T",
                                   "Steiger"),
                 colwidths = c(1,3,3,3)) %>% 
  set_header_labels(values=c("Measurement",
                             "ICC","p","q",#"MAPE",
                             "ICC","p","q",#"MAPE",
                             "ICC","p","q")) %>%
  italic(i=2,part="header") %>% merge_v(part="header") %>%# merge_v(j=1,part="body") %>% 
  # hline(i=1:4,j=1,border=fp_border(color = "black",width=2)) %>%
  align(align="center",part="all") %>% #italic(i=2,j=c(5,9),italic=FALSE,part="header") %>%
  fix_border_issues()%>%
  add_header_lines(values="STable 2. ICC of regional measurements from standard versus SynthSR-processed axial 64mT scans with 3T scans.")

save_as_html(flexAim1,path="ICC_STable2.html")


flexAim2<-left_join(all_merged, all_ref) %>% unique() %>%
  mutate(brain_region=ifelse(!is.na(replacewith),replacewith,tolower(str_replace_all(brain_region,"[-_]"," "))),
         brain_region=str_remove(brain_region,"midline "),
         measure=ifelse(is.na(measure),"subcortical volume",measure),
         brain_region=paste(brain_region,measure,sep=" "),
         brain_region=str_replace(brain_region,"cc ","corpus callosum "),
         brain_region=str_replace(brain_region,"DC ","diencephalon "),
         brain_region=str_replace(brain_region,"bankssts","banks of superior temporal sulcus")) %>%
  mutate(across(ends_with("_r"),DigFormat_V,cut=0.01,fixed="%.2f",scientific="%.1e"),
         across(ends_with("_Z"),DigFormat_V,cut=0.01,fixed="%.2f",scientific="%.1e"),
         across(ends_with("val"),DigFormat_V,cut=0.001,fixed="%.3f",scientific="%.2e")) %>%
  select(brain_region,
         raw3T_r,raw3T_pval,raw3T_qval,#rawaxi13T_mape,
         rawall3T_r,rawall3T_pval,rawall3T_qval,#synthaxi13T_mape,
         aim2steigers_Z,aim2steigers_pval,aim2steigers_qval) %>%
  flextable()  %>% theme_booktabs() %>% fontsize(size=9,part="all") %>%
  autofit() %>%
  vline(j=c(1,4,7),border=fp_border(color = "black",width=2)) %>% 
  hline(i=seq(from=68,to=242,by=68),border=fp_border(color = "black",width=2)) %>%
  border_outer(border=fp_border(color="black",width=2),part="all") %>%
  bold(i= ~ as.numeric(`raw3T_qval`)<0.05, j=2:4,part="body") %>%
  bold(i= ~ as.numeric(`rawall3T_qval`)<0.05, j=5:7,part="body") %>%
  bold(i= ~ as.numeric(`aim2steigers_qval`)<0.05, j=8:10,part="body") %>%
  bold(part="header") %>%
  add_header_row(top=TRUE,values=c("Measurement",
                                   "Standard Axial 64mT Correlations with 3T",
                                   "Standard Multi-Orientation 64mT Correlations with 3T",
                                   "Steiger"),
                 colwidths = c(1,3,3,3)) %>% 
  set_header_labels(values=c("Measurement",
                             "ICC","p","q",#"MAPE",
                             "ICC","p","q",#"MAPE",
                             "z","p","q")) %>%
  italic(i=2,part="header") %>% merge_v(part="header") %>%# merge_v(j=1,part="body") %>% 
  # hline(i=1:4,j=1,border=fp_border(color = "black",width=2)) %>%
  align(align="center",part="all") %>% #italic(i=2,j=c(5,9),italic=FALSE,part="header") %>%
  fix_border_issues()%>%
  add_header_lines(values="STable 3. ICC of regional measurements from standard axial versus standard multi-orientation 64mT scans with 3T scans.")

save_as_html(flexAim2,path="ICC_STable3.html")