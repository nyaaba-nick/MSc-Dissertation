# 	Load packages
library(netmeta)
library(meta)
library(metafor)


library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)
library(gridExtra)

library(dmetar)


library(BiocManager)
library(hasseDiagram)

library(mosaic)

##	Set working directory
setwd("Z:/Docs/i/o/Project X/m/Dissertation/Analyses")

##import data

data<-read.csv("main analyses dataset.csv")
datax <- read.csv("sensitivity analyses.csv")
datasub <- read.csv("subgroup analyses.csv")






#####Prepare Data####

########Transform meta-analysis data from two arm-based formats into contrast-based format#########


netdeath <- pairwise(treat = T,
                     event = mortality_e,
                     n = mortality_n,
                     data = data,
                     allstudies = T,
                     studlab = author_year,
                     sm = "RR")



netcoma<- pairwise(treat=T,
                   n = crt_n,
                   mean=crt_mean,
                   sd=crt_sd,
                   data = data,
                   studlab = author_year,
                   sm = "MD")


netpara<- pairwise(treat=T,
                   n = pct_n,
                   mean=pct_mean,
                   sd=pct_sd,
                   data = data,
                   studlab = author_year,
                   sm = "MD")



netneuro <- pairwise(treat = T,
                     event = neuro_e,
                     n = neuro_n,
                     data = data,
                     allstudies = T,
                     studlab = author_year,
                     sm = "RR")

nethypo <- pairwise(treat = T,
                    event = hypoglycaemia_e,
                    n = hypoglycaemia_n,
                    data = data,
                    allstudies = T,
                    studlab = author_year,
                    sm = "RR")



##################Create Factor variable for each pair of treatment########


netdeath<-as.data.frame.matrix(netdeath)


netdeath<- 
  mutate(netdeath, Pairs= derivedFactor (
    "AMI vs QN" =(t1=="AMI"&t2=="QN"),
    "AME vs QN" =(t1=="AME"&t2=="QN"),
    "ATE vs QN" =(t1=="ATE"& t2=="QN"),
    "ASU vs QN" =(t1=="ASU"& t2=="QN"),
    "ASU vs AME" =(t1=="ASU"& t2=="AME"),
    "ASU vs AMI" =(t1=="ASU"& t2=="AMI"),
    "AMI vs AME" =(t1=="AMI"& t2=="AME"),
    .method="first",
    .default=0))


#####Coma Recovery Time

netcoma<-as.data.frame.matrix(netcoma)

netcoma<- 
  mutate(netcoma, Pairs= derivedFactor (
    "AMI vs QN" =(t1=="AMI"&t2=="QN"),
    "AME vs QN" =(t1=="AME"&t2=="QN"),
    "ATE vs QN" =(t1=="ATE"& t2=="QN"),
    "ASU vs QN" =(t1=="ASU"& t2=="QN"),
    "ASU vs AME" =(t1=="ASU"& t2=="AME"),
    "ASU vs AMI" =(t1=="ASU"& t2=="AMI"),
    "AMI vs AME" =(t1=="AMI"& t2=="AME"),
    .method="first",
    .default=0))


#####Parasite Clearance Time

netpara<-as.data.frame.matrix(netpara)


netpara<- 
  mutate(netpara, Pairs= derivedFactor (
    "AMI vs QN" =(t1=="AMI"&t2=="QN"),
    "AME vs QN" =(t1=="AME"&t2=="QN"),
    "ATE vs QN" =(t1=="ATE"& t2=="QN"),
    "ASU vs QN" =(t1=="ASU"& t2=="QN"),
    "ASU vs AME" =(t1=="ASU"& t2=="AME"),
    "ASU vs AMI" =(t1=="ASU"& t2=="AMI"),
    "AMI vs AME" =(t1=="AMI"& t2=="AME"),
    .method="first",
    .default=0))


#######Neurological Sequela events


netneuro<-as.data.frame.matrix(netneuro)


netneuro<- 
  mutate(netneuro, Pairs= derivedFactor (
    "AMI vs QN" =(t1=="AMI"&t2=="QN"),
    "AME vs QN" =(t1=="AME"&t2=="QN"),
    "ATE vs QN" =(t1=="ATE"& t2=="QN"),
    "ASU vs QN" =(t1=="ASU"& t2=="QN"),
    "ASU vs AME" =(t1=="ASU"& t2=="AME"),
    "ASU vs AMI" =(t1=="ASU"& t2=="AMI"),
    "AMI vs AME" =(t1=="AMI"& t2=="AME"),
    .method="first",
    .default=0))




########Hypoglycemia events

nethypo<-as.data.frame.matrix(nethypo)

nethypo<- 
  mutate(nethypo, Pairs= derivedFactor (
    "AMI vs QN" =(t1=="AMI"&t2=="QN"),
    "AME vs QN" =(t1=="AME"&t2=="QN"),
    "ATE vs QN" =(t1=="ATE"& t2=="QN"),
    "ASU vs QN" =(t1=="ASU"& t2=="QN"),
    "ASU vs AME" =(t1=="ASU"& t2=="AME"),
    "ASU vs AMI" =(t1=="ASU"& t2=="AMI"),
    "AMI vs AME" =(t1=="AMI"& t2=="AME"),
    .method="first",
    .default=0))






########################Network Analyses#######


##Mortality
nma1 <- netmeta(TE,
                seTE,
                treat1,
                treat2,
                studlab,
                netdeath,
                reference.group = "Quinine",
                sm = "RR", 
                sep.trt=" vs ",
                comb.fixed = FALSE,
                prediction = T)


##Coma Recovery Time
nma2<-netmeta(TE,
              seTE,
              treat1,
              treat2,
              studlab,
              netcoma,
              sm = "MD",
              comb.fixed = FALSE,
              prediction=T,
              reference.group = "Quinine",
              sep.trt=" vs ",
              details.chkmultiarm = T)



##Parasite clearance time
nma3<-netmeta(TE,
              seTE,
              treat1,
              treat2,
              studlab,
              netpara,
              sm = "MD",
              comb.fixed = FALSE,
              prediction=T,
              reference.group = "Quinine",
              sep.trt=" vs " ,
              details.chkmultiarm = T)


##Neurological Sequela Events
nma4 <- netmeta(TE,
                seTE,
                treat1,
                treat2,
                studlab,
                netneuro,
                sm = "RR",
                comb.fixed = FALSE,
                prediction=T,
                reference.group = "Quinine",
                sep.trt="vs",
                details.chkmultiarm = T)












 ########################Using par(mfrow ) for forest plots############################



#par(mfrow ) seems not to be working with the forest plots######
par(mfrow = c(2, 2))

###Mortality##

forest(nma1,
       xlim = c(0.3,2),
       smlab = "Mortality",
       small.values = "good",
       reference.group = "Quinine",
       just.smlab = "left",
       leftcols =c("studlab","Pscore"),
       sortvar = -Pscore,
       col.square= "blue",
       col.square.lines = "blue",
       fontsize =12,
       plotwidth = "8cm",
       colgap.left="4mm",
       spacing= 2,
       print.I2= F,
       print.pval.Q=F,
       hetlab= " ",
       label.left = "favours treatment",
       label.right= "favours quinine",
       fs.smlab = 13, ff.smlab = "bold",
       fs.study.labels= 12, ff.study.labels = "bold",
       fs.axis = 9,
       fs.addline = 9, ff.addline = "italic",
       text.addline1="*I^2=0%, p=0.67"
)


###CRT###
forest(nma2,
       xlim = c(-40,30),
       smlab = "Coma Recovery Time",
       small.values = "good",
       reference.group = "Quinine",
       just.smlab = "left",
       leftcols =c("studlab","Pscore"),
       sortvar = -Pscore,
       col.square= "blue",
       col.square.lines = "blue",
       fontsize =12,
       plotwidth = "8cm",
       colgap.left="4mm",
       spacing= 2,
       print.I2= F,
       print.pval.Q=F,
       hetlab= " ",
       label.left = "favours treatment",
       label.right= "favours quinine",
       fs.smlab = 13, ff.smlab = "bold",
       fs.study.labels= 12, ff.study.labels = "bold",
       fs.axis = 9,
       fs.addline = 9, ff.addline = "italic",
       text.addline1=" *I^2=53.0%, p=0.006"
)

###PCT####
forest(nma3,
       xlim = c(-40,20),
       smlab = "Parasite Clearance Time",
       small.values = "good",
       reference.group = "Quinine",
       just.smlab = "left",
       leftcols =c("studlab","Pscore"),
       sortvar = -Pscore,
       col.square= "blue",
       col.square.lines = "blue",
       fontsize =12,
       plotwidth = "8cm",
       colgap.left="4mm",
       spacing= 2,
       print.I2= F,
       print.pval.Q=F,
       hetlab= " ",
       label.left = "favours treatment",
       label.right= "favours quinine",
       fs.smlab = 13, ff.smlab = "bold",
       fs.study.labels= 12, ff.study.labels = "bold",
       fs.axis = 9,
       fs.addline = 9, ff.addline = "italic",
       text.addline1= "*I^2 == 87.7% p <0.001"
       
)



###NSE####
forest(nma4,
       xlim = c(0.05,80),
       smlab = "Neurological Sequela Events",
       small.values = "good",
       reference.group = "Quinine",
       just.smlab = "left",
       leftcols =c("studlab","Pscore"),
       sortvar = -Pscore,
       col.square= "blue",
       col.square.lines = "blue",
       fontsize =12,
       plotwidth = "8cm",
       colgap.left="4mm",
       spacing= 2,
       print.I2= F,
       print.pval.Q=F,
       hetlab= " ",
       label.left = "favours treatment",
       label.right= "favours quinine",
       fs.smlab = 13, ff.smlab = "bold",
       fs.study.labels= 12, ff.study.labels = "bold",
       fs.axis = 9,
       fs.addline = 9, ff.addline = "italic",
       text.addline1=" *I^2=11.7%, p=0.33"
)



#################################Using par(mfrow ) for net heat plots ############################


#subset data
child<- subset(netpara,age_group=="c")
#Fit network consistency model  for children
nma3c<-netmeta(TE,
               seTE,
               treat1,
               treat2,
               studlab,
               child,
               sm = "MD",
               comb.fixed = FALSE,
               prediction=T,
               reference.group = "Quinine",
               sep.trt=" vs ",
               details.chkmultiarm = T)



#subset data
adult<- subset(netpara,age_group=="a")

#Fit network consistency model for adult
nma3a<-netmeta(TE,
               seTE,
               treat1,
               treat2,
               studlab,
               adult,
               sm = "MD",
               comb.fixed = FALSE,
               prediction=T,
               reference.group = "Quinine",
               sep.trt=" vs ",
               details.chkmultiarm = T)




#par(mfrow ) seems not to be working with the forest plots######


par(mfrow = c(2, 2))


#adults
netheat(nma3a,
        showall = T,
        random = T)

#children
netheat(nma3c,
        showall = T,
        random = T)


#lumped population
netheat(nma3,
        showall = T,
        random = T)






#################adjusting the resolution of plots using tiff()#########



tiff("plot5.tiff",
     height = 480,
     width=480,
     units = "px")


netgraph(nma1,
         plastic = FALSE,
         start.layout = "random",
         labels = c("Artemisinin; ASN" ,"Artemether; ATM","Artesunate; ATS","Arteether; ATT","Quinine; QN"),
         points = TRUE,
         cex.points = nodesize,
         offset = 0.068,
         col.points = "blue",
         col = "black",
         number.of.studies = T,
         col.number.of.studies = "grey",
         thickness = "number.of.studies")
