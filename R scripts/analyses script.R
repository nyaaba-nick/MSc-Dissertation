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







                                   #######Traditional Paiwise Meta_Analyses using "metagen" with REML and QP######

##Mortality
madeath<- metagen(TE,
               seTE,
               studlab,
               data = netdeath,
               byvar = netdeath$Pairs, 
               sm = "RR", 
               comb.fixed = FALSE,
               prediction = T,
               method.tau= "REML",
               method.tau.ci = "QP")


#display model output
summary(madeath)

#forest plot

forest(madeath, 
       xlim = c(0.01, 50), 
       smlab = "Direct Comparison of Mortality Rates",
       label.left = "reduces mortality",
       label.right= "increase mortality",
       
       
       test.subgroup.random=F,
       test.effect.subgroup.random=F,
       small.values = "good",
       rightcols = c("effect","ci"),
       leftcols = "studlab",
       just.smlab = "left",
       col.square.lines = "blue",
       fontsize =8,
       plotwidth = "10cm",
       colgap.left="8mm",
       spacing= 0.75,
       overall = F,
       overall.hetstat = F,
       prediction = F,
       studlab = T,
       study.results = T,
       pooled.totals = T,
       hetstat = F,
       
       fs.smlab = 10, ff.smlab = "bold",
       fs.study.labels= 7, ff.study.labels = "bold",
       fs.axis = 8,
       fs.xlab = 10, ff.xlab = "bold.italic",
       fs.heading= 10, ff.heading = "bold",
       fs.random = 10, ff.random = "bold",
       text.random.w=("")
)


##Coma Recovery Time
macrt<- metagen(TE,
               seTE,
               studlab,
               data = netcoma,
               byvar = netcoma$Pairs, 
               sm = "MD", 
               comb.fixed = FALSE,
               prediction = T,
               method.tau= "REML",
               method.tau.ci = "QP")

#display model output
summary(macrt)

#forest plot

forest(macrt, 
       xlim = c(-100, 150), 
       smlab = "Direct Comparison of Coma Recovery Time ",
       label.left = "reduces recovery time",
       label.right= "increases recovery time",
       
       
       test.subgroup.random=F,
       test.effect.subgroup.random=F,
       small.values = "good",
       rightcols = c("effect","ci"),
       leftcols = "studlab",
       just.smlab = "left",
       col.square.lines = "blue",
       fontsize =8,
       plotwidth = "10cm",
       colgap.left="8mm",
       spacing= 0.75,
       overall = F,
       overall.hetstat = F,
       prediction = F,
       studlab = T,
       study.results = T,
       pooled.totals = T,
       hetstat = F,
       
       fs.smlab = 10, ff.smlab = "bold",
       fs.study.labels= 7, ff.study.labels = "bold",
       fs.axis = 8,
       fs.xlab = 10, ff.xlab = "bold.italic",
       fs.heading= 10, ff.heading = "bold",
       fs.random = 10, ff.random = "bold",
       text.random.w=("")
)


##Parasite clearance time
mapct<- metagen(TE,
               seTE,
               studlab,
               data = netpara,
               byvar = netpara$Pairs, 
               sm = "MD", 
               comb.fixed = FALSE,
               prediction = T,
               method.tau= "REML",
               method.tau.ci = "QP")

#display model output
summary(mapct)

#forest plot
forest(mapct, 
       xlim = c(-80, 20), 
       smlab = "Direct Comparison of Parasite Clearance Time ",
       label.left = "reduces clearance time",
       label.right= "increases clearance time",
       
       
       test.subgroup.random=F,
       test.effect.subgroup.random=F,
       small.values = "good",
       rightcols = c("effect","ci"),
       leftcols = "studlab",
       just.smlab = "left",
       col.square.lines = "blue",
       fontsize =8,
       plotwidth = "10cm",
       colgap.left="8mm",
       spacing= 0.75,
       overall = F,
       overall.hetstat = F,
       prediction = F,
       studlab = T,
       study.results = T,
       pooled.totals = T,
       hetstat = F,
       
       fs.smlab = 10, ff.smlab = "bold",
       fs.study.labels= 7, ff.study.labels = "bold",
       fs.axis = 8,
       fs.xlab = 10, ff.xlab = "bold.italic",
       fs.heading= 10, ff.heading = "bold",
       fs.random = 10, ff.random = "bold",
       text.random.w=("")
)


##Neurological Sequela Events
maneuro<- metagen(TE,
               seTE,
               studlab,
               data = netneuro,
               byvar = netneuro$Pairs, 
               sm = "RR", 
               comb.fixed = FALSE,
               prediction = T,
               method.tau= "REML",
               method.tau.ci = "QP")

#display model output
summary(maneuro)

#forest plot
forest(maneuro, 
       xlim = c(0.01, 70), 
       smlab = "Direct Comparison of Neurological Sequela Events ",
       label.left = "reduces events",
       label.right= "increases events",
       
       
       test.subgroup.random=F,
       test.effect.subgroup.random=F,
       small.values = "good",
       rightcols = c("effect","ci"),
       leftcols = "studlab",
       just.smlab = "left",
       col.square.lines = "blue",
       fontsize =8,
       plotwidth = "10cm",
       colgap.left="8mm",
       spacing= 0.75,
       overall = F,
       overall.hetstat = F,
       prediction = F,
       studlab = T,
       study.results = T,
       pooled.totals = T,
       hetstat = F,
       
       fs.smlab = 10, ff.smlab = "bold",
       fs.study.labels= 7, ff.study.labels = "bold",
       fs.axis = 8,
       fs.xlab = 10, ff.xlab = "bold.italic",
       fs.heading= 10, ff.heading = "bold",
       fs.random = 10, ff.random = "bold",
       text.random.w=("")
)






                                                       #######Network Meta-Analyses#################




                                        ######Mortality#####

######Children

#Fit network consistency model 
nmac <- netmeta(TE,
                seTE,
                treat1,
                treat2,
                studlab,
                netdeath,
                subset = netdeath$age_group=="c",
                reference.group = "Quinine",
                sm = "RR", 
                comb.fixed = FALSE,
                sep.trt=" vs ",
                prediction = T)
#model output
 summary(nmac)
 
 
 #Netgraph
 child<- as.data.frame.matrix(netdeath)
 child<- subset(netdeath,age_group=="c")
 #Calculating the nodesize
 
 nQN <- sum(subset(child, t2=="QN")$n2)
 
 nASU <- sum(subset(child, t1=="ASU")$n1)
 
 nAMI <- sum(subset(child, t1=="AMI")$n1)
 nAMI <- nAMI + sum(subset(child, t2=="AMI")$n2)
 
 nAME<- sum(subset(child, t2=="AME")$n2)
 nAME <- nAME + sum(subset(child, t1=="AME")$n1)
 
 nATE<- sum(subset(child, t1=="ATE")$n1)
 
 nodesize <- c(nATE,nAME,nAMI,nASU,nQN)/sum(c(nQN,nAME,nASU,nATE,nAMI))*100
 
 #netgraph
 netgraph(nmac,
          plastic = FALSE,
          labels = c("Arteether; ATE","Artemether; AME","Artemisinin; AMI","Artesunate; ASU","Quinine; QN"),
          points = TRUE,
          cex.points = nodesize,
          offset = 0.14,
          col.points= "blue",
          col = "black",
          number.of.studies = T,
          scale = 0.95,
          col.number.of.studies = "grey",
          thickness = "number.of.studies")
 

  #League table
 lcm <- netleague(nmac, bracket = "(", digits=2)
 write.csv(lcm$random, "lcm.csv")
 
 
 ##Inconsistencies and heterogeneity###
 
 #   Investigating heterogeneity and inconsistency:
 decomp.design(nmac)
 
 #   Agreement between direct and indirect evidence using node splitting:
 netsplit(nmac)
 #   Notice you can also see contribution of the direct evidence to the
 #   estimation of the RR for each comparison of the network.
 
 #   Use a forest plot to graphically contrast the results using the
 #   direct evidence vs indirect evidence:
 forest(netsplit(nmac))
 
 #Detailed output of direct evidence vs indirect evidence
 print(netsplit(nmac), digits = 2, ci = TRUE, test = T)
 
 #Netheat Plot
 netheat(nmac,
         showall = T,
         tau.preset=T,
         random = T)
 


 
 
 ########Adults

 #Fit network consistency model 
nmaa <- netmeta(TE,
                seTE,
                treat1,
                treat2,
                studlab,
                netdeath,
                subset = netdeath$age_group=="a",
                reference.group = "Quinine",
                sm = "RR", 
                sep.trt=" vs ",
                comb.fixed = FALSE,
                prediction = F)
#display model output
summary(nmaa)

##Netgraph
adult<- as.data.frame.matrix(netdeath)
adult<- subset(netdeath,age_group=="a")
#Calculating the nodesize

nQN <- sum(subset(adult, t2=="QN")$n2)

nASU <- sum(subset(adult, t1=="ASU")$n1)


nAMI <- sum(subset(adult, t1=="AMI")$n1)
nAMI <- nAMI + sum(subset(adult, t2=="AMI")$n2)

nAME<- sum(subset(adult, t2=="AME")$n2)
nAME <- nAME + sum(subset(adult, t1=="AME")$n1)

nodesize <- c(nAME,nAMI,nASU,nQN)/sum(c(nQN,nAME,nASU,nAMI))*70

#netgraph
netgraph(nmaa,
         plastic = FALSE,
         labels = c("Artemether; AME","Artemisinin; AMI","Artesunate; ASU","Quinine; QN"),
         points = TRUE,
         cex.points = nodesize,
         col.points= "blue",
         offset = 0.09,
         col = "black",
         number.of.studies = T,
         col.number.of.studies = "grey",
         thickness = "number.of.studies",
         scale = 0.95)



# Produce League table in CSV format
(lam<- netleague(nmaa, bracket = "(", digits=2))
write.csv(lam$random, "lam.csv")


##Inconsistencies and heterogeneity###

#   Investigating heterogeneity and inconsistency:
decomp.design(nmaa)

#Detailed output of direct evidence vs indirect evidence
print(netsplit(nmaa), digits = 2, ci = TRUE, test = T)

#   Use a forest plot to graphically contrast the results using the
#   direct evidence vs indirect evidence:
forest(netsplit(nmaa))








                            ######Coma Recovery Time#########

         ###Children
child<- subset(netcoma,age_group=="c")

#Fit network consistency model 
nma2c<-netmeta(TE,
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

# Produce League table in CSV format
lccrt <- netleague(nma2c, bracket = "(", digits=2)
write.csv(lccrt$random, "lccrt.csv")


##Inconsistencies and heterogeneity###

#   Investigating heterogeneity and inconsistency:
decomp.design(nma2c)

#Detailed output of direct evidence vs indirect evidence
print(netsplit(nma2c), digits = 2, ci = TRUE, test = T)

#   Use a forest plot to graphically contrast the results using the
#   direct evidence vs indirect evidence:
forest(netsplit(nma2c))




          ###Adults


adult<- subset(netcoma,age_group=="a")
#Fit network consistency model 
nma2a<-netmeta(TE,
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
# Produce League table in CSV format
lacrt <- netleague(nma2a, bracket = "(", digits=2)
write.csv(lacrt$random, "lacrt.csv")

##Inconsistencies and heterogeneity###

#   Investigating heterogeneity and inconsistency:
decomp.design(nma2a)

#Detailed output of direct evidence vs indirect evidence
print(netsplit(nma2a), digits = 2, ci = TRUE, test = T)


#   Use a forest plot to graphically contrast the results using the
#   direct evidence vs indirect evidence:
forest(netsplit(nma2a))








                                   #######Parasite clearance time#######

         ###Children
child<- subset(netpara,age_group=="c")
#Fit network consistency model 
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

# Produce League table in CSV format
lcpct <- netleague(nma3c, bracket = "(", digits=2)
write.csv(lcpctt$random, "lacrt.csv")

##Inconsistencies and heterogeneity###

#   Investigating heterogeneity and inconsistency
decomp.design(nma3c)

#Detailed output of direct evidence vs indirect evidence
print(netsplit(nma3c), digits = 2, ci = TRUE, test = T)


#   Use a forest plot to graphically contrast the results using the
#   direct evidence vs indirect evidence:
nsc<- netsplit(nma3c)
forest(nsc,
       prediction = F,
       overall = F,
       colgap.left="18mm",
       leftcols = c("studlab","k"))

#netheat plot
netheat(nma3c,
        showall = T,
        tau.preset=T,
        random = T)



         ###Adults
adult<- subset(netpara,age_group=="a")

#Fit network consistency model 
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

# Produce League table in CSV format
lapct <- netleague(nma3a, bracket = "(", digits=2)
write.csv(lapct$random, "lapct.csv")

##Inconsistencies and heterogeneity###

#   Investigating heterogeneity and inconsistency:

#Detailed output of direct evidence vs indirect evidence
print(netsplit(nma3a), digits = 2, ci = TRUE, test = T)

help(decomp.design)


#   Use a forest plot to graphically contrast the results using the
#   direct evidence vs indirect evidence:
nsa<- netsplit(nma3a)
forest(nsa,
       prediction = F,
       overall = F,
       colgap.left="18mm",
       leftcols = c("studlab","k"))

#heatplot
netheat(nma3a,
        showall = F,
        random = T)


                                                      





                                        ####Neurological Sequela Events#######

##Children
#Fit network consistency model 
nma4c <- netmeta(TE,
                 seTE,
                 treat1,
                 treat2,
                 studlab,
                 netneuro,
                 subset = netneuro$age_group=="c",
                 sm = "RR",
                 comb.fixed = FALSE,
                 prediction=T,
                 reference.group = "Quinine",
                 sep.trt=" vs ",
                 details.chkmultiarm = T)


# Produce League table in CSV format
lcneuro <- netleague(nma4c, bracket = "(", digits=2)
write.csv(lcneuro$random, "lcneuro.csv")


##Inconsistencies and heterogeneity###

#   Investigating heterogeneity and inconsistency:
decomp.design(nma4c)

#Detailed output of direct evidence vs indirect evidence
print(netsplit(nma4c), digits = 2, ci = TRUE, test = T)

#   Use a forest plot to graphically contrast the results using the
#   direct evidence vs indirect evidence:
forest(netsplit(nma4c))




##Adults
#Fit network consistency model 
nma4a<- netmeta(TE,
                seTE,
                treat1,
                treat2,
                studlab,
                netneuro,
                subset = netneuro$age_group=="a",
                sm = "RR",
                comb.fixed = FALSE,
                prediction=T,
                reference.group = "Quinine",
                sep.trt=" vs ",
                details.chkmultiarm = T)


# Produce League table in CSV format
laneuro <- netleague(nma4a, bracket = "(", digits=2)
write.csv(laneuro$random, "laneuro.csv")

##Inconsistencies and heterogeneity###

#   Investigating heterogeneity and inconsistency:
decomp.design(nma4a)

#Detailed output of direct evidence vs indirect evidence
print(netsplit(nma4a), digits = 2, ci = TRUE, test = T)

#   Use a forest plot to graphically contrast the results using the
#   direct evidence vs indirect evidence:
forest(netsplit(nma4a))




                                               #######Ranking and Forest plots ########



                          ##Mortality

##Children
#Rank
netrank(nmac, small.values = "good")

#forest plot with P-Scores and ranking

forest(nmac,
       xlim = c(0.04,4),
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
       fs.study.labels= 13, ff.study.labels = "bold",
       fs.axis = 9,
       fs.addline = 9, ff.addline = "italic",
       text.addline1=" *I^2=0%, p=0.93"
)



##Adult
#Rank
netrank(nmaa, small.values = "good")


#forest plot with P-Scores and ranking
forest(nmaa,
       xlim = c(0.25,4),
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
       text.addline1="*I^2=23.1%, p=0.24"
)




     
                                ##Coma Recovery Time

##Children
#Rank
netrank(nma2c, small.values = "good")

#forest plot with P-Scores and ranking
forest(nma2c,
       xlim = c(-86,140),
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
       text.addline1=" *I^2=61.7%, p=0.005"
)

##Adult
#Rank
netrank(nma2a, small.values = "good")

#forest plot with P-Scores and ranking
forest(nma2a,
       xlim = c(-40,40),
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
       text.addline1=" *I^2=0%, p=0.54"
)




                          ##Parasite clearance time

##Children
#Rank
netrank(nma3c, small.values = "good")

#forest plot with P-Scores and ranking
forest(nma3c,
       xlim = c(-50,13),
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
       text.addline1=" *I^2=77.6%, p<0.0001"
)

##Adult
#Rank
netrank(nma3a, small.values = "good")

#forest plot with P-Scores and ranking
forest(nma3a,
       xlim = c(-50,10),
       smlab = "Parasite Clearance Time",
       small.values = "good",
       reference.group = "Quinine",
       just.smlab = "left",
       leftcols =c("studlab","Pscore"),
       sortvar = -Pscore,
       col.square= "blue",
       col.square.lines = "blue",
       fontsize=12,
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
       text.addline1=" *I^2=89.7%, p<0.0001"
)





                          ##Neurological Sequela Events



##Children

#Rank
netrank(nma4c, small.values = "good")

#forest plot with P-Scores and ranking
forest(nma4c,
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
       text.addline1=" *I^2=24.4%, p=0.23"
)



##Adult

#Rank
netrank(nma4a, small.values = "good")

#forest plot with P-Scores and ranking

forest(nma4a,
       xlim = c(0.2,15),
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
       text.addline1=" *I^2=0%, p=0.57"
)



                      ###Hasse diagram


#Children

## P-scores 
pscore.deathc  <-c(0.8517, 0.2472, 0.5986, 0.6468, 0.1557)
pscore.neuroc <- c(0.3740, 0.6469, 0.4617, 0.3551, 0.6596) 
pscore.pctc  <-  c(0.9945, 0.5610, 0.1330, 0.6646, 0.1469) 
pscore.crtc  <-  c(0.7048, 0.7048, 0.3893, 0.1148, 0.5071) 

pscore.matrixc<-data.frame(pscore.deathc, pscore.crtc,pscore.pctc)


# Partial order of treatment rankings of mortality, neuro sequela, crt and pct
netc<-netposet(pscore.matrixc)

#Construct a logical matrix

netsetc<- matrix(c(rep(0,5),1,rep(0,4),1,rep(0,4),1,rep(0,5),1,0,0,0),5,5)
netsetc

netsetc<-as.logical(netsetc)
netc<-matrix(netsetc,5,5)
netc

#Construct Hasse Diagram
hasse(netc, c("Artemisinin", "Artemether", "Arteether", "Artesunate","Quinine"))





#Adults


# P-scores 
pscore.deatha  <-c(0.2452, 0.7554, 0.8312, 0.1682)
pscore.neuroa <- c(0.000, 0.2431, 0.000, 0.7569) 
pscore.crta  <-  c(0.5305, 0.2610, 0.8889, 0.3196) 
pscore.pcta  <-  c(0.7302, 0.4247,0.8111, 0.0341) 


# Construct matrix with P-scores
pscore.matrixa<-data.frame(pscore.deatha, pscore.crta,pscore.pcta)



# Partial order of treatment rankings of mortality, neuro sequela, crt and pct
neta<-netposet(pscore.matrixa)
neta


#Convert to a logical matrix
netseta<- matrix(c(rep(0,2),1,rep(0,3),1,rep(0,5),1,rep(0,3)),4,4)

netseta<-as.logical(netseta)
neta<-matrix(netseta,4,4)

#Construct Hasse Diagram

hasse(neta, c("Artemisinin", "Artemether", "Artesunate","Quinine"))






                      #####Adverse Events####

####Network meta-analyses of hypoglycaemia events

#Fit network consistency model 
nma5 <- netmeta(TE,
                seTE,
                treat1,
                treat2,
                studlab,
                nethypo,
                sm = "RR",
                comb.fixed = FALSE,
                prediction=T,
                reference.group = "Quinine",
                sep.trt=" vs ",
                details.chkmultiarm = T)

#View output of model
summary(nma5)

#Forest plot of hypoglycaemia events
forest(nma5,
       xlim = c(0.04,4),
       smlab = "Hypoglycaemia Events",
       small.values = "good",
       reference.group = "Quinine",
       just.smlab = "left",
       col.square.lines = "blue",
       fontsize =10,
       plotwidth = "8cm",
       colgap.left="4mm",
       spacing= 2,
       print.I2= F,
       print.pval.Q=F,
       hetlab= " ",
       label.left = "favours treatment",
       label.right= "favours quinine",
       fs.smlab = 11, ff.smlab = "bold",
       fs.study.labels= 10, ff.study.labels = "bold",
       fs.axis = 8,
       fs.addline = 8, ff.addline = "italic",
       text.addline1=" Variability: I^2=0%, p=0.45"
)


##Inconsistencies and heterogeneity##

#   Investigating heterogeneity and inconsistency:
decomp.design(nma5)

#Detailed output of direct evidence vs indirect evidence
print(netsplit(nma5), digits = 2, ci = TRUE, test = T)

#   Use a forest plot to graphically contrast the results using the
#   direct evidence vs indirect evidence:
forest(netsplit(nma5))


#Ranking of hypoglycaemia events

netrank(nma5, small.values = "good")




##Analyses of EcG abnormalities 

#Transform data from two arm-based formats into contrast-based format
netecg <- pairwise(treat = T,
                   event = ecg_e,
                   n = ecg_n,
                   data = data,
                   allstudies = T,
                   studlab = author_year,
                   sm = "RR")


#Conduct traditional pairwise meta-analyses with "metagen"
maecg<- metagen(TE,
               seTE,
               studlab,
               data = netecg,
               sm = "RR", 
               comb.fixed = FALSE,
               prediction = T,
               method.tau= "REML",
               method.tau.ci = "QP")


#Construct forest plot

forest(maecg, 
       xlim = c(0.01, 70), 
       test.subgroup.random=F,
       test.effect.subgroup.random=F,
       label.test.effect.subgroup.random = "Test for effect:  " ,
       smlab = "ECG Abnormalities",
       small.values = "good",
       rightcols = c("effect","ci"),
       leftcols = c("studlab","event1","n1","event2","n2"),
      
       lab.e = "Artemether",
       lab.c = "Quinine",
       lab.e.attach.to.col = "event1",
       lab.c.attach.to.col = "event2",
       just.smlab = "left",
       leftlabs = c("Study","events","out of","events", "out of"),
       title= "ECG",
       
       col.square.lines = "blue",
       fontsize =12,
       plotwidth = "8cm",
       colgap.left="4mm",
       spacing= 3,
       overall = T,
       overall.hetstat = T,
       text.subgroup.nohet="NA",
       prediction = F,
       studlab = T,
       study.results = T,
       pooled.totals = T,
       print.tau2 = F,
       label.left = "favours artemether",
       label.right= "favours quinine",
       fs.smlab = 13, ff.smlab = "bold",
       fs.study.labels= 11, ff.study.labels = "bold",
       fs.axis = 9,
       fs.addline = 7, ff.addline = "bold.italic"
)



help("forest.meta")

                                  #############Sensitivity Analyses#######################


                #######Traditional Paiwise Meta_Analyses using "metagen" with DL and Jacksons######

##Mortality
madeath<- metagen(TE,
                  seTE,
                  studlab,
                  data = netdeath,
                  byvar = netdeath$Pairwise, 
                  sm = "RR", 
                  comb.fixed = FALSE,
                  prediction = T,
                  )

#display model output
summary(madeath)




  
                ######Lumped Population

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

#display model output
summary(nma1)


##Netgraph
####Calculating the nodesize


nQN <- sum(subset(netdeath, t2=="QN")$n2)

nASU <- sum(subset(netdeath, t1=="ASU")$n1)

nAMI <- sum(subset(netdeath, t1=="AMI")$n1)
nAMI <- nAMI + sum(subset(netdeath, t2=="AMI")$n2)

nAME<- sum(subset(netdeath, t2=="AME")$n2)
nAME <- nAME + sum(subset(netdeath, t1=="AME")$n1)

nATE<- sum(subset(netdeath, t1=="ATE")$n1)

nodesize <- c(nATE,nAME,nAMI,nASU,nQN)/sum(c(nQN,nAME,nASU,nATE,nAMI))*100



#netgraph
netgraph(nma1,
         plastic = FALSE,
         seq = c("Artemisinin","Artemether","Artesunate","Arteether","Quinine"),
         points = TRUE,
         offset = 0.136,
         col.points = "blue",
         col = "black",
         scale = 0.95,
         number.of.studies = T,
         col.number.of.studies = "grey",
         thickness = "number.of.studies",
         cex.points = nodesize,
         labels = c("Arteether; ATE" ,"Artemether; AME","Artemisinin; AMI","Artesunate; ASU","Quinine; QN"))



# Produce League table in CSV format
leaguedeath <- netleague(nma1, bracket = "(", digits=2)
write.csv(leaguedeath$random, "leaguedeath.csv")

##Inconsistencies and heterogeneity###

#   Investigating heterogeneity and inconsistency:
decomp.design(nma1)

#Detailed output of direct evidence vs indirect evidence
print(netsplit(nma1), digits = 2, ci = TRUE, test = T)

#   Use a forest plot to graphically contrast the results using the
#   direct evidence vs indirect evidence:
forest(netsplit(nma1))


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


# Produce League table in CSV format
leaguecoma <- netleague(nma2, bracket = "(", digits=2)
write.csv(leaguecoma$random, "leaguecoma2.csv")


##Inconsistencies and heterogeneity###

#   Investigating heterogeneity and inconsistency:
decomp.design(nma2)

#Detailed output of direct evidence vs indirect evidence
print(netsplit(nma2), digits = 2, ci = TRUE, test = T)

#   Use a forest plot to graphically contrast the results using the
#   direct evidence vs indirect evidence:
forest(netsplit(nma2))



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

# Produce League table in CSV format
leaguepara <- netleague(nma3, bracket = "(", digits=2)
write.csv(leaguepara$random, "leaguepara2.csv")


##Inconsistencies and heterogeneity###

#   Investigating heterogeneity and inconsistency:
decomp.design(nma3)

#Detailed output of direct evidence vs indirect evidence
print(netsplit(nma3), digits = 2, ci = TRUE, test = T)


#   Use a forest plot to graphically contrast the results using the
#   direct evidence vs indirect evidence:
ns<- netsplit(nma3)
forest(ns,
       prediction = F,
       overall = F,
       colgap.left="18mm",
       leftcols = c("studlab","k"))

#heatplot
netheat(nma3,
        showall = T,
        random = T)


nma3

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

# Produce League table in CSV format
leagueneuro <- netleague(nma4, bracket = "(", digits=2)
write.csv(leagueneuro$random, "leagueneuro.csv")

##Inconsistencies and heterogeneity###

#   Investigating heterogeneity and inconsistency:
decomp.design(nma4a)

#Detailed output of direct evidence vs indirect evidence
print(netsplit(nma4a), digits = 2, ci = TRUE, test = T)

#   Use a forest plot to graphically contrast the results using the
#   direct evidence vs indirect evidence:
forest(netsplit(nma4a))




#Ranking and Forest plots #

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
       text.addline1=" *I^2=87.7%, p<0.001"
)

###Neurological sequela events####
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


###Hasse diagram###

#P-scores 
pscore.death  <-c(0.3145, 0.5060, 0.7045, 0.8058, 0.1693)
pscore.neuro <- c(0.3721, 0.6577, 0.4599, 0.3529, 0.6573) 
pscore.crt  <-  c(0.8412, 0.5852, 0.1576, 0.6323, 0.2837) 
pscore.pct  <-  c(0.9417, 0.4824, 0.1425, 0.7964, 0.1370) 



#Death, Parasite Clerance Time And Coma Recovery Time


# Construct matrix with P-scores
pscore.matrix<-data.frame(pscore.death, pscore.crt, pscore.pct)


# Partial order of treatment rankings of mortality, neuro sequela, crt and pct
net<-netposet(pscore.matrix)


#Convert to a logical matrix
netset<- matrix(c(rep(0,8),1,rep(0,4),1,rep(0,6),1,1,0,0,0),5,5)

netset<-as.logical(netset)
net<-matrix(netset,5,5)

#Construct a Hasse daigram
hasse(net, labels = c("Artemisinin", "Artemether", "Arteether", "Artesunate","Quinine"))







     ##########Hasse diagram for all Outcomes (mortality, neuro sequela, crt, pct and hypoglycaemia)#####

#P-scores 
pscore.death  <-c(0.3145, 0.5060, 0.7045, 0.8058, 0.1693)
pscore.neuro <- c(0.3721, 0.6577, 0.4599, 0.3529, 0.6573) 
pscore.crt  <-  c(0.8412, 0.5852, 0.1576, 0.6323, 0.2837) 
pscore.pct  <-  c(0.9417, 0.4824, 0.1425, 0.7964, 0.1370) 
pscores.hypo<-  c(0.8518, 0.6821, 0.0000, 0.4589, 0.0071)


#Construct matrix with P-scores

pscore.matrix<-data.frame(pscore.death,pscore.neuro, pscore.crt, pscore.pct,pscores.hypo)

# Partial order of treatment rankings of mortality, neuro sequela, crt, pct and hypoglycaemia
net<-netposet(pscore.matrix)
net



#Convert to a logical matrix
netset<- matrix(c(rep(0,21),1,rep(0,3)),5,5)

netset

netset<-as.logical(netset)
net<-matrix(netset,5,5)
net


 #Construct a Hasse daigram
hasse(net, labels = c("Artemisinin", "Artemether", "Arteether","Artesunate","Quinine"))



########################### Anlayses with extra data from systematic reviews ################################



#Fit a NMA model for dataset that includes extra data
netdeathx <- pairwise(treat = T,
                      event = mortality_e,
                      n = mortality_n,
                      data = datax,
                      allstudies = T,
                      studlab = author_year,
                      sm = "RR")



nma7 <- netmeta(TE,
                seTE,
                treat1,
                treat2,
                studlab,
                netdeathx,
                reference.group = "Quinine",
                sm = "RR", 
                comb.fixed = FALSE,
                prediction = T)

#display model output
summary(nma7)



#Compare and display results with a forest plot
nb2 <- netbind(nma1, nma7,
               name = c(" Available RCTs only ", " Available RCTs with added RCTs"),
               col.study = c("red", "black"),
               col.square = c("red", "black"))
#forest plot
forest(nb2,
       xlim = c(0.3, 2),
       col.by="black",
       addrow.subgroups = T,
       plotwidth = "8cm",
       colgap.left="6mm",
       colgap.right="4mm",
       leftlabs = "Treatment",
       smlab = "Mortality",
       fontsize = 13, 
       spacing = 1.5, 
       squaresize = 0.7,
       label.left = "favours treatment",
       label.right = "favours quinine",
       fs.smlab=14)





############# Mean SD vs Median  iqr/range#################################



#Fit a NMA model for dataset that includes only reported means and sds

netcrt<- pairwise(treat=T,
                    n = crt_n,
                    mean=crt_mean,
                    sd=crt_sd,
                    data = datax,
                    studlab = author_year,
                    sm = "MD",
                    allstudies =T )



nma8<-netmeta(TE,
                seTE,
                treat1,
                treat2,
                studlab,
                netcrt,
                sm = "MD",
                comb.fixed = FALSE,
                prediction=T,
                reference.group = "Quinine",
                sep.trt="vs",
                details.chkmultiarm = T)

#Compare and display results with a forest plot

nb1 <- netbind(nma8,nma2,
               name = c(" Reported means and sds only", " Reported means and sds with approximations"),
               col.study = c("red", "black"),
               col.square = c("red", "black"))

#forest plot
forest(nb1, 
       col.by="black",
       addrow.subgroups = T,
       plotwidth = "8cm",
       colgap.left="4mm",
       colgap.right="4mm",
       leftlabs = "Treatment",
       smlab = "Coma Recovery Time",
       fontsize = 13, 
       spacing = 1.5, 
       squaresize = 0.7,
       label.left = "favours treatment",
       label.right = "favours quinine",
       fs.smlab=15)



######## Parasite Clearance Time #######


#Fit a NMA model for dataset that includes only reported means and sds

netpct<- pairwise(treat=T,
                    n = pct_n,
                    mean=pct_mean,
                    sd=pct_sd,
                    data = datax,
                    studlab = author_year,
                    sm = "MD",
                    allstudies =T )


nma9<-netmeta(TE,
               seTE,
               treat1,
               treat2,
               studlab,
               netpct,
               sm = "MD",
               comb.fixed = FALSE,
               prediction=T,
               reference.group = "Quinine",
               sep.trt="vs",
               details.chkmultiarm = T)


#Compare and display results with a forest plot

nb2 <- netbind(nma9, nma3,
               name = c(" Reported means and sds only ", " Reported means and sds with approximations"),
               col.study = c("red", "black"),
               col.square = c("red", "black"))


#forest plot
forest(nb2, 
       col.by="black",
       addrow.subgroups = T,
       plotwidth = "8cm",
       colgap.left="4mm",
       colgap.right="4mm",
       leftlabs = "Treatment",
       smlab = "Parasite Clearance Time",
       fontsize = 13, 
       spacing = 1.5, 
       squaresize = 0.7,
       label.left = "favours treatment",
       label.right = "favours quinine",
       fs.smlab=15)



              #########SUBGROUP ANALYSES######



 ###Cerebral Malaria vs non specified severe malria

#
netcereb<- pairwise(treat = T,
                    event = cereb_e,
                    n = cereb_n,
                    data = datasub,
                    allstudies = T,
                    studlab = author_year,
                    sm = "RR")



nma10 <- netmeta(TE,
                seTE,
                treat1,
                treat2,
                studlab,
                netcereb,
                subset = netcereb$c_malaria.only =="Yes",
                reference.group = "Quinine",
                sm = "RR", 
                comb.fixed = FALSE,
                prediction = T)

# Produce League table in CSV format
leagueyes<- netleague(nma10, bracket = "(", digits=2)
write.csv(leagueyes$random, "leagueayes.csv")

#display model output
summary(nma10)

#   Investigating heterogeneity and inconsistency:
decomp.design(nmay)

#forest plot
forest(nmay, 
       reference.group = "QN",
       col.square = "blue")




nma11 <- netmeta(TE,
                seTE,
                treat1,
                treat2,
                studlab,
                netcereb,
                subset = netcereb$c_malaria.only =="No",
                reference.group = "Quinine",
                sm = "RR", 
                comb.fixed = FALSE,
                prediction = T)

# Produce League table in CSV format
leagueno <- netleague(nma11, bracket = "(", digits=2)
write.csv(leagueno$random, "leagueno.csv")

#display model output
summary(nma11)


#   Investigating heterogeneity and inconsistency:
decomp.design(nman)

#forest plot
forest(nman, 
       reference.group = "QN",
       col.square = "blue")




########Study Continent

nmAf <- netmeta(TE,
                seTE,
                treat1,
                treat2,
                studlab,
                netdeath,
                subset = netdeath$study.continent =="Africa",
                reference.group = "Quinine",
                sm = "RR", 
                comb.fixed = FALSE,
                prediction = T)
#display model output
summary(nmAf)


# Produce League table in CSV format
leagueAf<- netleague(nmAf, bracket = "(", digits=2)
write.csv(leagueAf$random, "leagueAf.csv")


#forest plot
forest(nmAf, 
       reference.group = "QN",
       col.square = "blue")



nmAs <- netmeta(TE,
                seTE,
                treat1,
                treat2,
                studlab,
                netdeath,
                subset = netdeath$study.continent =="Asia",
                reference.group = "Quinine",
                sm = "RR", 
                comb.fixed = FALSE,
                prediction = T)

#display model output
summary(nmAs)

# Produce League table in CSV format

leagueAs <- netleague(nmAs, bracket = "(", digits=2)
write.csv(leagueAs$random, "leagueAs.csv")


#forest plot
forest(nmAs, 
       reference.group = "QN",
       col.square = "blue")




#####Neurological sequela events

netneurolo<- pairwise(treat = T,
                      event = neuro_e,
                      n = neuro_n,
                      data = datasub,
                      allstudies = T,
                      studlab = author_year,
                      sm = "RR")


nmacute <- netmeta(TE,
                   seTE,
                   treat1,
                   treat2,
                   studlab,
                   netneurolo,
                   subset = netneurolo$neuro =="a",
                   reference.group = "Quinine",
                   sm = "RR", 
                   comb.fixed = FALSE,
                   prediction = T)

#display model output
summary(nmacute)

# Produce League table in CSV format
(leagueacute<- netleague(nmacute, bracket = "(", digits=2))
write.csv(leagueacute$random, "leagueacute.csv")

#   Investigating heterogeneity and inconsistency:
decomp.design(nmacute)

#forest plot
forest(nmacute, 
       reference.group = "QN",
       col.square = "blue")





nmaper <- netmeta(TE,
                  seTE,
                  treat1,
                  treat2,
                  studlab,
                  netneurolo,
                  subset = netneurolo$neuro =="p",
                  reference.group = "Quinine",
                  sm = "RR", 
                  comb.fixed = FALSE,
                  prediction = T)



# Produce League table in CSV format
leagueper<- netleague(nmaper, bracket = "(", digits=2)
write.csv(leagueper$random, "leagueaper.csv")


summary(nmaper)

#   Investigating heterogeneity and inconsistency:
decomp.design(nmaper)

#forest plot
forest(nmaper, 
       reference.group = "QN",
       col.square = "blue")




########Comparison Adjusted Funnel PLot 

#create order for funnel lot
ord<- c('AMI',"AME","ATE","ASU","QN")



#funnel plot 
oldpar <- par(mfrow = c(2, 1))
funnel.netmeta(nma1, 
               order = ord,
               fontsize=12,
               linreg = F,
               spacing=1.2,
               legend=F,
               lwd=2, cex = 1, pch = 16, cex.studlab = 1.25,
               rank=F)


#Funnel plot with comparisons 
funnel.netmeta(nma1, 
               order = ord,
               linreg = T,
               spacing=1.2,
               fontsize=12,
               lty.ref.triangle=3,
               col = c("Red","Blue","Green","Yellow","Orange","Brown", "Pink"))









#Children


#Adults

##Mortality
##Coma Recovery Time
##Parasite clearance time
##Neurological Sequela Events