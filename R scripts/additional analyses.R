
#### Network Analyses descriptives###
#How I got the number of participants and confirmed numbers of pairwise comparisons, studies and treatments.
#Using Coma Recovery time analyses as example

#convert output "pairwise" command to dataframe

netcoma<-as.data.frame.matrix(netcoma)


#Number of studiee
length( unique(netcoma$studlab))

#take studies omitted from the model out ; Vihn et al 1997

netcoma<-slice(netcoma,-c(9,10,11))


#create factor variable with each pairwise comaprison

netcoma<- mutate(netcoma, Pairwise= derivedFactor (
  "ASN vs QN" =(t1=="ASN"&t2=="QN"),
  "ATM vs QN" =(t1=="ATM"&t2=="QN"),
  "ATT vs QN" =(t1=="ATT"& t2=="QN"),
  "ATS vs QN" =(t1=="ATS"& t2=="QN"),
  "ATS vs ATM" =(t1=="ATS"& t2=="ATM"),
  "ATS vs ASN" =(t1=="ATS"& t2=="ASN"),
  "ASN vs ATM" =(t1=="ASN"& t2=="ATM"),
  .method="first",
  .default=0))

#Number of pairwise comparisons
length(netcoma$Pairwise)

#Number of participants
         

sum(unique(data$N))

#Number of children

unique(r<- data%>%
         select(age_group,author_year,N,c_malaria.only))

sum(unique(r$N))



child<- data%>%
         filter(age_group=="a")

sum(unique(child$N))

N<-sum(data$n1)
N<-N + sum(data$n2)

View(data)

unique(r<- data%>%
         filter(age_group=="c")%>%
         select(age_group,N,author_year,))
length (unique(r$author_year))

nc<- sum(unique(subset(r, age_group=="c")$N))
nc


unique
(r<- data%>%
         filter(age_group=="a")%>%
         select(age_group,N,author_year))
length (r$author_year)

na <- sum(unique(subset(r, age_group=="a")$N))
na


unique(r<- data%>%
         select(age_group,N,author_year))
length (unique(r$author_year))

N <- sum(unique(r$N))
N

unique(r<- desc%>%
         select(age_group,author_year,N,c_malaria.only))
length (unique(r$author_year))

r<-slice(r,-c(42,43))
desc

N <- sum(unique(r$N))
N


unique(r<- desc%>%
         filter(age_group=="c")%>%
         select(age_group,author_year,N,c_malaria.only))
length (unique(r$author_year))

nc<- sum(unique(subset(r, age_group=="c")$N))
nc


unique(r<- desc%>%
         filter(age_group=="a")%>%
         select(age_group,author_year,N,c_malaria.only))
length (unique(r$author_year))

na <- sum(unique(subset(r, age_group=="a")$N))
na
