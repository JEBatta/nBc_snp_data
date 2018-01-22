####Borrar todas las variables y grÃ¡ficos previos#####
rm(list = setdiff(ls(), lsf.str()))
graphics.off()

####################################################
####################################################

groups <- function(x,n){ #divide x in n categories from 1 to n dividing range almost equally 
  d<-length(unique(x))
  if(d < n){n <- d}
  d<-ceiling(d/n)
  a<-min(x)
  b<-a+d
  res<-x
  for(i in 1:(n-1)){
    res[which(x >= a & x < b)]<-i
    a<- b
    b<-b+d
  }
  res[which(x >= a)]<-n
  return(res)
}

groupDatabase<-function(db,n){
  res<-db
  for(i in 1:length(db)){
    res[,i]<-groups(db[,i],n)
  }
  return(res)
}
####################################################

setwd("C:/Users/Jesús/Desktop/DataMining-SNPs/6_110118_pathways")
md<-read.csv('path_database_notZero.csv')
md$X<-NULL
n<-nrow(md)
path<-md[,352:358]
pav<-sapply(1:length(path), function(x) mean(path[,x]))
psd<-sapply(1:length(path), function(x) sd(path[,x]))

#sapply(1:length(path), function(x) sum(path[,x]<(pav[x]-psd[x])))
#sapply(1:length(path), function(x) sum(path[,x]>(pav[x]+psd[x])))
#sapply(1:length(path), function(x) sum(path[,x]<=(pav[x]+psd[x]) & path[,x]>=(pav[x]-psd[x])))
pathA<-path
pathB<-path
for(i in 1:length(path)){
  idxA1<-which(path[,i]<(pav[i]-psd[i]))
  idxA2<-which(path[,i]<=(pav[i]+psd[i]) & path[,i]>=(pav[i]-psd[i]))
  idxA3<-which(path[,i]>(pav[i]+psd[i]))
  idxB0<-which(path[,i]<pav[i])
  idxB1<-which(path[,i]>=pav[i])
  pathA[idxA1,i]<-1
  pathA[idxA2,i]<-2
  pathA[idxA3,i]<-3
  pathB[idxB0,i]<-0
  pathB[idxB1,i]<-1
}

clinic<-md[,1:4]
SNP<-md[,5:351]

write.csv(cbind(clinic,pathA,SNP),'pathA.csv')
write.csv(cbind(clinic,pathB,SNP),'pathB.csv')

#############################################################################
## 
path3g<-groupDatabase(path,3)
write.csv(cbind(clinic,path3g,SNP),'path3g.csv')
path4g<-groupDatabase(path,4)
write.csv(cbind(clinic,path4g,SNP),'path4g.csv')
path5g<-groupDatabase(path,5)
write.csv(cbind(clinic,path5g,SNP),'path5g.csv')
