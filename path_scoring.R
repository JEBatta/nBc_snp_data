####Borrar todas las variables y grÃ¡ficos previos#####
rm(list = setdiff(ls(), lsf.str()))
graphics.off()

####################################################
####################################################

setwd("C:/Users/Jesús/Desktop/DataMining-SNPs/6_110118_pathways") #set working directory
paths<-read.csv('pathways_dictionary.csv')
md<-read.csv('database_snp.csv')
md$X<-NULL
n<-nrow(md)
for(i in 5:length(md)){
  names(md)[i]<-strsplit(names(md)[i],"_")[[1]][1]   
}
paths$rsIDs<-gsub("\t","",paths$rsIDs)
paths$rsIDs<-gsub(" ","",paths$rsIDs)

md$path2<-rep(0,n)
md$path4<-rep(0,n)
md$path5<-rep(0,n)
md$path6<-rep(0,n)
md$path7<-rep(0,n)
md$path9<-rep(0,n)
md$path11<-rep(0,n)
pval<-c(2,4,5,6,7,9,11)
idx<-list()
for(i in 1:7){
  idx[[i]] <-which(names(md) %in% unique(paths$rsIDs[paths$path==pval[i]]))
}

for(i in 1:n){
  for(j in 1:7){
    md[i,j+351]<-sum(md[i,idx[[j]]]>0,na.rm = T)
  }
}
write.csv(md,'path_database_notZero.csv')

