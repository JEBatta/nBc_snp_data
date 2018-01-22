####Borrar todas las variables y grÃ¡ficos previos#####
rm(list = setdiff(ls(), lsf.str()))
graphics.off()

####################################################
####################################################

setwd("C:/Users/Jesús/Desktop/DataMining-SNPs/6_110118_pathways")
#md<-read.csv('path_database_notZero.csv')
md<-read.csv('pathB.csv')
#md<-read.csv('path3g.csv')
#md<-read.csv('path4g.csv')
#md<-read.csv('path5g.csv')

md$X<-NULL
########################################################
#only paths (this part is only for verify if path alone are 
# predictive. Don't use otherwise)
#md<-md[,1:9]
#########################################################

n<-nrow(md)
K<-10
MCK<-10
set.seed(123) # seed is fixed for reproducibility purposes
source("sourcesNBC_gen2.R")
lpaths<-7
df <- data.frame(feature=character(),
                 category=numeric(), 
                 epsilon=numeric(),
                 score=numeric(),
                 nx=numeric(),
                 nxc=numeric(),
                 n=numeric(),
                 nc=numeric(),
                 assemble=numeric(),
                 stringsAsFactors=FALSE)
AUC<-data.frame(ob=rep(0,K*MCK),ao=rep(0,K*MCK),tgb=rep(0,K*MCK))


for(i in 2:4){

for(jj in 1:MCK){
  
  idx<-split(1:n, sample(ceiling(seq_along(1:n)/(n/K))))
  
  for(ii in 1:K){
    
    idx_test<-idx[[ii]]
    idx_train<-which(!is.element(1:n,idx_test))
    md_train<-md[idx_train,]
    md_test<-md[idx_test,] 
    n_train<-nrow(md_train)
    n_test<-nrow(md_test)
    
    X<-md_train[,5:length(md_train)]
    
    #i<-2 #index: 2. BMI >= 30, 3. WC High, 4. TGB high
    
    ncs<-sum(md_train[i] == 1, na.rm = T)
    
    TPR<-rep(0,n_test)
    FPR<-rep(0,n_test)
    
    nam<-names(md_train)[i]
    y<-md_train[,i]
    v<-1
    lst <- sapply(1:length(X), function(x) funcEpsilon(X,x,y,v))
    mxX<-max(sapply(lst,length))
    eps <- as.data.frame(do.call(rbind,lapply(lst, `length<-`,
                                              mxX)))
    
    lst <- sapply(1:length(X), function(x) funcNCX(X,x,y,v))
    nxc <- as.data.frame(do.call(rbind,lapply(lst, `length<-`,
                                              mxX)))
    
    lst <- sapply(1:length(X), function(x) funcNX(X,x,y,v))
    nx <- as.data.frame(do.call(rbind,lapply(lst, `length<-`,
                                             mxX)))
    
    lst <- sapply(1:length(X), function(x) funcScore(X,x,y,v))
    score <- as.data.frame(do.call(rbind,lapply(lst, `length<-`,
                                                mxX)))
    
    relevantFeatures<-funResultsEps1(X,eps,score,nx,nxc)
    
    nn<-nrow(relevantFeatures)
    
    rft<-data.frame(relevantFeatures, n = rep(n,nn), nc = rep(ncs,nn), assemble = rep(ii,nn))
    df<-merge(df,rft,all = T)
    
    Pc<-sum(y,na.rm = T)/n_train
    iSc<-sapply(1:n_test, function(x) funIndivScore(md_test,x,relevantFeatures,Pc))
    iRes<-md_test[,i]
    ORes<-iRes[order(iSc,decreasing = T)]
    for(j in 1:n_test){
      TP<-sum(ORes[1:j]==1, na.rm = T)
      FN<-sum(ORes[(j+1):n_test]==1, na.rm = T)
      FP<-sum(ORes[1:j]==0, na.rm = T)
      TN<-sum(ORes[(j+1):n_test]==0, na.rm = T)
      TPR[j]<-TP/(TP+FN)
      FPR[j]<-FP/(FP+TN)
    }
    temp<-0
    for (j in 1:(n_test-1)){
      AUC[(jj-1)*10 + ii,i-1]<-AUC[(jj-1)*10 + ii,i-1]+TPR[j]*(FPR[j]-temp)
      temp<-FPR[j]
    }  
  }
}
}
#write.csv(df,'relevant_SNPs_HighTGB_all_list.csv')
write.csv(AUC,'auc_pathB.csv')
library(plotly)
plot_ly(y = AUC$ob, type = "box")%>%
  add_boxplot(y = AUC$ao)%>%
  add_boxplot(y = AUC$tgb)
####################################################
## get frequencies of relevant SNPs (in K-fold lists)
## together with average score and epsilon (and StdDev)

freqSNP<-data.frame(feature=character(),
                    category=numeric(), 
                    AVepsilon=numeric(),
                    SDepsilon=numeric(),
                    AVscore=numeric(),
                    SDepsilon=numeric(),
                    freq=numeric(),
                    stringsAsFactors=FALSE)
for( ft in unique(df$feature)){
  ts1<-subset(df, feature == ft)
  for(ct in unique(ts1$category)){
    ts2<-subset(ts1, category == ct)
    freqSNP<-rbind(freqSNP,data.frame(feature=ft,
                                      category=ct,
                                      AVepsilon=mean(ts2$epsilon),
                                      SDepsilon=sd(ts2$epsilon),
                                      AVscore=mean(ts2$score),
                                      SDscore=sd(ts2$score),
                                      freq=nrow(ts2)))
  }
}
