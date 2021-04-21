library(MASS)
library(parallel)

# geno_n is genotype vector, including 0, 1 and 2
# y_n is the phenotype matrix, each row is a subject, each column represents a phenotype
geno_2n=do.call('c',lapply(geno_n,function(x) if(x==2){c(1,1)} else if(x==0) {c(0,0)} else{c(0,1)}))
y_2n<-sapply(as.data.frame(y_n),function(x) rbind(x,x))

# B is the number of permutation test
# Q is the number of thread you used, within the range of the comuputer, the more threads, the faster you get the result
# The p-value of WCHC is given by
WCHC<-function(geno_2n,y_2n,B,Q){
  HCDC<-function(y_n){
    sumcol=ncol(y_n)
    colnames(y_n)=1:sumcol
    lis=list()
    h=rep(NA,sumcol-1)
    comb=rep(NA,sumcol-1)
    comb_NULL=rep(NA,sumcol-1)
    for (i in 1:(sumcol-1)) {
      if(i==1) {
        lis[i]=list(matrix(NA,choose(sumcol,2),3))
        lis[[i]][,1:2]=t(combn(sumcol,2))
        lis[[i]][,3]=sapply(1:choose(sumcol,2), function(x) {
          (1-abs(cor(y_n[,lis[[i]][x,1]],y_n[,lis[[i]][x,2]])))
        })
        lis_one=lis[[i]][which(lis[[i]][,3]==min(lis[[i]][,3])),]
        h[i]=lis_one[3]
        comb_NULL[i]=comb[i]=paste(lis_one[1],lis_one[2],sep='_')
      }
      else {
        lis[i]=list(matrix(NA,choose((sumcol-i+1),2),3))
        lis[[i]][,1:2]=t(combn(c(unique(comb_NULL[1:(i-1)]),(1:sumcol)[
          -as.numeric(unlist(strsplit(unique(comb_NULL[1:(i-1)]),split='_')))]),2))
        lis[[i]][,3]=sapply(1:choose(sumcol-i+1,2), function(x) {
          s1=unlist(strsplit(lis[[i]][x,1],split='_'))
          s2=unlist(strsplit(lis[[i]][x,2],split='_'))
          1-(cancor(y_n[,s1],y_n[,s2])$cor[1])
        })
        lis_two=lis[[i]][which(lis[[i]][,3]==min(lis[[i]][,3])),]
        h[i]=lis_two[3]
        comb_NULL[i]=comb[i]=paste(lis_two[1],lis_two[2],sep='_')
        index=sapply(1:(i-1),function(x) {
          (unlist(strsplit(comb[x],split = '_')) %in% unlist(strsplit(comb[i],split = '_')))[1]})
        comb_NULL[which(index)]=comb[i]
      }
    }
    h=as.numeric(h)
    D_value=do.call('c',sapply(1:(sumcol-1), function(x) {
      if (x>1) {h[x]-h[x-1]}
    }))
    lis_o=lis[[which(D_value==max(D_value))+1]]
    cluster_o=unique(rbind(matrix(lis_o[,1]),matrix(lis_o[,2])))
    clus=list()
    for (i in 1:nrow(cluster_o)) {
      if (length(grep("_",cluster_o[i,1])==1)) {clus[[i]]=as.numeric(unlist(strsplit(cluster_o[i,1],split = '_')))}
      else {clus[[i]]=as.numeric(cluster_o[i,1])}
    }
    return(list(n=nrow(cluster_o),cluster=clus))
  }
  HCDC_res=HCDC(y_2n)
  cluster=HCDC_res$cluster
  n=HCDC_res$n
  T=c()
  T=sapply(1:n,function(x) {TT(geno_2n,as.matrix(y_2n[,cluster[[x]]]))})
  func<-function(XXXX){
    TB=c()
    cluster=HCDC_res$cluster
    n=HCDC_res$n
    geno_2n_b=sample(geno_2n)
    TB=sapply(1:n,function(x) {TT(geno_2n_b,as.matrix(y_2n[,cluster[[x]]]))})
    return(TB)
  }
  cl <- makeCluster(Q)
  clusterExport(cl, c('TT','y_2n','geno_2n'))
  TB<-parLapply(cl,1:B, func)
  TB<-do.call("rbind",TB)
  stopCluster(cl)
  T<-rbind(T,TB)
  p_matrix=apply(as.matrix(T),2,function(x) {sapply(x,function(y) {(sum(y<=x)-1)/B})})
  p_vector=apply(p_matrix,1,min)
  p=(sum(p_vector<=p_vector[1])-1)/B
  return(p)
}
