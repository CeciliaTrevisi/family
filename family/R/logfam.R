logfam<-function(Y,D,X){
  # Multi-dimensional gamma function (dim=q)
  loggammaq<-function(q,x){
    j<-(1:q)
    ((((q*(q-1))/4))*log(pi))+sum(lgamma(x+((1-j)/2)))
  }

  # Log-marginal likelihood on Y|X
  logmar<- function(Y,X){
    if(dim(Y)[2]==0){
      return(0)
    }
    else{
      X=cbind(X,1)
      q=as.numeric(dim(Y)[2])
      n=as.numeric(dim(Y)[1])
      p=as.numeric(dim(X)[2]-1)
      n0=p+2
      aD=q-1
      E=Y-X%*%(solve((t(X)%*%X)))%*%t(X)%*%Y
      z0=-(((n-n0)*q)/2)
      z1=(aD+n-p-1)/2
      z2=(aD+n0-p-1)/2
      z3=(q*(aD+n0))/2
      z4=-((n-n0)/2)
      return ((z0*log(pi))+loggammaq(q,z1)-loggammaq(q,z2)+(z3*log(n0/n))+(z4*log(abs(det(t(E)%*%E)))))
    }
  }
  q=as.numeric(dim(Y)[2])
  n=as.numeric(dim(Y)[1])
  fad=c()
  fam=c()
  for (i in 1:q){
    pad=c()
    for (j in 1:q){
      if (D[i,j]==1){
        pad=append(pad,j)}
    }

    figlio=i
    fad=append(pad,figlio)

    if(dim(matrix(Y[,pad],n))[2]==0){
      fam=append(fam,logmar(matrix(Y[,fad],n),X))
    }
    else{
      fam=append(fam,logmar(matrix(Y[,fad],n),X)-logmar(matrix(Y[,pad],n),X))
    }
  }
  return(sum(fam))
}
