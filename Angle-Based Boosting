trsf=function(kvec, K) {
    res=matrix(-1/(K-1),length(kvec),K)
    for ( i in 1:length(kvec) ) res[i,kvec[i]]=1
    return(res)
}

trsf1=function(kvec, K) {
    res=matrix(-1,length(kvec),K)
    for ( i in 1:length(kvec) ) res[i,kvec[i]]=1
    return(res)
}

trsf2=function(kvec, K) {
    res=matrix(0,length(kvec),K)
    for ( i in 1:length(kvec) ) res[i,kvec[i]]=1
    return(res)
}

std<-function(x) { return((x-mean(x))/sd(x)) }
max.ind=function(x) return(order(-x)[1])

my.optimize=function(fn,lower,upper,n.try=100) {
   xseq=seq(lower,upper,length=n.try)
   res=rep(0,n.try)
   for ( i in 1:n.try) res[i]=fn(xseq[i])
   min.x=xseq[order(res)[1]]
   obj=min(res)
   return(list(minimum=min.x,objective=obj))
}

split.m=function(xlearn, ylearn, xtest, mat) {
    n=nrow(xlearn)
    n.test=nrow(xtest)
    p=ncol(xlearn)
    K=ncol(mat)
    
    res=matrix(0,p,2)
    for ( j in 1:p ) {
        f.foo=function(aa) {
            aa.ind=which(xlearn[,j]>=aa)
            mat.ind=mat[aa.ind,]
            if (!is.matrix(mat.ind)) mat.ind=t(mat.ind)
            mat.ind.n=mat[-aa.ind,]
            if (!is.matrix(mat.ind.n)) mat.ind.n=t(mat.ind.n)
            return( min(apply(mat.ind,2,sum))+min(apply(mat.ind.n,2,sum)) ) 
        }
        if (sd(xlearn[,j])>10^(-3)) {
           foo=my.optimize(f.foo,min(xlearn[,j]),max(xlearn[,j]))
           res[j,1]=foo$minimum
           res[j,2]=foo$objective
        } else {
           res[j,1]=0
           res[j,2]=10^6
        }
    }
    j.best=order(res[,2])[1]
    learn.ind=as.numeric(xlearn[,j.best]>=res[j.best,1])+1
    test.ind=as.numeric(xtest[,j.best]>=res[j.best,1])+1
    
    return(list(learn.ind=learn.ind,test.ind=test.ind))
}

learner.m=function(xlearn, ylearn, xtest, wgt=NULL, mat=NULL, level=2) {
    n=nrow(xlearn)
    n.test=nrow(xtest)
    p=ncol(xlearn)
    K=length(table(ylearn))
    if (is.vector(wgt)) {
       mat=matrix(rep(wgt/(K-1),K),ncol=K) 
       mat[cbind(1:n,ylearn)]=0    
    }

    learn.subset=rep(1,n)   
    test.subset=rep(1,n.test) 
    for ( l in 1:level ) {
        learn.foo=learn.subset
        test.foo=test.subset
        for ( ll in 1:(2^(l-1)) ) {
            if (sum(learn.subset==ll)>1) {
                split.foo=split.m(xlearn[learn.subset==ll,],ylearn[learn.subset==ll],xtest[test.subset==ll,],mat[learn.subset==ll,])
                learn.tmp=split.foo$learn.ind
                test.tmp=split.foo$test.ind
            } else {
                learn.tmp=1
                test.tmp=rep(1,sum(test.subset==ll))
            }
            learn.foo[learn.subset==ll]=2*(ll-1)+learn.tmp
            test.foo[test.subset==ll]=2*(ll-1)+test.tmp
        }
        learn.subset=learn.foo
        test.subset=test.foo
    }
    
    ylearn.pred=rep(0,n)
    ytest.pred=rep(0,n.test)
    for ( lll in 1:(2^level) ) {
        if (sum(learn.subset==lll)==0) next
        pred.tmp=ifelse(sum(learn.subset==lll)==1,ylearn[learn.subset==lll],order(apply(mat[learn.subset==lll,],2,sum))[1])
        ylearn.pred[learn.subset==lll]=pred.tmp
        ytest.pred[test.subset==lll]=pred.tmp
    }
    
    return(list(learn=ylearn.pred,test=ytest.pred))
}


adaboost.my <- function(x.train, y.train, x.test, y.test, cost=NULL, mfinal = 200) {
    ## Initialization
    n=nrow(x.train)
    p=ncol(x.train)
    K=length(unique(y.train))
    n.test=nrow(x.test)         
    if (is.matrix(cost)) c.mat=cost else c.mat=1-diag(K)
    
    Flearn=matrix(0,n,K)
    Ftest=matrix(0,n.test,K)

    a.mat=matrix(0,n,K)
    for ( i in 1:n ) a.mat[i,]=c.mat[y.train[i],]
    
    ## Boosting Iterations
    
    err=matrix(0,mfinal,2)
    for (m in 1:mfinal) {
        ## Fitting the tree
        update=learner.m(x.train, y.train, x.test, mat=a.mat, level=1)
        flearn=trsf(update$learn,K)
        ftest=trsf(update$test,K)
        f.err=sum(a.mat[cbind(1:n,update$learn)])
         
        ## Updating
        if (f.err>0) {
            ## f.foo=function(aa) return(sum(a.mat*exp(aa*flearn))) 
            ## f.coef=my.optimize(f.foo,0,100)$minimum   ### Line search        
            f.coef=(K-1)/K*(log((sum(a.mat)-f.err)/f.err)-log(K-1))
            ## f.coef=1/n
        
            Flearn=Flearn + f.coef*flearn
            Ftest=Ftest + f.coef*ftest
            a.mat=a.mat*exp(f.coef*flearn)
            if (sum(a.mat)==0) break 
            a.mat=pmin(a.mat/sum(a.mat),1e24)
        }
        if (f.err==0) {
            a.mat=a.mat
            Flearn=Flearn + flearn
            Ftest=Ftest + ftest
        }
        
        err[m,1]=sum(c.mat[cbind(y.train,apply(Flearn,1,max.ind))])/n
        err[m,2]=sum(c.mat[cbind(y.test,apply(Ftest,1,max.ind))])/n.test
    }
    
    return(list(ytr=apply(Flearn,1,max.ind),yte=apply(Ftest,1,max.ind),error=err[1:m,]))   
}



#####################################################################################
################################## Existing methods #################################
#####################################################################################

learner <- function(xlearn, ylearn, xtest, w) {
    ## Definitions
    require(rpart)
    
    learn  <- dim(xlearn)[1]
    test   <- dim(xtest)[1]
    K = length(unique(ylearn))
    
    ## Currently only stumps as learners are supported, no choice of args!!!
    cntrl <- rpart.control(maxdepth=1)
    
    ## Bagging stumps/trees
   	tmp=data.frame(y=ylearn,xlearn)
    fit=rpart(factor(y)~.,data=tmp, weights = w/mean(w), control = cntrl)
	blearn=predict(fit, type="class")    
	tmp=data.frame(y=rep(0,nrow(xtest)),xtest)
	btest=predict(fit, tmp, type="class")

    ## Output
    list(learn = as.numeric(as.character(blearn)), test = as.numeric(as.character(btest)))
}

adaboost <- function(x.train, y.train, x.test, y.test, mfinal = 200, wgt=0.5) {
    ## Initialization
    n=nrow(x.train)
    if (min(y.train)<0) y.train=(y.train+1)/2
    n.test=nrow(x.test)         
    if (length(wgt)==1) { w=rep(wgt,n); w[y.train==1]=1-wgt } else w=wgt
    Flearn=rep(0,n)
    Ftest=rep(0,n.test)

    err=matrix(0,mfinal,2)
    Flearn.mat=matrix(0,n,mfinal)
    Ftest.mat=matrix(0,n.test,mfinal)
    ## Boosting Iterations
    for (m in 1:mfinal) {
        ## Fitting the tree
        update=learner(x.train, y.train, x.test, w)
        flearn=update$learn
        ftest=update$test
         
        ## Updating
        f.err=sum(w*(flearn!=y.train))/sum(w)
        if (f.err>0) {
            f.coef=log((1-f.err)/f.err)
            w=w*exp(f.coef*(flearn!=y.train))     
            w=pmax(w/sum(w), 1e-24)
            Flearn=Flearn + f.coef*flearn
            Ftest=Ftest + f.coef*ftest
        }
        if (f.err==0) {
            w=w
            Flearn=Flearn + flearn
            Ftest=Ftest + ftest
        }
        err[m,1]=sum(sign(Flearn)!=y.train)/n
        err[m,2]=sum(sign(Ftest)!=y.test)/n.test
        Flearn.mat[,m]=Flearn
        Ftest.mat[,m]=Ftest
    }
       
    return(list(Flearn=Flearn,Ftest=Ftest,Flearn.mat=Flearn.mat,Ftest.mat=Ftest.mat,error=err))   
}


adaboost.mh <- function(x.train, y.train, x.test, y.test, mfinal = 200) {
    ## Initialization
    n=nrow(x.train)
    p=ncol(x.train)
    K=length(unique(y.train))
    n.test=nrow(x.test)         
    
    Flearn=matrix(0,n,K)
    Ftest=matrix(0,n.test,K)
    w=rep(1,n)

    err=matrix(0,mfinal,2)
    Flearn.array=array(0,c(K,n,mfinal))
    Ftest.array=array(0,c(K,n.test,mfinal))
    for ( j in 1:K ) {
        y.train.j=2*as.numeric(y.train==j)-1
        y.test.j=2*as.numeric(y.test==j)-1
    
        foo=adaboost(x.train, y.train.j, x.test, y.test.j, mfinal)
        Flearn[,j]=foo$Flearn
        Ftest[,j]=foo$Ftest
    
        Flearn.array[j,,]=foo$Flearn.mat
        Ftest.array[j,,]=foo$Ftest.mat
    }
    
    for ( m in 1:mfinal) {
        Flearn.tmp=t(Flearn.array[,,m])
        Ftest.tmp=t(Ftest.array[,,m])
        err[m,1]=sum(apply(Flearn.tmp,1,max.ind)!=y.train)/n #error rate on the training set
        err[m,2]=sum(apply(Ftest.tmp,1,max.ind)!=y.test)/n.test #error rate on the test set
    }
    
    return(list(ytr=apply(Flearn,1,max.ind),yte=apply(Ftest,1,max.ind),error=err))   
}


adaboost.mg <- function(x.train, y.train, x.test, y.test, cost=NULL, mfinal = 200) {
    ## Initialization
    n=nrow(x.train)
    p=ncol(x.train)
    K=length(unique(y.train))
    n.test=nrow(x.test)         
    if (is.matrix(cost)) c.mat=cost else c.mat=matrix(1,K,K)
    
    Flearn=matrix(0,n,K)
    Ftest=matrix(0,n.test,K)
    w=rep(1,n)

    err=matrix(0,mfinal,2)
    ## Boosting Iterations
    for (m in 1:mfinal) {
    
        ## Fitting the tree
###        wmat=w*matrix(1,n,K); wmat[cbind(1:n,y.train)]=0
###        update=learner.m(x.train, y.train, x.test, mat=wmat, level=1)
        update=learner.m(x.train, y.train, x.test, wgt=w, level=1)
        flearn=matrix(0,n,K)
        flearn[cbind(1:n,update$learn)]=1
        ftest=matrix(0,n.test,K)
        ftest[cbind(1:n.test,update$test)]=1
         
        ## Updating
        f.err=sum(w*(update$learn!=y.train))/sum(w)
        if (f.err>0) {
            f.coef=(log((1-f.err)/f.err)+log(K-1))
            Flearn=Flearn + f.coef*flearn
            Ftest=Ftest + f.coef*ftest
            w=w*exp(f.coef*(update$learn!=y.train))
            w=pmax(w/sum(w), 1e-24)
        }
        if (f.err==0) {
            Flearn=Flearn + flearn
            Ftest=Ftest + ftest
        }
        
        err[m,1]=sum(apply(Flearn,1,max.ind)!=y.train)/n
        err[m,2]=sum(apply(Ftest,1,max.ind)!=y.test)/n.test
    
    }
    return(list(ytr=apply(Flearn,1,max.ind),yte=apply(Ftest,1,max.ind),error=err))   
}


adaboost.m2 <- function(x.train, y.train, x.test, y.test, cost=NULL, mfinal = 200) {
    ## Initialization
    n=nrow(x.train)
    p=ncol(x.train)
    K=length(unique(y.train))
    n.test=nrow(x.test)         
    if (is.matrix(cost)) c.mat=cost else c.mat=1-diag(K)
    
    Flearn=matrix(0,n,K)
    Ftest=matrix(0,n.test,K)

    a.mat=matrix(0,n,K)
    for ( i in 1:n ) a.mat[i,]=c.mat[y.train[i],]
    
    ## Boosting Iterations
    
    err=matrix(0,mfinal,2)
    for (m in 1:mfinal) {
        ## Fitting the tree
        update=learner.m(x.train, y.train, x.test, mat=a.mat, level=1)
        flearn=trsf1(update$learn,K)
        ftest=trsf1(update$test,K)
        f.err=sum(a.mat[cbind(1:n,update$learn)])
        flearn.y=flearn
        for ( i in 1:n ) flearn.y[i,]=rep(flearn[i,y.train[i]],K)
         
        ## Updating
        if (f.err>0) {
            ## f.foo=function(aa) return(sum(a.mat*exp(aa*flearn))) 
            ## f.coef=my.optimize(f.foo,0,100)$minimum   ### Line search        
            w.pos=w.neg=0
            for ( i in 1:n ) {
               w.pos=w.pos+sum(a.mat[i,flearn[i,]>flearn.y[i,]])
               w.neg=w.neg+sum(a.mat[i,flearn[i,]<flearn.y[i,]])
            }
            f.coef=0.5*log(w.pos/w.neg)
            f.coef=1/n
        
            Flearn=Flearn + f.coef*flearn
            Ftest=Ftest + f.coef*ftest
            a.mat=a.mat*exp(0.5*f.coef*(flearn-flearn.y))
            if (sum(a.mat)==0) break 
            a.mat=pmin(a.mat/sum(a.mat),1e24)
        }
        if (f.err==0) {
            a.mat=a.mat
            Flearn=Flearn + flearn
            Ftest=Ftest + ftest
        }
        
        err[m,1]=sum(c.mat[cbind(y.train,apply(Flearn,1,max.ind))])/n
        err[m,2]=sum(c.mat[cbind(y.test,apply(Ftest,1,max.ind))])/n.test
    }
    
    return(list(ytr=apply(Flearn,1,max.ind),yte=apply(Ftest,1,max.ind),error=err[1:m,]))   
}


adaboost.lp <- function(x.train, y.train, x.test, y.test, pnorm=2, cost=NULL, mfinal = 200) {
    ## Initialization
    n=nrow(x.train)
    p=ncol(x.train)
    K=length(unique(y.train))
    n.test=nrow(x.test)         
    if (is.matrix(cost)) c.mat=cost else c.mat=1-diag(K)
    
    Flearn=matrix(1/K,n,K)
    Ftest=matrix(1/K,n.test,K)

    a.mat=matrix(0,n,K)
    for ( i in 1:n ) a.mat[i,]=c.mat[y.train[i],]
    
    ## Boosting Iterations
    
    err=matrix(0,mfinal,2)
    for (m in 1:mfinal) {
        ## Fitting the tree
        update=learner.m(x.train, y.train, x.test, mat=a.mat, level=1)
        flearn=trsf2(update$learn,K)
        ftest=trsf2(update$test,K)
        f.err=sum(a.mat[cbind(1:n,update$learn)])
         
        ## Updating
        if (f.err>0) {
            ## f.foo=function(aa) return(sum(a.mat*exp(aa*flearn))) 
            ## f.coef=my.optimize(f.foo,0,100)$minimum   ### Line search        
            f.coef=1/n
        
            Flearn=(1-f.coef)*Flearn + f.coef*flearn
            Ftest=(1-f.coef)*Ftest + f.coef*ftest
            for ( i in 1:n ) {
               a.mat[i,]=c.mat[y.train[i],]*(Flearn[i,]^pnorm)
            }
            if (sum(a.mat)==0) break 
            a.mat=pmin(a.mat/sum(a.mat),1e24)
        }
        if (f.err==0) {
            a.mat=a.mat
            Flearn=(1-1/n)*Flearn + flearn/n
            Ftest=(1-1/n)*Ftest + ftest/n
        }
        
        err[m,1]=sum(c.mat[cbind(y.train,apply(Flearn,1,max.ind))])/n
        err[m,2]=sum(c.mat[cbind(y.test,apply(Ftest,1,max.ind))])/n.test
    }
    
    return(list(ytr=apply(Flearn,1,max.ind),yte=apply(Ftest,1,max.ind),error=err[1:m,]))   
}





#The following code are developped by YANG Yi

##############################################################################################
###################Cost-sensitive error extension of existing methods#########################
##############################################################################################


adaboost.mhCS <- function(x.train, y.train, x.test, y.test, cost=NULL, mfinal = 200) {
  ## Initialization
  n=nrow(x.train)
  p=ncol(x.train)
  K=length(unique(y.train))
  n.test=nrow(x.test)
  if (is.matrix(cost)) c.mat=cost else c.mat=1-diag(K)
  
  Flearn=matrix(0,n,K)
  Ftest=matrix(0,n.test,K)
  w=rep(1,n)
  
  err=matrix(0,mfinal,2)
  err_CS=matrix(0,mfinal,2)
  
  Flearn.array=array(0,c(K,n,mfinal))
  Ftest.array=array(0,c(K,n.test,mfinal))
  for ( j in 1:K ) {
    y.train.j=2*as.numeric(y.train==j)-1
    y.test.j=2*as.numeric(y.test==j)-1
    
    foo=adaboost(x.train, y.train.j, x.test, y.test.j, mfinal)
    Flearn[,j]=foo$Flearn
    Ftest[,j]=foo$Ftest
    
    Flearn.array[j,,]=foo$Flearn.mat
    Ftest.array[j,,]=foo$Ftest.mat
  }
  
  for ( m in 1:mfinal) {
    Flearn.tmp=t(Flearn.array[,,m])
    Ftest.tmp=t(Ftest.array[,,m])
    err[m,1]=sum(apply(Flearn.tmp,1,max.ind)!=y.train)/n #error rate on the training set
    err[m,2]=sum(apply(Ftest.tmp,1,max.ind)!=y.test)/n.test #error rate on the test set
    
    err_CS[m,1]=sum(c.mat[cbind(y.train,apply(Flearn.tmp,1,max.ind))])/n
    err_CS[m,2]=sum(c.mat[cbind(y.test,apply(Ftest.tmp,1,max.ind))])/n.test
  }
  
  return(list(ytr=apply(Flearn,1,max.ind),yte=apply(Ftest,1,max.ind),error=err_CS))   
}


adaboost.mgCS <- function(x.train, y.train, x.test, y.test, cost=NULL, mfinal = 200) {
  ## Initialization
  n=nrow(x.train)
  p=ncol(x.train)
  K=length(unique(y.train))
  n.test=nrow(x.test)         
  if (is.matrix(cost)) c.mat=cost else c.mat=matrix(1,K,K)
  
  Flearn=matrix(0,n,K)
  Ftest=matrix(0,n.test,K)
  w=rep(1,n)
  
  err=matrix(0,mfinal,2)
  ## Boosting Iterations
  for (m in 1:mfinal) {
    
    ## Fitting the tree
    ###        wmat=w*matrix(1,n,K); wmat[cbind(1:n,y.train)]=0
    ###        update=learner.m(x.train, y.train, x.test, mat=wmat, level=1)
    update=learner.m(x.train, y.train, x.test, wgt=w, level=1)
    flearn=matrix(0,n,K)
    flearn[cbind(1:n,update$learn)]=1
    ftest=matrix(0,n.test,K)
    ftest[cbind(1:n.test,update$test)]=1
    
    ## Updating
    f.err=sum(w*(update$learn!=y.train))/sum(w)
    if (f.err>0) {
      f.coef=(log((1-f.err)/f.err)+log(K-1))
      Flearn=Flearn + f.coef*flearn
      Ftest=Ftest + f.coef*ftest
      w=w*exp(f.coef*(update$learn!=y.train))
      w=pmax(w/sum(w), 1e-24)
    }
    if (f.err==0) {
      Flearn=Flearn + flearn
      Ftest=Ftest + ftest
    }
    
    # err[m,1]=sum(apply(Flearn,1,max.ind)!=y.train)/n
    # err[m,2]=sum(apply(Ftest,1,max.ind)!=y.test)/n.test
    err[m,1]=sum(c.mat[cbind(y.train,apply(Flearn,1,max.ind))])/n
    err[m,2]=sum(c.mat[cbind(y.test,apply(Ftest,1,max.ind))])/n.test
    
  }
  return(list(ytr=apply(Flearn,1,max.ind),yte=apply(Ftest,1,max.ind),error=err))   
}



###################################################################################################
################################## My Angle-based Boosting ########################################
###################################################################################################


## Convert phi(x)(int class label) to g(x)(angle-based vertex)
trsf_g=function(kvec, K) {
  #K:the No. of the classes
  #kvec: the class vector of the input instances
  
  res=matrix(0,length(kvec),K-1)
  for ( i in 1:length(kvec) ) res[i,]=W[kvec[i],]
  return(res)
}

#trsf_g(c(2,2,2),3)

## Convert g(x)(angle-based vertex) to phi(x)(int class label)
max_anglevec.ind=function(x) return(order(-x%*%t(W))[1])




##################### My angle-based cost-sensitive adaboost classifier############################

adaboost.my_CSangle <- function(x.train, y.train, x.test, y.test, cost=NULL, mfinal = 200) {
  ## Initialization
  n=nrow(x.train)
  p=ncol(x.train)
  K=length(unique(y.train))
  n.test=nrow(x.test) 
  
  if (is.matrix(cost)) c.mat=cost else c.mat=1-diag(K)
  
  Flearn=matrix(0,n,K-1)
  Ftest=matrix(0,n.test,K-1)
  
  a.mat=matrix(0,n,K)
  for ( i in 1:n ) a.mat[i,]=c.mat[y.train[i],]
  
  ## Boosting Iterations
  
  err=matrix(0,mfinal,2)
  for (m in 1:mfinal) {
    ## Fitting the tree
    update=learner.m(x.train, y.train, x.test, mat=a.mat, level=1)
    flearn=trsf_g(update$learn,K)
    ftest=trsf_g(update$test,K)
    f.err=sum(a.mat[cbind(1:n,update$learn)])
    
    ## Updating
    if (f.err>0) {
      ## f.foo=function(aa) return(sum(a.mat*exp(aa*flearn))) 
      ## f.coef=my.optimize(f.foo,0,100)$minimum   ### Line search        
      f.coef=(K-1)/K*(log((sum(a.mat)-f.err)/f.err)-log(K-1))
      ## f.coef=1/n
      
      Flearn=Flearn + f.coef*flearn
      Ftest=Ftest + f.coef*ftest
      
      a.mat=a.mat*exp(f.coef*(flearn%*%t(W)))
      
      if (sum(a.mat)==0) break 
      a.mat=pmin(a.mat/sum(a.mat),1e24)
    }
    if (f.err==0) {
      a.mat=a.mat
      Flearn=Flearn + flearn
      Ftest=Ftest + ftest
    }
    
    err[m,1]=sum(c.mat[cbind(y.train,apply(Flearn,1,max_anglevec.ind))])/n
    err[m,2]=sum(c.mat[cbind(y.test,apply(Ftest,1,max_anglevec.ind))])/n.test
  }
  
  return(list(ytr=apply(Flearn,1,max_anglevec.ind),yte=apply(Ftest,1,max_anglevec.ind),error=err[1:m,]))   
}


##################### My angle-based cost-sensitive logitboost classifier############################

logitboost.my_CSangle <- function(x.train, y.train, x.test, y.test, cost=NULL, mfinal = 200) {
  ## Initialization
  n=nrow(x.train)
  p=ncol(x.train)
  K=length(unique(y.train))
  n.test=nrow(x.test) 
  
  if (is.matrix(cost)) c.mat=cost else c.mat=1-diag(K)
  
  Flearn=matrix(0,n,K-1)
  Ftest=matrix(0,n.test,K-1)
  
  a.mat=matrix(0,n,K)
  Cyik=matrix(0,n,K)
  gamma.mat=matrix(1,n,K)# new added
  
  for ( i in 1:n ) {
    a.mat[i,]=c.mat[y.train[i],]
    Cyik[i,]=c.mat[y.train[i],]
    }
  
  
  ## Boosting Iterations
  
  err=matrix(0,mfinal,2)
  for (m in 1:mfinal) {
    ## Fitting the tree
    update=learner.m(x.train, y.train, x.test, mat=a.mat, level=1)#caculate the optimal phi(m+1)
    flearn=trsf_g(update$learn,K)# convert phi(m+1) to g(m+1) for learning set
    ftest=trsf_g(update$test,K)# convert phi(m+1) to g(m+1) for testing set
    f.err=sum(a.mat[cbind(1:n,update$learn)])
    
    ## Updating
    if (TRUE) {
      ## f.foo=function(aa) return(sum(a.mat*exp(aa*flearn))) 
      ## f.coef=my.optimize(f.foo,0,100)$minimum   ### Line search        
      #f.coef=(K-1)/K*(log((sum(a.mat)-f.err)/f.err)-log(K-1))
      
      f_beta<-function(beta_b){
        sum(Cyik*log(1+gamma.mat*exp(beta_b*(flearn%*%t(W)))))
      }
      
      f.coef=optim(1, f_beta, method = "BFGS")[[1]]
      
      ## f.coef=1/n
      
      Flearn=Flearn + f.coef*flearn
      Ftest=Ftest + f.coef*ftest
      
      gamma.mat=gamma.mat*exp(f.coef*(flearn%*%t(W)))
      a.mat=(Cyik*gamma.mat)/(1+gamma.mat)

      
      if (sum(a.mat)==0) break 
      
      a.mat=pmin(a.mat/sum(a.mat),1e24)
    }
    

    err[m,1]=sum(c.mat[cbind(y.train,apply(Flearn,1,max_anglevec.ind))])/n
    err[m,2]=sum(c.mat[cbind(y.test,apply(Ftest,1,max_anglevec.ind))])/n.test
  }
  
  return(list(ytr=apply(Flearn,1,max_anglevec.ind),yte=apply(Ftest,1,max_anglevec.ind),error=err[1:m,]))   
}
