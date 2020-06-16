K=4# pls input the No. of classes 

one_vec=rep(1,K-1)

e_j=function(K,j){
  if(j>K-1){
    print("Error: j must <=K-1")
  }
  else{
    return(c(rep(0,j-1),1,rep(0,K-1-j)))
  }  
}


W=matrix(1,K,K-1)

for (j in 1:K){
  if(j==1){
    W[j,]=one_vec*(K-1)^(-0.5)
  }
  else{
    W[j,]=-(1+sqrt(K))/{(K-1)^1.5}+e_j(K,j-1)*{K/(K-1)}^0.5
  }
}

print(W)
