library(rgl)

K=4#No. of classes 
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

W


open3d()

plot3d(W[,1], W[,2], W[,3], col = rainbow(1000))


a<-rbind(c(1,2,3),c(4,5,6))
a
WW<-rbind(c(1,1),c(2,2),c(3,3))
WW
g_x<-rbind(c(1,1),c(2,2))
g_x
a*exp(g_x%*%t(WW))

(order(-g_x[2,]%*%t(WW))[1])

max_anglevec.ind(g_x[2,],WW)

help(plot3d)
# for (seedd in 1:100){
#   set.seed(seedd)
#   if(sample(1:100, size = 1)==3){
#     print(seedd)
#     break
#   }
# }

# gamma.mat=matrix(1,2,3)
# flearn=matrix(1:10,2,2)
# 
# Cyik=matrix(2,2,3)
# 
# f_beta<-function(beta_b){
#   sum(Cyik*log(1+gamma.mat*exp(beta_b*(flearn%*%t(W)))))
# }
# optim(1, f_beta, method = "BFGS")
# gg=optim(1, f_beta, method = "BFGS")[[1]]
# gg
# 
# AAA=matrix(1:6,2,3)
# BBB=matrix(2:7,2,3)
# BBB/AAA
# BBB
# AAA
