###################### Asymptotic Bias and Variance of the Estimators #####################

options(digits=3); options(scipen=200)
x<-c(seq(.1,.6,.1),seq(.75,2,.25),3,5)

asyVar<-function(x){
  q<-exp(-x); p<-1-q
  q/(p-x*q)
}
asyBias.c<-function(x){
  q<-exp(-x); p<-1-q
  x*p*q/(2*(p-x*q)^2)-.5
}
asyBias.u<-function(x){
  q<-exp(-x); p<-1-q
  (x*q)^2/(2*(p-x*q)^2)-.5
}
asyBias.m<-function(x){
  q<-exp(-x); p<-1-q
  x*q*(x*q+2*(p-x*q)*(1/3-1/x))/(2*(p-x*q)^2)-.5
}

table2<-data.frame(lambda=x, sigma2=asyVar(x),
                   Nu=asyBias.u(x), Nc=asyBias.c(x), Nm=asyBias.m(x))



################# Preparing for Simulation #################

gpois<-function(n,lambda){
  # generating poisson rvs
  N<-numeric(n)
  for(i in 1:n){
    X<-0; P<-1
    while(P>=exp(-lambda)){
      X<-X+1
      U<-runif(1); P<-U*P
    }
    N[i]<-X-1
  }
  return(N)
}

constant<-function(Data){
  N=length(Data)
  n=sum(Data!=0)
  S=sum(Data)
  return(list(N=N,n=n,S=S))
}

gpois2<-function(n,lambda,seed=NULL){
  set.seed(seed)
  temp1<-temp2<-gpois(n,lambda)
  cons1<-constant(temp1)
  temp.S<-cons1$S; temp.n<-cons1$n
  
  while(temp.S<=(temp.n+1) | temp.n==0){
    temp2<-gpois(n,lambda)
    cons2<-constant(temp2)
    temp.S<-cons2$S; temp.n<-cons2$n
  }
  return(temp2)
}

N.est<-function(Data,method="m"){
  # computation of estimate of N
  # method is chosen among {u,c,m}
  
  const<-constant(Data)
  N<-const$N; n<-const$n; S<-const$S
  
  if(method=="c"){
    g<-function(x) return(x/(1-exp(-x))-S/n)
    lambda_c<-uniroot(g,lower=.1,upper=1000)$root
    p<-1-exp(-lambda_c); 
    Vc<-n/p; Nc<-floor(Vc)
    return(Nc)
  }
  
  if(method=="u"){
    D<-function(V) return(log(1-n/V)/log(1-1/V)-S)
    Vu<-uniroot(D,lower=n,upper=1000)$root
    Nu<-floor(Vu)
    return(Nu)
  }
  
  if(method=="m"){
    D<-function(V) return(log(1-n/V)/log(1-1/(V+1/3))-(S+1))
    Vm<-uniroot(D,lower=n,upper=1000)$root
    Nm<-floor(Vm)
    return(Nm)
  }
}

SingleMethod<-function(N,lambda,method="m",Nsim=1000,seeds=NULL){
  # N, lambda and method all take only one value

  estimate<-numeric(Nsim)
  
  for(i in 1:Nsim){
    Data<-gpois2(N,lambda,seed=seeds[i])
    estimate[i]<-N.est(Data,method=method)
  }
  Bias<-mean(estimate)-N
  MSE<-mean((estimate-N)^2)
  return(list(method=method,Bias=Bias,"MSE/N"=MSE/N))
}

CombinedMethod<-function(N,lambda,Nsim=1000,seeds=NULL){
  # N and lambda take only one value
  # all of the three methods are applied
  
  method<-c("u","c","m"); m<-length(method)
  Bias<-MSE.N<-numeric(m)
  for(i in 1:m){
    temp<-SingleMethod(N,lambda,method=method[i],Nsim,seeds)
    Bias[i]<-temp$Bias
    MSE.N[i]<-temp$"MSE/N"
  }
  result<-c(N,lambda,Bias,MSE.N)
  return(result)
}



######################  Estimated Biases and Mean Squared Errors #####################

simul1<-function(N,lambda,Nsim=1000,seeds=NULL){
  # N and lambda both are vectors
  
  N.len<-length(N); lambda.len<-length(lambda)
  
  result<-NULL
  for(i in 1:N.len){
    for(j in 1:lambda.len){
      result<-rbind(result,
                    CombinedMethod(N[i],lambda[j],Nsim,seeds))
    }
  }
  colnames(result)=c("N","lambda","Bias.u","Bias.c","Bias.m","MSE/N.u","MSE/N.c","MSE/N.m")
  return(as.data.frame(result))
}



###################### Estimated Bias and Mean Squared Error #####################
simul2<-function(N,lambda,Nsim=1000,seeds=NULL){
  # N and lambda both are vectors
  
  N.len<-length(N); lambda.len<-length(lambda)
  
  result<-NULL
  Bias<-MSE.N<-matrix(nr=lambda.len,nc=N.len)
  for(i in 1:lambda.len){
    for(j in 1:N.len){
      temp<-SingleMethod(N[j],lambda[i],seeds=seeds)
      Bias[i,j]<-temp$Bias
      MSE.N[i,j]<-temp$"MSE/N"
    }
  }
  
  result<-cbind(lambda,Bias,MSE.N)
  colnames(result)<-c("lambda",rep(as.character(N),2))
  return(as.data.frame(result))
}



######################## Simulation ########################

seeds<-13001:14000
table3.temp<-simul1(N=c(25,50,100),
                    lambda=c(seq(.75,2,.25),3,5),Nsim=1000,seeds=seeds)

table3<-cbind(table3.temp,'sigma2'=table2[7:14,2])

table4<-simul2(N=c(50,100),lambda=1:6/10,Nsim=1000,seeds=seeds)




