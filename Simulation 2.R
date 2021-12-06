library(R.matlab)
data=readMat("./example_data.mat")
X0=data$X
y0=data$y
for(i in 1:1427){
  if(y0[i]==-1){
    y0[i]=0
  }
}
X0=matrix(X0,nrow=1427,ncol=22178,byrow=TRUE)
X=matrix(nrow=1000,ncol=22178)
X_test=matrix(nrow=427,ncol=22178)
y=rep(0,1000)
y_test=rep(427)
X[1:500,]=X0[1:500,]
X[501:1000,]=X0[800:1299,]
X_test[1:299,]=X0[501:799,]
X_test[300:427,]=X0[1300:1427,]
y[1:500]=y0[1:500]
y[501:1000]=y0[800:1299]
y_test[1:299]=y0[501:799]
y_test[300:427]=y0[1300:1427]
for(i in 1:22178){
  if(sqrt(sum(X[,i]^2))==0){
    next
  }
  X[,i]=X[,i]/sqrt(sum(X[,i]^2))
}
for(i in 1:22178){
  if(sqrt(sum(X_test[,i]^2))==0){
    next
  }
  X_test[,i]=X_test[,i]/sqrt(sum(X_test[,i]^2))
}
fr=function(beta){
  p=1/(1+exp(-X%*%beta))
  temp=0
  for(i in 1:22178){
    temp=temp+sum(p*(1-p)*X[,i]^2*beta[i]^2)
  }
  -sum(y*(X%*%beta))+sum(log(1+exp(X%*%beta)))+temp
}
grr=function(beta){
  p=1/(1+exp(-X%*%beta))
  temp0=-exp(X%*%beta)*(exp(2*X%*%beta)-1)/(1+exp(X%*%beta))^4
  temp1=exp(X%*%beta)/(1+exp(X%*%beta))
  temp=rep(0,22178)
  for(i in 1:22178){
    temp[i]=-sum(y*X[,i])+sum(temp1*X[,i])
    temp[i]=temp[i]+2*sum(p*(1-p)*X[,i]^2*beta[i])+sum(temp0*X[,i]^3*beta[i]^2)
  }
  temp
}
beta_record=matrix(nrow=30,ncol=22178)
beta=rep(0,22178)
for(i in 1:30){
  beta_record[i,]=beta
  res=optim(beta,fr,grr,method="BFGS",control=list(trace=6,maxit=1))
  beta=res$par
}
loss=rep(0,30)
penalty_d=rep(0,30)
penalty_q=rep(0,30)
for(i in 1:30){
  beta_temp=beta_record[i,]
  loss[i]=fr(beta_temp)
  p=1/(1+exp(-X%*%beta_temp))
  temp=0
  for(j in 1:22178){
    temp=temp+sum(p*(1-p)*X[,j]^2*beta_temp[j]^2)
  }
  penalty_q[i]=temp
}
index=list(which(X[1,]!=0))
for(i in 2:1000){
  index[[i]]=which(X[i,]!=0)
}
for(i in 1:30){
  temp=0
  beta_temp=beta_record[i,]
  for(k in 1:1000){
    temp1=X[k,][index[[k]]]
    temp2=beta_temp[index[[k]]]
    noise=rbinom(length(index[[k]])*500,1,0.3333)
    noise=matrix(noise,nrow=500,ncol=length(index[[k]]),byrow=TRUE)
    noise=noise/0.3333
    temp1=sweep(noise,MARGIN=2,temp1,"*")
    temp=temp+sum(log(1+exp(temp1%*%temp2)))
  }
  temp=temp-500*sum(log(1+exp(X%*%beta_temp)))
  penalty_d[i]=temp/500
}
plot(log(loss+1),type="l",ylim=c(0,7),yaxt="n",ylab="Loss",xlab="Training Iteration")
axis(side=2,at=log(c(10,20,50,100,200,500)),labels=c(10,20,50,100,200,500))
lines(log(penalty_q+1),col="red",lty=3)
lines(log(penalty_d+1),col="blue",lty=3)
legend("topright",c("Dropout Penalty","Quadratic Penalty","Negative Log-Likelihood"), col=c("blue","red","black"), lty=c(3,3,1))