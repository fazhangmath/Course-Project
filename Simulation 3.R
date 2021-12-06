count1=0
count2=0
count3=0
count4=0
for(ii in 1:100){
  library(MASS)
  X0=matrix(0,nrow=150,ncol=1050)
  y0=rep(0,150)
  delta=0.9
  lambda=32
  beta0=c(rep(0.057,50),rep(0,1000))
  for(index in 1:6){
    for(g in 1:25){
      i=g+(index-1)*25
      sgn=(rbinom(1,1,0.5)-0.5)*2
      if(g<=5){
        for(j in (10*(g-1)+1):(10*(g-1)+10)){
          X0[i,j]=sgn*abs(4.3976*rnorm(1,mean=0,sd=1))
        }
      }
      X0[i,51:1050]=rnorm(1000,mean=0,sd=1)
      y0[i]=rbinom(1,1,1/(1+exp(-X0[i,]%*%beta0)))
    }
  }
  X=X0[1:75,]
  y=y0[1:75]
  for(i in 1:1050){
    X[,i]=X[,i]/sqrt(sum(X[,i]^2))
  }
  X_test=X0[76:150,]
  y_test=y0[76:150]
  fr=function(beta){
    p=1/(1+exp(-X%*%beta))
    temp=0
    for(i in 1:1050){
      temp=temp+sum(p*(1-p)*X[,i]^2*beta[i]^2)
    }
    -sum(y*(X%*%beta))+sum(log(1+exp(X%*%beta)))+temp*0.5*delta/(1-delta)
  }
  grr=function(beta){
    p=1/(1+exp(-X%*%beta))
    temp0=-exp(X%*%beta)*(exp(2*X%*%beta)-1)/(1+exp(X%*%beta))^4
    temp1=exp(X%*%beta)/(1+exp(X%*%beta))
    temp=rep(0,1050)
    for(i in 1:1050){
      temp[i]=-sum(y*X[,i])+sum(temp1*X[,i])
      temp[i]=temp[i]+2*sum(p*(1-p)*X[,i]^2*beta[i])*0.5*delta/(1-delta)+sum(temp0*X[,i]^3*beta[i]^2)*0.5*delta/(1-delta)
    }
    temp
  }
  beta=rep(0,1050)
  res=optim(beta,fr,grr,method="BFGS")
  beta=res$par
  beta_ridge=solve(t(X)%*%X+diag(rep(lambda,1050)))%*%t(X)%*%y
  p_test1=1/(1+exp(-X_test%*%beta))
  p_test2=1/(1+exp(-X_test%*%beta_ridge))
  for(i in 1:75){
    if(p_test1[i]>0.5&y_test[i]==1){
      count1=count1+1
    }
    if(p_test1[i]<=0.5&y_test[i]==0){
      count1=count1+1
    }
    if(p_test2[i]>0.5&y_test[i]==1){
      count2=count2+1
    }
    if(p_test2[i]<=0.5&y_test[i]==0){
      count2=count2+1
    }
  }
  X_active=matrix(nrow=15,ncol=1050)
  y_active=rep(0,15)
  for(index in 4:6){
    for(g in 1:5){
      i=g+(index-1)*25
      i0=(index-4)*5+g
      X_active[i0,]=X0[i,]
      y_active[i0]=y0[i]
    }
  }
  p_test1=1/(1+exp(-X_active%*%beta))
  p_test2=1/(1+exp(-X_active%*%beta_ridge))
  for(i in 1:15){
    if(p_test1[i]>0.5&y_active[i]==1){
      count3=count3+1
    }
    if(p_test1[i]<=0.5&y_active[i]==0){
      count3=count3+1
    }
    if(p_test2[i]>0.5&y_active[i]==1){
      count4=count4+1
    }
    if(p_test2[i]<=0.5&y_active[i]==0){
      count4=count4+1
    }
  }
}
cat(" The accuracy of L2-regularization for active instances is: ", count4/15/100,"\n")
cat(" The accuracy of dropout training for active instances is: ", count3/15/100,"\n")
cat(" The accuracy of L2-regularization for all instances is: ", count2/75/100,"\n")
cat(" The accuracy of dropout training for all instances is: ", count1/75/100,"\n")