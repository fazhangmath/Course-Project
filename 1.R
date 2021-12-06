library(MASS)
library(statmod)
penalty_r=matrix(nrow=5,ncol=50)
penalty_q=matrix(nrow=5,ncol=50)
p_list=c(0.5,0.73,0.82,0.88,0.95)
for(i1 in 1:5){
  p=p_list[i1]
  for(i2 in 1:50){
    sigma=i2/50*1.5
    R_q=1/2*sigma^2*p*(1-p)
    R=-log(1+p/(1-p))
    g=function(z) log(1+exp(log(p/(1-p))+sqrt(2*sigma^2)*z))
    hermite=gauss.quad(50,kind="hermite")
    R=R+sum(hermite$weights*g(hermite$nodes))/sqrt(pi)
    penalty_r[i1,i2]=R
    penalty_q[i1,i2]=R_q
  }
}
plot(c(1:50)/50*1.5,type="l",penalty_r[1,],col="black",lty=1,xlim=c(0,1.5),ylim=c(0,0.3),xlab="Sigma",ylab="Noising Penalty")
lines(c(1:50)/50*1.5,penalty_q[1,],col="black",lty=3)
lines(c(1:50)/50*1.5,penalty_r[2,],col="red",lty=1)
lines(c(1:50)/50*1.5,penalty_q[2,],col="red",lty=3)
lines(c(1:50)/50*1.5,penalty_r[3,],col="green",lty=1)
lines(c(1:50)/50*1.5,penalty_q[3,],col="green",lty=3)
lines(c(1:50)/50*1.5,penalty_r[4,],col="blue",lty=1)
lines(c(1:50)/50*1.5,penalty_q[4,],col="blue",lty=3)
lines(c(1:50)/50*1.5,penalty_r[5,],col="blueviolet",lty=1)
lines(c(1:50)/50*1.5,penalty_q[5,],col="blueviolet",lty=3)
legend("topleft",c("p=0.5","p=0.73","p=0.82","p=0.88","p=0.95"), col=c("black","red","green","blue","blueviolet"), lty=c(1,1,1,1,1))