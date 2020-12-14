# This function is adapted from Chang, Jinyuan, Eric D. Kolaczyk, and Qiwei Yao. "Estimation of subgraph densities in noisy networks." Journal of the American Statistical Association just-accepted (2020): 1-40.

getNetNoise<- function(alpha0, P, Y, Ystar, Ystar2, Nb=500, tiny=0.0001){
  
  # This program calculates point-estimates for two error rates Alpha, Beta.
  
  # P: number of nodes
  # Y, Ystar, Ystar2: 3 observed PxP adjacent matrices
  # Nb: No. bootstrap replications
  
  P012=P*(P-1)*(P-2); PP=P*P; P01=P*(P-1); P01h=P01/2; rP01h=sqrt(P01h)
  H=matrix(1:6, nrow=2)
  G=matrix(1:6, nrow=2)
  DD=matrix(1:4, nrow=2)
  S2=matrix(1:9, ncol=3)
  Sv=matrix(nrow=Nb, ncol=2)
  
  # Part 3: Point estimate for Alpha, Beta and Delta
  #     NO Part 1, Part 2!!!
  u1=sum(Y)/P01
  u2=sum(abs(Y-Ystar))/(2*P01)
  u3=length(Y[((Ystar2-2*Ystar+Y)==1)|((Ystar2-2*Ystar+Y)==-2)])/(3*P01)
  alpha=alpha0; alpha0=alpha+10*tiny
  while(abs(alpha-alpha0)>tiny) { alpha0=alpha
  beta=(u2-alpha0+alpha0*u1)/(u1-alpha0)
  delta=((u1-alpha0)^2)/(u1-u2-2*u1*alpha0+alpha0^2)
  alpha=(u3-delta*(1-beta)*beta^2)/((1-delta)*(1-alpha0)^2)
  }
  eS=rep(0,2)
  eS[1]=max(0,alpha); eS[2]=max(0,beta)
  return(eS)
} # The end
