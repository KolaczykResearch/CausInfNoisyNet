# Part of the getEst function is adapted from Kolaczyk, Eric D. Topics at the Frontier of Statistics and Network Analysis:(re) visiting the Foundations. Cambridge University Press, 2017.

# This program compute 3 types of causal effect estimators (A&S estimators in the true network, A&S estimators in the noisy network, and MME in the noisy network), A&S variance estimators in the true network and MME variance estimators in the noisy network.

# data.true: adjacency matrix of the true network
# data.obs: adjacency matrix of the noisy network
# p: the probability for Bernoulli random assignment of treatment
# Nv: number of nodes
# alphaEst: point-estimate for alpha (Type I error)
# betaEst: point-estimate for beta (Type II error)
# z: treatment assignment vector
# O.c: values of potential outcomes

getEst<-function(data.true,data.obs,p,Nv,alphaEst,betaEst,z,O.c)
{
    O.c11 <- O.c[1,]
    O.c10 <- O.c[2,]
    O.c01 <- O.c[3,]
    O.c00 <- O.c[4,]
    # true exposure level
    my.prod.true <- z%*%data.true
    my.I.true <- as.numeric(my.prod.true > 0)
    c11.true <- z*my.I.true
    c10.true <- z*(1-my.I.true)
    c01.true <- (1-z)*my.I.true
    c00.true <- (1-z)*(1-my.I.true)
    # observed exposure level
    my.prod.obs <- z%*%data.obs
    my.I.obs <- as.numeric(my.prod.obs > 0)
    c11.obs <- z*my.I.obs
    c10.obs <- z*(1-my.I.obs)
    c01.obs <- (1-z)*my.I.obs
    c00.obs <- (1-z)*(1-my.I.obs)
    # get c.obs*c.true
    c11.c11 <- c11.obs*c11.true
    c11.c10 <- c11.obs*c10.true
    c10.c10 <- c10.obs*c10.true
    c10.c11 <- c10.obs*c11.true
    c01.c01 <- c01.obs*c01.true
    c01.c00 <- c01.obs*c00.true
    c00.c00 <- c00.obs*c00.true
    c00.c01 <- c00.obs*c01.true
    # get A&S estimators in true network
    prob.true <- getProb(data.true,p)
    est.c11.true <- mean(O.c11*c11.true/prob.true[1,])
    est.c10.true <- mean(O.c10*c10.true/prob.true[2,])
    est.c01.true <- mean(O.c01*c01.true/prob.true[3,])
    est.c00.true <- mean(O.c00*c00.true/prob.true[4,])
    est.true <- c(est.c11.true,est.c10.true,est.c01.true,est.c00.true)
    # get A&S estimators in noisy network
    prob.obs <- getProb(data.obs,p)
    est.c11.obs <- sum(O.c11*c11.c11/prob.obs[1,]+O.c10*c11.c10/prob.obs[1,],na.rm = T)/Nv
    est.c10.obs <- sum(O.c10*c10.c10/prob.obs[2,]+O.c11*c10.c11/prob.obs[2,],na.rm = T)/Nv
    est.c01.obs <- sum(O.c01*c01.c01/prob.obs[3,]+O.c00*c01.c00/prob.obs[3,],na.rm = T)/Nv
    est.c00.obs <- sum(O.c00*c00.c00/prob.obs[4,]+O.c01*c00.c01/prob.obs[4,],na.rm = T)/Nv
    est.obs <- c(est.c11.obs,est.c10.obs,est.c01.obs,est.c00.obs)
    # get MME estimators in noisy network
    if (alphaEst>1 | betaEst >1) {
      est.new <- rep(0,4)
    } else{
      weight <- getWeight(data.obs,p,Nv,alphaEst,betaEst)
      ind <- which(weight[1,]==0)
      est.c11.new <- sum(O.c11*c11.c11*weight[1,]+O.c10*c11.c10*weight[1,]+O.c10*c10.c10*weight[2,]+O.c11*c10.c11*weight[2,],na.rm = T)/Nv+sum(O.c11[ind]*c11.c11[ind]/prob.obs[1,ind]+O.c10[ind]*c11.c10[ind]/prob.obs[1,ind],na.rm = T)/Nv
      est.c10.new <- sum(O.c10*c10.c10*weight[4,]+O.c11*c10.c11*weight[4,]+O.c11*c11.c11*weight[3,]+O.c10*c11.c10*weight[3,],na.rm = T)/Nv+sum(O.c10[ind]*c10.c10[ind]/prob.obs[2,ind]+O.c11[ind]*c10.c11[ind]/prob.obs[2,ind],na.rm = T)/Nv
      est.c01.new <- sum(O.c01*c01.c01*weight[5,]+O.c00*c01.c00*weight[5,]+O.c00*c00.c00*weight[6,]+O.c01*c00.c01*weight[6,],na.rm = T)/Nv+sum(O.c01[ind]*c01.c01[ind]/prob.obs[3,ind]+O.c00[ind]*c01.c00[ind]/prob.obs[3,ind],na.rm = T)/Nv
      est.c00.new <- sum(O.c00*c00.c00*weight[8,]+O.c01*c00.c01*weight[8,]+O.c01*c01.c01*weight[7,]+O.c00*c01.c00*weight[7,],na.rm = T)/Nv+sum(O.c00[ind]*c00.c00[ind]/prob.obs[4,ind]+O.c01[ind]*c00.c01[ind]/prob.obs[4,ind],na.rm = T)/Nv
      est.new <- c(est.c11.new,est.c10.new,est.c01.new,est.c00.new)
    }
    est <- c(est.true,est.obs,est.new)
    # get A&S estimators variance in true network
    prob.true.joint <- getProbJoint(as.matrix(data.true),p)
    est.c11.true.var <- sum(c11.true*(1-prob.true[1,])*O.c11^2/prob.true[1,]^2)/Nv^2+
      sum((O.c11*c11.true/prob.true[1,])%o%(O.c11*c11.true/prob.true[1,])*(1-prob.true[1,]%o%prob.true[1,]/prob.true.joint[[1]])*(1-diag(Nv))*(prob.true.joint[[1]]>0) ,na.rm = T)/Nv^2+
      sum(outer(O.c11^2*c11.true/prob.true[1,]/2,O.c11^2*c11.true/prob.true[1,]/2,FUN = function(x,y) x+y)*(1-diag(Nv))*(prob.true.joint[[1]]==0),na.rm = T)/Nv^2
    est.c10.true.var <- sum(c10.true*(1-prob.true[2,])*O.c10^2/prob.true[2,]^2)/Nv^2+
      sum((O.c10*c10.true/prob.true[2,])%o%(O.c10*c10.true/prob.true[2,])*(1-prob.true[2,]%o%prob.true[2,]/prob.true.joint[[2]])*(1-diag(Nv))*(prob.true.joint[[2]]>0) ,na.rm = T)/Nv^2+
      sum(outer(O.c10^2*c10.true/prob.true[2,]/2,O.c10^2*c10.true/prob.true[2,]/2,FUN = function(x,y) x+y)*(1-diag(Nv))*(prob.true.joint[[2]]==0),na.rm = T)/Nv^2
    est.c01.true.var <- sum(c01.true*(1-prob.true[3,])*O.c01^2/prob.true[3,]^2)/Nv^2+
      sum((O.c01*c01.true/prob.true[3,])%o%(O.c01*c01.true/prob.true[3,])*(1-prob.true[3,]%o%prob.true[3,]/prob.true.joint[[3]])*(1-diag(Nv))*(prob.true.joint[[3]]>0) ,na.rm = T)/Nv^2+
      sum(outer(O.c01^2*c01.true/prob.true[3,]/2,O.c01^2*c01.true/prob.true[3,]/2,FUN = function(x,y) x+y)*(1-diag(Nv))*(prob.true.joint[[3]]==0),na.rm = T)/Nv^2
    est.c00.true.var <- sum(c00.true*(1-prob.true[4,])*O.c00^2/prob.true[4,]^2)/Nv^2+
      sum((O.c00*c00.true/prob.true[4,])%o%(O.c00*c00.true/prob.true[4,])*(1-prob.true[4,]%o%prob.true[4,]/prob.true.joint[[4]])*(1-diag(Nv))*(prob.true.joint[[4]]>0) ,na.rm = T)/Nv^2+
      sum(outer(O.c00^2*c00.true/prob.true[4,]/2,O.c00^2*c00.true/prob.true[4,]/2,FUN = function(x,y) x+y)*(1-diag(Nv))*(prob.true.joint[[4]]==0),na.rm = T)/Nv^2
    est.true.var <- c(est.c11.true.var,est.c10.true.var,est.c01.true.var,est.c00.true.var)
    
    # get MME estimators variance in noisy network
    if (alphaEst>1 | betaEst >1) {
      est.new.var <- rep(0,4)
    } else{
      prob_joint= getProbJointTrueObs(data.obs,p,Nv,alphaEst,betaEst)
      prob_joint0=prob_joint$prob_joint0
      prob_joint1=prob_joint$prob_joint1
      
      est.c11.new.var.tmp3 = est.c10.new.var.tmp3 = est.c01.new.var.tmp3 = est.c00.new.var.tmp3 = matrix(0,Nv,Nv)
          
      est.c11.new.var.tmp1 <- O.c11*c11.c11*weight[1,]+O.c10*c11.c10*weight[1,]+O.c10*c10.c10*weight[2,]+O.c11*c10.c11*weight[2,]
      est.c11.new.var.tmp2 <- O.c11^2*c11.c11*weight[1,]+O.c10^2*c11.c10*weight[1,]+O.c10^2*c10.c10*weight[2,]+O.c11^2*c10.c11*weight[2,]
      est.c10.new.var.tmp1 <- O.c10*c10.c10*weight[4,]+O.c11*c10.c11*weight[4,]+O.c11*c11.c11*weight[3,]+O.c10*c11.c10*weight[3,]
      est.c10.new.var.tmp2 <- O.c10^2*c10.c10*weight[4,]+O.c11^2*c10.c11*weight[4,]+O.c11^2*c11.c11*weight[3,]+O.c10^2*c11.c10*weight[3,]
      est.c01.new.var.tmp1 <- O.c01*c01.c01*weight[5,]+O.c00*c01.c00*weight[5,]+O.c00*c00.c00*weight[6,]+O.c01*c00.c01*weight[6,]
      est.c01.new.var.tmp2 <- O.c01^2*c01.c01*weight[5,]+O.c00^2*c01.c00*weight[5,]+O.c00^2*c00.c00*weight[6,]+O.c01^2*c00.c01*weight[6,]
      est.c00.new.var.tmp1 <- O.c00*c00.c00*weight[8,]+O.c01*c00.c01*weight[8,]+O.c01*c01.c01*weight[7,]+O.c00*c01.c00*weight[7,]
      est.c00.new.var.tmp2 <- O.c00^2*c00.c00*weight[8,]+O.c01^2*c00.c01*weight[8,]+O.c01^2*c01.c01*weight[7,]+O.c00^2*c01.c00*weight[7,]
      
      degree.obs <- colSums(data.obs)
      degree.obs.sum <- outer(degree.obs,degree.obs,"+")
      degree.obs.i <- matrix(degree.obs,Nv,Nv)
      degree.obs.j <- matrix(degree.obs,Nv,Nv,byrow = TRUE)
      C.obs <- data.obs %*% data.obs  
      C.hat <- (C.obs-alphaEst*(1-alphaEst-betaEst)*(degree.obs.sum-2*data.obs)+alphaEst^2*(Nv-2))/(1-alphaEst-betaEst)^2
      
      A.hat <- (data.obs-alphaEst)/(1-alphaEst-betaEst)
      data.true = as.matrix(data.true)
      for (i in setdiff(1:(Nv-1),ind)) {
        for (j in setdiff(((i+1):Nv),ind)) {
          if(C.hat[i,j]>0 & (!is.null(j))) {
            prob_joint1_inv = solve((matrix(prob_joint1[i,j,],4)))
            prob_joint0_inv = solve((matrix(prob_joint0[i,j,],4)))
 
 
            est.c11.new.var.tmp3[i,j] = est.c11.new.var.tmp1[i]*est.c11.new.var.tmp1[j]+(1-A.hat[i,j])*(-
                                                                                                          prob_joint1_inv[1,1]*(O.c11[i]*c11.c11[i]+O.c10[i]*c11.c10[i])*(O.c11[j]*c11.c11[j]+O.c10[j]*c11.c10[j])-
                                                                                                          prob_joint1_inv[1,2]*(O.c11[i]*c11.c11[i]+O.c10[i]*c11.c10[i])*(O.c10[j]*c10.c10[j]+O.c11[j]*c10.c11[j])-
                                                                                                          prob_joint1_inv[1,3]*(O.c10[i]*c10.c10[i]+O.c11[i]*c10.c11[i])*(O.c11[j]*c11.c11[j]+O.c10[j]*c11.c10[j])-
                                                                                                          prob_joint1_inv[1,4]*(O.c10[i]*c10.c10[i]+O.c11[i]*c10.c11[i])*(O.c10[j]*c10.c10[j]+O.c11[j]*c10.c11[j]))+
              A.hat[i,j]*(- (O.c11[i]*c11.c11[i]+O.c10[i]*c11.c10[i])*(O.c11[j]*c11.c11[j]+O.c10[j]*c11.c10[j])/((p^2*(1-(1-p)^(degree.obs[j]-data.obs[i,j])-(1-p)^(degree.obs[i]-data.obs[i,j])+(1-p)^(degree.obs[j]+degree.obs[i]-2*data.obs[i,j]-C.obs[i,j])) )*betaEst+(1-betaEst)*p^2 ) )
            
            
            est.c10.new.var.tmp3[i,j] = est.c10.new.var.tmp1[i]*est.c10.new.var.tmp1[j]+(1-A.hat[i,j])*(-
                                                                                                          prob_joint1_inv[4,1]*(O.c11[i]*c11.c11[i]+O.c10[i]*c11.c10[i])*(O.c11[j]*c11.c11[j]+O.c10[j]*c11.c10[j])-
                                                                                                          prob_joint1_inv[4,2]*(O.c11[i]*c11.c11[i]+O.c10[i]*c11.c10[i])*(O.c10[j]*c10.c10[j]+O.c11[j]*c10.c11[j])-
                                                                                                          prob_joint1_inv[4,3]*(O.c10[i]*c10.c10[i]+O.c11[i]*c10.c11[i])*(O.c11[j]*c11.c11[j]+O.c10[j]*c11.c10[j])-
                                                                                                          prob_joint1_inv[4,4]*(O.c10[i]*c10.c10[i]+O.c11[i]*c10.c11[i])*(O.c10[j]*c10.c10[j]+O.c11[j]*c10.c11[j]) )+
              A.hat[i,j]*(est.c10.new.var.tmp2[i]/2+ est.c10.new.var.tmp2[j]/2)
            
            
            est.c01.new.var.tmp3[i,j] = est.c01.new.var.tmp1[i]*est.c01.new.var.tmp1[j]-
              prob_joint0_inv[1,1]*(O.c01[i]*c01.c01[i]+O.c00[i]*c01.c00[i])*(O.c01[j]*c01.c01[j]+O.c00[j]*c01.c00[j])-
              prob_joint0_inv[1,2]*(O.c01[i]*c01.c01[i]+O.c00[i]*c01.c00[i])*(O.c00[j]*c00.c00[j]+O.c01[j]*c00.c01[j])-
              prob_joint0_inv[1,3]*(O.c00[i]*c00.c00[i]+O.c01[i]*c00.c01[i])*(O.c01[j]*c01.c01[j]+O.c00[j]*c01.c00[j])-
              prob_joint0_inv[1,4]*(O.c00[i]*c00.c00[i]+O.c01[i]*c00.c01[i])*(O.c00[j]*c00.c00[j]+O.c01[j]*c00.c01[j])
            
            est.c00.new.var.tmp3[i,j] = est.c00.new.var.tmp1[i]*est.c00.new.var.tmp1[j]-
              prob_joint0_inv[4,1]*(O.c01[i]*c01.c01[i]+O.c00[i]*c01.c00[i])*(O.c01[j]*c01.c01[j]+O.c00[j]*c01.c00[j])-
              prob_joint0_inv[4,2]*(O.c01[i]*c01.c01[i]+O.c00[i]*c01.c00[i])*(O.c00[j]*c00.c00[j]+O.c01[j]*c00.c01[j])-
              prob_joint0_inv[4,3]*(O.c00[i]*c00.c00[i]+O.c01[i]*c00.c01[i])*(O.c01[j]*c01.c01[j]+O.c00[j]*c01.c00[j])-
              prob_joint0_inv[4,4]*(O.c00[i]*c00.c00[i]+O.c01[i]*c00.c01[i])*(O.c00[j]*c00.c00[j]+O.c01[j]*c00.c01[j])
            
          }
        }
      }
      est.c11.new.var.tmp3 = est.c11.new.var.tmp3+t(est.c11.new.var.tmp3)
      est.c10.new.var.tmp3 = est.c10.new.var.tmp3+t(est.c10.new.var.tmp3)
      est.c01.new.var.tmp3 = est.c01.new.var.tmp3+t(est.c01.new.var.tmp3)
      est.c00.new.var.tmp3 = est.c00.new.var.tmp3+t(est.c00.new.var.tmp3)
       
      est.c11.new.var <-  sum( est.c11.new.var.tmp1^2-est.c11.new.var.tmp2,na.rm = T )/Nv^2 + sum(est.c11.new.var.tmp3,na.rm = T)/Nv^2
      est.c10.new.var <-  sum( est.c10.new.var.tmp1^2-est.c10.new.var.tmp2,na.rm = T )/Nv^2+sum(est.c10.new.var.tmp3,na.rm = T)/Nv^2
      est.c01.new.var <- sum(est.c01.new.var.tmp1^2-est.c01.new.var.tmp2,na.rm = T )/Nv^2+sum(est.c01.new.var.tmp3,na.rm = T)/Nv^2
      est.c00.new.var <-  sum(est.c00.new.var.tmp1^2-est.c00.new.var.tmp2 ,na.rm = T)/Nv^2+sum(est.c00.new.var.tmp3,na.rm = T)/Nv^2
      
      est.new.var <- c(est.c11.new.var,est.c10.new.var,est.c01.new.var,est.c00.new.var)
    }
    return(c(est,sqrt(est.true.var),sqrt(est.new.var)))
}

# Get exposure probabilities for the four-level exposure model with Bernoulli  random assignment of treatment
getProb<-function(data,p)
{
    g.degree <- colSums(data)
    prob.c11 <- p*(1-(1-p)^g.degree)
    prob.c10 <- p*(1-p)^g.degree
    prob.c01 <- (1-p)*(1-(1-p)^g.degree)
    prob.c00 <- (1-p)^(g.degree+1)
    prob <- rbind(prob.c11,prob.c10,prob.c01,prob.c00)
    return(prob)
}
getProbJoint<-function(data,p)
{
    g.degree <- colSums(data)
    Nv <- ncol(data)
    common.neigh <- data %*% data
    prob.c11 <-  prob.c10 <-  prob.c01 <- prob.c00 <- matrix(0,nrow = Nv,ncol=Nv)
    for (i in 1:Nv) {
        for (j in 1:Nv) {
            prob.c11[i,j] <- ifelse(data[i,j]==1,p^2,p^2*(1-(1-p)^g.degree[i]-(1-p)^g.degree[j] +(1-p)^(g.degree[i]+g.degree[j]-common.neigh[i,j]) ))
            prob.c10[i,j] <- ifelse(data[i,j]==1,0,p^2*(1-p)^(g.degree[i]+g.degree[j]-common.neigh[i,j]) )
            prob.c01[i,j] <- (1-p)^2*(1-(1-p)^(g.degree[i]-data[i,j])-(1-p)^(g.degree[j]-data[i,j])+
                                          (1-p)^(g.degree[i]+g.degree[j]-common.neigh[i,j]-2*data[i,j]))
            prob.c00[i,j] <- (1-p)^(2+g.degree[i]+g.degree[j]-common.neigh[i,j]-2*data[i,j])
        }
    }
    prob <- list(prob.c11,prob.c10,prob.c01,prob.c00)
    return(prob)
}


# Get weights for MME estimators in noisy networks
getWeight<-function(data.obs,p,Nv,alpha,beta)
{
    degree.obs <- colSums(data.obs)
    degree.hat <- (degree.obs-(Nv-1)*alpha)/(1-alpha-beta)
    ind <- which(degree.hat<1)
    b <- (1-(1-alpha*p)^(Nv-1-degree.hat))*(1-p)^degree.hat
    e <- (1-alpha*p)^(Nv-1-degree.hat)*(1-p)^degree.hat
    s <- (1-alpha*p)^(Nv-1-degree.hat)*(1-(1-beta)*p)^degree.hat-e
    a <- 1-(1-p)^degree.hat-s
    down <- a*e-b*s
    weight <- matrix(0,nrow=8,ncol=Nv)
    weight[1,] <- e/(down*p)
    weight[2,] <- -b/(down*p)
    weight[3,] <- -s/(down*p)
    weight[4,] <- a/(down*p)
    weight[5,]  <- weight[1,]*p/(1-p)
    weight[6,]  <- weight[2,]*p/(1-p)
    weight[7,]  <- weight[3,]*p/(1-p)
    weight[8,]  <- weight[4,]*p/(1-p)
    weight[,ind] <- 0
    return(weight)
}

# Get estimators of the expected confusion matrices B and H
getProbJointTrueObs<-function(data.obs,p,Nv,alpha,beta)
{
  prob_joint1 <- prob_joint0 <- array(0,c(Nv,Nv,16))
  degree.obs <- colSums(data.obs)
  degree.obs.sum <- outer(degree.obs,degree.obs,"+")
  degree.obs.i <- matrix(degree.obs,Nv,Nv)
  degree.obs.j <- matrix(degree.obs,Nv,Nv,byrow = TRUE)
  degree.hat <- (degree.obs-(Nv-1)*alpha)/(1-alpha-beta)
  degree.hat.sum <- outer(degree.hat,degree.hat,"+")
  degree.hat.i <- matrix(degree.hat,Nv,Nv)
  degree.hat.j <- matrix(degree.hat,Nv,Nv,byrow = TRUE)
  A.hat <- (data.obs-alpha)/(1-alpha-beta)
  
  C.obs <- data.obs %*% data.obs
  degree.obs.sum <- outer(degree.obs,degree.obs,"+")
  C.hat <- (C.obs-alpha*(1-alpha-beta)*(degree.obs.sum-2*data.obs)+alpha^2*(Nv-2))/(1-alpha-beta)^2
  
  prob_joint30 <- array(0,c(Nv,Nv,8))
  prob_joint30[,,1] <- (1-p)^(Nv)*(1+p/(1-p)*(1-alpha))^(Nv-2-degree.hat.sum+C.hat+2*A.hat)
  prob_joint30[,,2] <- (1-p)^(Nv)*(1+p/(1-p)*(1-alpha)^2)^(Nv-2-degree.hat.sum+C.hat+2*A.hat)*(1+p/(1-p)*(1-alpha)*beta)^(degree.hat.i-C.hat-A.hat)
  prob_joint30[,,3] <- prob_joint30[,,1]
  prob_joint30[,,4] <- (1-p)^(Nv)*(1+p/(1-p)*(1-alpha)^2)^(Nv-2-degree.hat.sum+C.hat+2*A.hat)*(1+p/(1-p)*(1-alpha)*beta)^(degree.hat.j-C.hat-A.hat)
  prob_joint30[,,5] <- (1-p)^(Nv)*(1+p/(1-p)*(1-alpha))^(Nv-2-degree.hat.sum+C.hat+2*A.hat)*(1+p/(1-p)*beta)^(degree.hat.i-C.hat-A.hat) - prob_joint30[,,2]
  prob_joint30[,,6] <- (1-p)^(Nv)*(1+p/(1-p)*(1-alpha))^(Nv-2-degree.hat.i+A.hat)- prob_joint30[,,4]
  prob_joint30[,,7] <- (1-p)^(Nv)*(1+p/(1-p)*(1-alpha))^(Nv-2-degree.hat.j+A.hat)- prob_joint30[,,2]
  prob_joint30[,,8] <- (1-p)^(Nv)*(1+p/(1-p)*(1-alpha))^(Nv-2-degree.hat.sum+C.hat+2*A.hat)*(1+p/(1-p)*beta)^(degree.hat.j-C.hat-A.hat) - prob_joint30[,,4]
 
  prob_00_00_obs =(1-p)^2*(1-p)^(degree.obs.sum-C.obs-2*data.obs)
  prob_00_01_obs =(1-p)^2*(1-p)^(degree.obs.i-data.obs)*(1-(1-p)^(degree.obs.j-C.obs-data.obs))
  prob_01_00_obs =(1-p)^2*(1-p)^(degree.obs.j-data.obs)*(1-(1-p)^(degree.obs.i-C.obs-data.obs))
  prob_01_01_obs =(1-p)^2*(1-(1-p)^(degree.obs.j-data.obs)-(1-p)^(degree.obs.i-data.obs)+(1-p)^(degree.obs.sum-C.obs-2*data.obs))
  
  prob_00_00_hat =(1-p)^2*(1-p)^(degree.hat.sum-C.hat-2*A.hat)
  prob_00_01_hat =(1-p)^2*(1-p)^(degree.hat.i-A.hat)*(1-(1-p)^(degree.hat.j-C.hat-A.hat))
  prob_01_00_hat =(1-p)^2*(1-p)^(degree.hat.j-A.hat)*(1-(1-p)^(degree.hat.i-C.hat-A.hat))
  prob_01_01_hat =(1-p)^2*(1-(1-p)^(degree.hat.j-A.hat)-(1-p)^(degree.hat.i-A.hat)+(1-p)^(degree.hat.sum-C.hat-2*A.hat))
  
  
  prob_joint0[,,16] = (1-p)^(Nv)*(1+p/(1-p)*(1-alpha)^2)^(Nv-2-degree.hat.sum+C.hat+2*A.hat)
  prob_joint0[,,12] =prob_joint30[,,4]-prob_joint0[,,16]
  prob_joint0[,,8] = prob_joint30[,,2]-prob_joint0[,,16]
  prob_joint0[,,4] =  prob_00_00_obs - prob_joint0[,,8]- prob_joint0[,,12] -prob_joint0[,,16]
  
  prob_joint0[,,15] = prob_joint30[,,3] -prob_joint0[,,16]
  prob_joint0[,,11] = prob_joint30[,,6] -prob_joint0[,,15]
  prob_joint0[,,7] = prob_joint30[,,5] -prob_joint0[,,15]
  prob_joint0[,,3] = prob_00_01_obs - prob_joint0[,,7] - prob_joint0[,,11]-prob_joint0[,,15]
  
  prob_joint0[,,14] = prob_joint30[,,1] -prob_joint0[,,16]
  prob_joint0[,,10] = prob_joint30[,,8] -prob_joint0[,,14]
  prob_joint0[,,6] = prob_joint30[,,7] -prob_joint0[,,14]
  prob_joint0[,,2] = prob_01_00_obs - prob_joint0[,,6]- prob_joint0[,,10]- prob_joint0[,,14]
  
  prob_joint0[,,13] = prob_00_00_hat - prob_joint0[,,14] -prob_joint0[,,15] -prob_joint0[,,16]
  prob_joint0[,,9] = prob_00_01_hat - prob_joint0[,,10] -prob_joint0[,,11] -prob_joint0[,,12]
  prob_joint0[,,5] = prob_01_00_hat - prob_joint0[,,6] -prob_joint0[,,7] -prob_joint0[,,8]
  prob_joint0[,,1] = prob_01_01_hat - prob_joint0[,,2] -prob_joint0[,,3] -prob_joint0[,,4]
  
  prob_joint1=prob_joint0/(1-p)^2*p^2*(1-alpha)
  prob_joint1[,,1] = prob_joint1[,,1] + alpha*prob_01_01_hat/(1-p)^2*p^2
  prob_joint1[,,5] = prob_joint1[,,5] + alpha*prob_01_00_hat/(1-p)^2*p^2
  prob_joint1[,,9] = prob_joint1[,,9] + alpha*prob_00_01_hat/(1-p)^2*p^2
  prob_joint1[,,13] = prob_joint1[,,13] + alpha*prob_00_00_hat/(1-p)^2*p^2
  
  return(list(prob_joint0=prob_joint0,prob_joint1=prob_joint1) )
}
