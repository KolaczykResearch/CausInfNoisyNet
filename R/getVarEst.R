# This program compute varinces and confidence intervals of 2 types of estimators: A&S estimators in the noisy network and MME in the noisy network.

# data.true: adjacency matrix of the true network
# data.obs1, data.obs2, data.obs3: adjacency matrix of the noisy network
# p: the probability for Bernoulli random assignment of treatment
# Nv: number of nodes
# alphaEst: point-estimate for alpha (Type I error)
# betaEst: point-estimate for beta (Type II error)
# z: treatment assignment vector
# O.c: values of potential outcomes
# result: A&S estimators in the noisy network and MME in the noisy network
# Nb2: number of bootstrapping
# Zcv: Z critical value

getVarEst<-function(data.true,data.obs1,data.obs2,data.obs3,p,Nv,alphaEst,betaEst,z,O.c,result,Nb2,Zcv)
{
    est_ck <- getCk(data.true,data.obs1,data.obs2,data.obs3,p,Nv,alphaEst,betaEst,z,O.c)
    est_boot <- replicate(Nb2,getVarOne(data.true,data.obs1,data.obs2,data.obs3,p,Nv,alphaEst,betaEst,z,O.c,est_ck))
    res <- apply(est_boot,1,function(x) sd(x) )
    res_mean <- apply(est_boot,1,function(x) mean(x) )
    res <- rbind(res,2*result[5:12]-res_mean-Zcv*res,2*result[5:12]-res_mean+Zcv*res)
    RF_normal <- O.c>res[2,] & O.c<res[3,]
    return(c(res[1,],RF_normal))
}

# get a bootstrap resample adjacency matrix
getBootOne <- function(Nv, data.obs1, data.obs2, data.obs3)
{
    Y_boot <- data.obs1
    sample_ind <- sample(3,Nv*Nv,replace = T)
    ind2 <- sample_ind==2
    ind3 <- sample_ind==3
    Y_boot[ind2] <- data.obs2[ind2]
    Y_boot[ind3] <- data.obs3[ind3]
    Y_boot <-  Matrix::forceSymmetric(Y_boot,uplo="U")
    return(Y_boot)
}

# get an empirical cumulative distribution function
getcdf <- function(x_seq)
{
    x <- sort(unique(x_seq))
    if(length(x)<1) {cdf.fun <- function(y){runif(1)}
     } else { cdf.fun <- ecdf(x_seq) }
    return(cdf.fun)
}

# get an inverse empirical cumulative distribution function
getinvcdf <- function(x_seq,prob) {
  x <- sort(unique(x_seq))
  if(length(x)<1) { res <- NA 
  } else if (length(x)==1) {
    res <- x
  } else {
    e_cdf <- c(sum(x_seq==x[1])/length(x_seq),1)
    res <- x[which(e_cdf >= prob)[1]]
  }
  return(res)
}

# impute potential values for each individual
getCk<-function(data.true,data.obs1,data.obs2,data.obs3,p,Nv,alphaEst,betaEst,z,O.c)
{
    O.c11 <- O.c[1]
    O.c10 <- O.c[2]
    O.c01 <- O.c[3]
    O.c00 <- O.c[4]
    # true exposure level
    my.prod.true <- z%*%data.true
    my.I.true <- as.numeric(my.prod.true > 0)
    c11.true <- z*my.I.true
    c10.true <- z*(1-my.I.true)
    c01.true <- (1-z)*my.I.true
    c00.true <- (1-z)*(1-my.I.true)
    # observed exposure level
    my.prod.obs <- z%*%data.obs1
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
    # get obs y_i(c_k)
    y.obs.c11 <- O.c11*c11.c11+O.c10*c11.c10
    y.obs.c10 <- O.c10*c10.c10+O.c11*c10.c11
    y.obs.c01 <- O.c01*c01.c01+O.c00*c01.c00
    y.obs.c00 <- O.c00*c00.c00+O.c01*c00.c01
    y.obs.c11_seq <- y.obs.c11[y.obs.c11!=0]
    y.obs.c10_seq <- y.obs.c10[y.obs.c10!=0]
    y.obs.c01_seq <- y.obs.c01[y.obs.c01!=0]
    y.obs.c00_seq <- y.obs.c00[y.obs.c00!=0]  
    y.obs.ck <- rbind(y.obs.c11,y.obs.c10,y.obs.c01,y.obs.c00)
    y.obs.ck.seq <- list(y.obs.c11_seq,y.obs.c10_seq,y.obs.c01_seq,y.obs.c00_seq)
    # get cdf of obs y_i(c_k)
    cdf.c11 <- getcdf(y.obs.c11_seq)
    cdf.c10 <- getcdf(y.obs.c10_seq)
    cdf.c01 <- getcdf(y.obs.c01_seq)
    cdf.c00 <- getcdf(y.obs.c00_seq)
    cdf.ck <- list(cdf.c11,cdf.c10,cdf.c01,cdf.c00)
    # get boot y_i(c_k)
    y.boot.ck <- y.obs.ck 
    for (i in 1:4) { 
        ind <- which(y.obs.ck[i,]==0)
        for (j in ind) {
                ind_class <- which(y.obs.ck[,j]!=0)
                y.boot.ck[i,j] <- getinvcdf(y.obs.ck.seq[[i]],prob=cdf.ck[[ind_class]](y.obs.ck[ind_class,j]) ) 
                if(is.na(y.boot.ck[i,j]))  y.boot.ck[i,j] <- sum(y.obs.ck[,j])
        }
    }
  
    return(y.boot.ck)
}


# compute A&S estimators and MME in a bootstrap resample network
getVarOne<-function(data.true,data.obs1,data.obs2,data.obs3,p,Nv,alphaEst,betaEst,z,O.c,est_ck)
{
    data.boot <- getBootOne(Nv, data.obs1,data.obs2,data.obs3)
    z.boot <- rbinom(Nv, 1, p)
    O.c11 <- O.c[1]
    O.c10 <- O.c[2]
    O.c01 <- O.c[3]
    O.c00 <- O.c[4]
    # true exposure level
    my.prod.true <- z%*%data.true
    my.I.true <- as.numeric(my.prod.true > 0)
    c11.true <- z*my.I.true
    c10.true <- z*(1-my.I.true)
    c01.true <- (1-z)*my.I.true
    c00.true <- (1-z)*(1-my.I.true)
    # observed exposure level
    my.prod.obs <- z%*%data.obs1
    my.I.obs <- as.numeric(my.prod.obs > 0)
    c11.obs <- z*my.I.obs
    c10.obs <- z*(1-my.I.obs)
    c01.obs <- (1-z)*my.I.obs
    c00.obs <- (1-z)*(1-my.I.obs)
    # bootstrap exposure level
    my.prod.boot <- z.boot%*%data.boot
    my.I.boot <- as.numeric(my.prod.boot > 0)
    c11.boot <- z.boot*my.I.boot
    c10.boot <- z.boot*(1-my.I.boot)
    c01.boot <- (1-z.boot)*my.I.boot
    c00.boot <- (1-z.boot)*(1-my.I.boot)
    ck.boot <- rbind(c11.boot,c10.boot,c01.boot,c00.boot)
    # get boot y_i(c_k)
    y.boot.ck <- ck.boot*est_ck
    y.boot.c11 <- y.boot.ck[1,]
    y.boot.c10 <- y.boot.ck[2,]
    y.boot.c01 <- y.boot.ck[3,]
    y.boot.c00 <- y.boot.ck[4,]
    # get MME estimators in noisy network
    if (alphaEst>1 | betaEst >1) { 
        est.new <- rep(0,4) 
    } else{
        prob.boot <- getProb(data.boot,p)
        weight <- getWeight(data.boot,p,Nv,alphaEst,betaEst)
        ind <- which(weight[1,]==0)
        est.c11.new <- sum(y.boot.c11*weight[1,]+y.boot.c10*weight[2,],na.rm = T)/Nv+sum(y.boot.c11[ind]/prob.boot[1,ind],na.rm = T)/Nv
        est.c10.new <- sum(y.boot.c10*weight[4,]+y.boot.c11*weight[3,],na.rm = T)/Nv+sum(y.boot.c10[ind]/prob.boot[2,ind],na.rm = T)/Nv
        est.c01.new <- sum(y.boot.c01*weight[5,]+y.boot.c00*weight[6,],na.rm = T)/Nv+sum(y.boot.c01[ind]/prob.boot[3,ind],na.rm = T)/Nv
        est.c00.new <- sum(y.boot.c00*weight[8,]+y.boot.c01*weight[7,],na.rm = T)/Nv+sum(y.boot.c00[ind]/prob.boot[4,ind],na.rm = T)/Nv
        est.new <- c(est.c11.new,est.c10.new,est.c01.new,est.c00.new) 
    }
    est.c11.obs <- sum(y.boot.c11/prob.boot[1,],na.rm = T)/Nv
    est.c10.obs <- sum(y.boot.c10/prob.boot[2,],na.rm = T)/Nv
    est.c01.obs <- sum(y.boot.c01/prob.boot[3,],na.rm = T)/Nv
    est.c00.obs <- sum(y.boot.c00/prob.boot[4,],na.rm = T)/Nv
    est.obs <- c(est.c11.obs,est.c10.obs,est.c01.obs,est.c00.obs)
    return(c(est.obs,est.new))
}





 
