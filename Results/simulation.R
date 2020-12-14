library(igraph)
library(Matrix)
library(ggplot2)
library(dplyr)

################ functions ##################
 
# get estimators for a Monte Carlo simulation trial
getOne <- function(data.network,alpha,beta,p,O.c,z,alpha0, Nb, tiny)
{
    data.true <- get.adjacency(data.network)
    Nv <- ncol(data.true)
    data.obs1 <- as.matrix(getNoisyNet(data.true,Nv,alpha,beta) )
    data.obs2 <- as.matrix(getNoisyNet(data.true,Nv,alpha,beta) )
    data.obs3 <- as.matrix(getNoisyNet(data.true,Nv,alpha,beta) )
    NetNoise <- getNetNoise(alpha0,Nv, data.obs1, data.obs2, data.obs3,Nb,tiny)
    alphaEst <- NetNoise[1]
    betaEst <- NetNoise[2]
    result <- getEst(data.true,data.obs1,p,Nv,alphaEst,betaEst,z,O.c)
    return(result)
}

# get 95% confidence intervals and mean for biases and standard deviations of estimators 
getCI<-function(boot_case,boot_index,boot_n,O.c)
{
  boot_mean <- sapply(1:boot_n,function(i) apply(boot_case[boot_index[i,],],2,mean))
  boot_sd <- sapply(1:boot_n,function(i) apply(boot_case[boot_index[i,],],2,sd))
  result <- matrix(c(apply(boot_mean,1,mean)-rep(O.c,3),apply(boot_mean,1,quantile,probs=0.025)-rep(O.c,3),
                     apply(boot_mean,1,quantile,probs=0.975)-rep(O.c,3),apply(boot_sd,1,mean),
                     apply(boot_sd,1,quantile,probs=0.025),apply(boot_sd,1,quantile,probs=0.975)),nrow=4)
  return(result)
}


# transform data to desired format for plotting
getRes <- function(school,alpha,beta,res_boot)
{
  cbind(rep(school,12),rep(alpha,12),rep(beta,12),
        rbind(cbind(rep(1,4),seq(1,4),res_boot[,c(1,4,7)]),cbind(rep(2,4),seq(1,4),res_boot[,c(2,5,8)]),cbind(rep(3,4),seq(1,4),res_boot[,c(3,6,9)])))
}
 
source("../../R/getNoisyNet.R")
source("../../R/getNetNoise.R")
source("../../R/getEst.R")
################ load data ##################
indir <-  "../../Data/"
load( paste0(indir,"dataTrue.RData") )
load( paste0(indir,"z.RData") )

################ simulation ##################
# number of Monte Carlo simulation trials
m <- 10^4

# probability of treatment
p <- 0.1

# values of potential outcomes
O.c<- c(10,7,5,1)

# number of bootstrap resampling
boot_n <- 10^3

# parameters for algorithm 1 (estimate alpha and beta)
alpha0 <- 0.001
Nb <- 500 
tiny <- 0.0001

# set alpha and beta 
beta <- rep(c(0.05,0.1,0.15),2)
alpha <- rep(c(0.005,0.01),each=3)
 
# estimators for 4 schools
res_est1 <- res_est2 <- res_est3 <- res_est4 <- list()
for (i in 1:length(beta)) {
  res_est1[[i]] <- t(sapply(1:m, function(j) getOne(school1_undiCol_true,alpha[i],beta[i],p,O.c,z1[j,],alpha0, Nb, tiny)))
  res_est2[[i]] <- t(sapply(1:m, function(j) getOne(school2_undiCol_true,alpha[i],beta[i],p,O.c,z2[j,],alpha0, Nb, tiny)))
  res_est3[[i]] <- t(sapply(1:m, function(j) getOne(school3_undiCol_true,alpha[i],beta[i],p,O.c,z3[j,],alpha0, Nb, tiny)))
  res_est4[[i]] <- t(sapply(1:m, function(j) getOne(school4_undiCol_true,alpha[i],beta[i],p,O.c,z4[j,],alpha0, Nb, tiny)))
}
 
# biases and standard deviations of estimators 
boot_index <- t(replicate(boot_n,sample(1:m,m,replace = T)))
res_boot1 <- lapply(1:length(beta),function (i) getCI(res_est1[[i]],boot_index,boot_n,O.c))
res_boot2 <- lapply(1:length(beta),function (i) getCI(res_est2[[i]],boot_index,boot_n,O.c))
res_boot3 <- lapply(1:length(beta),function (i) getCI(res_est3[[i]],boot_index,boot_n,O.c))
res_boot4 <- lapply(1:length(beta),function (i) getCI(res_est4[[i]],boot_index,boot_n,O.c))
school_bias <- school_sd <- list()
school_bias[[1]] <- Reduce("rbind",lapply(1:length(beta),function (i) getRes(1,alpha[i],beta[i],res_boot1[[i]][,1:9])))
school_bias[[2]] <- Reduce("rbind",lapply(1:length(beta),function (i) getRes(2,alpha[i],beta[i],res_boot2[[i]][,1:9])))
school_bias[[3]] <- Reduce("rbind",lapply(1:length(beta),function (i) getRes(3,alpha[i],beta[i],res_boot3[[i]][,1:9])))
school_bias[[4]] <- Reduce("rbind",lapply(1:length(beta),function (i) getRes(4,alpha[i],beta[i],res_boot4[[i]][,1:9])))
school_bias <- Reduce("rbind",school_bias)
school_sd[[1]] <- Reduce("rbind",lapply(1:length(beta),function (i) getRes(1,alpha[i],beta[i],res_boot1[[i]][,10:18])))
school_sd[[2]] <- Reduce("rbind",lapply(1:length(beta),function (i) getRes(2,alpha[i],beta[i],res_boot2[[i]][,10:18])))
school_sd[[3]] <- Reduce("rbind",lapply(1:length(beta),function (i) getRes(3,alpha[i],beta[i],res_boot3[[i]][,10:18])))
school_sd[[4]] <- Reduce("rbind",lapply(1:length(beta),function (i) getRes(4,alpha[i],beta[i],res_boot4[[i]][,10:18])))
school_sd <- Reduce("rbind",school_sd)
colnames(school_bias) <- colnames(school_sd) <- c("school","alpha","beta","type","level","est","low","upper")
school_bias <- data.frame(school_bias)
school_sd <- data.frame(school_sd)

################ Figure 4.1 ##################
school_bias <- school_bias %>% mutate(level=case_when(level==1 ~ "c11",level==2 ~ "c10",
                                                      level==3 ~ "c01",level==4 ~ "c00") )
Fig4.1 <- school_bias %>% mutate(beta = as.numeric(as.character(beta)) + as.numeric(school)/400 - .0125*(type== 1)+ .0125*(type== 3) +0.00125) %>%
  ggplot(aes(x = beta, y = est,shape = factor(school),color=factor(type),group=factor(type)))+
  geom_point(alpha = 1)  +  
  facet_grid(level~alpha,labeller = label_both) + theme_bw() +    
  scale_x_continuous(breaks = c(0.05, .1, .15), name=expression(beta))+
  scale_y_continuous(name=c("Bias"),limits = c(-1.05,1.05)) +
  theme(plot.title = element_text(hjust = 0.5,size=16),axis.line = element_line(colour = "black"),panel.grid.minor = element_blank(),panel.background = element_blank())+
  scale_shape_manual(name="School",labels =c("School 1","School 2","School 3","School 4"),values=c(15,16,17,18))  +
  scale_color_manual(labels=c("A&S_TrueNetwork","A&S_NoiseNetwork","MME_NoiseNetwork"),name="Estimators",values=c("black", "grey35","grey70"))+ 
  geom_errorbar(aes(ymin=low, ymax=upper,color=factor(type),group=factor(type)), width=.001) 
 
################ Figure 4.2 ##################
school_sd <- school_sd %>% mutate(level=case_when(level==1 ~ "c11",level==2 ~ "c10",
                                                  level==3 ~ "c01",level==4 ~ "c00") )

Fig4.2 <- school_sd %>% mutate(beta = as.numeric(as.character(beta)) + as.numeric(school)/400 - .0125*(type== 1) + .0125*(type== 3) +0.00125) %>%
  ggplot(aes(x = beta, y = est,shape = factor(school),color=factor(type),group=factor(type)))+
  geom_point(alpha = 1)  +  
  facet_grid(level~alpha,labeller = label_both) + theme_bw() +
  scale_x_continuous(breaks = c(0.05, .1, .15), name=expression(beta))+
  scale_y_continuous(name=c("Standard deviation")) +
  theme(plot.title = element_text(hjust = 0.5,size=16),axis.line = element_line(colour = "black"),panel.grid.minor = element_blank(),panel.background = element_blank())+
  scale_shape_manual(name="School",labels =c("School 1","School 2","School 3","School 4"),values=c(15,16,17,18))  +
  scale_color_manual(labels=c("A&S_TrueNetwork","A&S_NoiseNetwork","MME_NoiseNetwork"),name="Estimators",values=c("black", "grey35","grey70"))+ 
  geom_errorbar(aes(ymin=low, ymax=upper,color=factor(type),group=factor(type)), width=.001) 
 




