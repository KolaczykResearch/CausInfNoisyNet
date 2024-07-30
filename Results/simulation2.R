library(igraph)
library(Matrix)
library(ggplot2)
library(dplyr)

################ functions ##################
 
# get estimators for a Monte Carlo simulation trial
getOne <- function(data.network,alpha,beta,p,O.c,z,alpha0, Nb, tiny,Zcv=1.96,Nb2=10000)
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
    VarEst <- getVarEst(data.true,data.obs1,data.obs2,data.obs3,p,Nv,alphaEst,betaEst,z,O.c,result,Nb2,Zcv)
    return(c(result,VarEst))
}
 
# transform data to desired format for plotting
getRes <- function(school,alpha,beta,res_sim,O.c,Zcv=1.96)
{
  est_bias <- colMeans(res_sim)[1:12]-rep(O.c,3)
  sim_sd <- apply(res_sim,2,sd)[1:12]
  est_se <- colMeans(res_sim)[13:24]
  mae_se <- rowMeans(abs(sim_sd-t(res_sim[,13:24])))
  est_Coverage <- rowMeans(sapply(1:nrow(res_sim), function(i) rep(O.c,3) < (res_sim[i,1:12]+Zcv*res_sim[i,13:24]) & rep(O.c,3) > (res_sim[i,1:12]-Zcv*res_sim[i,13:24])))
  res <- cbind(rep(school,12),rep(alpha,12),rep(beta,12),rep(c(1,2,3),each=4),rep(1:4,3),est_bias,sim_sd,mae_se,est_Coverage)
  return(res)
}
getRes_inhom <- function(school,alpha,beta,res_sim,O.c,Zcv=1.96)
{
  est_bias <- colMeans(res_sim)[1:12]-rep(O.c,3)
  sim_sd <- apply(res_sim,2,sd)[1:12]
  est_se <- colMeans(res_sim)[13:24]
  mae_se <- rowMeans(abs(sim_sd-t(res_sim[,13:24])))
  est_Coverage <- rowMeans(sapply(1:nrow(res_sim), function(i) rep(O.c,3) < (res_sim[i,1:12]+Zcv*res_sim[i,13:24]) & rep(O.c,3) > (res_sim[i,1:12]-Zcv*res_sim[i,13:24])))
  est_Coverage2 <- colMeans(res_sim)[25:32]
  res <- cbind(rep(school,12),rep(alpha,12),rep(beta,12),rep(c(1,2,3),each=4),rep(1:4,3),est_bias,sim_sd,mae_se,c(est_Coverage[1:4],est_Coverage2))
  return(res)
}


source("../R/getNoisyNet.R")
source("../R/getNetNoise.R")
source("../R/getEst.R")
source("../R/getVarEst.R")
################ load data ##################
indir <-  "../Data/"
load( paste0(indir,"dataTrue.RData") )
load( paste0(indir,"z.RData") )

################ simulation ##################
# number of Monte Carlo simulation trials
m <- 10^4

# probability of treatment
p <- .1

# values of potential outcomes
O.c<- c(10,7,5,1)

# parameters for algorithm 1 (estimate alpha and beta)
alpha0 <- .001
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

#  biases, standard deviations, mean absolute errors of standard errors, and coverage rates of 95% confidence intervals for 4 schools
school_res <- list()
school_res[[1]] <- Reduce("rbind",lapply(1:length(beta),function (i) getRes(1,alpha[i],beta[i],res_est1[[i]],O.c)))
school_res[[2]] <- Reduce("rbind",lapply(1:length(beta),function (i) getRes(2,alpha[i],beta[i],res_est2[[i]],O.c)))
school_res[[3]] <- Reduce("rbind",lapply(1:length(beta),function (i) getRes(3,alpha[i],beta[i],res_est3[[i]],O.c)))
school_res[[4]] <- Reduce("rbind",lapply(1:length(beta),function (i) getRes_inhom (4,alpha[i],beta[i],res_est4[[i]],O.c)))
for (i in 1:n_school) {
  school_res[[i]] <- cbind(rep(1:6,each=12),school_res[[i]])
}
school_res <- Reduce("rbind",school_res)
colnames(school_res) <- c("case","school","alpha","beta","type","level","Bias","SD","MAE","Coverage")
school_res <- data.frame(school_res)%>% mutate(level=case_when(level==1 ~ "c11",level==2 ~ "c10",
                                                               level==3 ~ "c01",level==4 ~ "c00") ) %>%
  mutate(Bias=as.numeric(as.character(Bias)),SD=as.numeric(as.character(SD)),MAE=as.numeric(as.character(MAE)),Coverage=as.numeric(as.character(Coverage)))%>%
  rename(c('MAE(SE)'='MAE'))
res.long <- melt(school_res, id=c("case","school","alpha","beta","type","level"),measure=c("Bias","SD","MAE(SE)","Coverage") )

 
################ Figure 4.1 ##################
Fig4.1 <- res.long %>% mutate(case = as.numeric(as.character(case)) + as.numeric(school)/40 - .125*(type== 1)+ .125*(type== 3) +0.0125) %>%
ggplot(aes(x = case, y = value,shape = factor(school),color=factor(type),group=factor(type)))+
geom_point(alpha = 1)  +
facet_grid(variable~level,scales= "free") +  ylab(" ")  +
scale_x_continuous(name="Edge noise settings",breaks = 1:6,labels = paste0("case",1:6))+theme_bw(base_size = 14)+
theme(legend.position="bottom",legend.box="vertical", legend.margin=margin(),plot.title = element_text(hjust = 0.5,size=16),axis.line = element_line(colour = "black"),panel.grid.minor = element_blank(),panel.background = element_blank())+
scale_shape_manual(name="School",labels =c("School 1","School 2","School 3","School 4"),values=c(15,16,17,18))  +
scale_color_manual(labels=c("A&S_TrueNetwork","A&S_NoiseNetwork","MME_NoiseNetwork"),name="Estimators",values=c("black", "grey35","grey70")) +
guides(fill=guide_legend(nrow=2,byrow=TRUE))
