library(igraph)
library(Matrix)
library(dplyr)
library(ggplot2)

################ functions ##################
# get observed networks
getObs<-function(data.true,Nv,alpha,beta)
{
  
  data.obs <- rsparsematrix(Nv, Nv, nnz = rbinom(1,Nv^2,alpha), rand.x = function(n) return(1))
  ind1 <- which(data.true==1)
  data.obs[ind1] <- rbinom(length(ind1),1,1-beta)
  diag(data.obs)<- 0
  data.obs <-  Matrix::forceSymmetric(data.obs,uplo="L")
  return(data.obs)
}

# get the value of branching factor in noisy networks
getRatio <- function(alpha,beta,Nv,data.true,m) 
{
  ratio.obs <- rep(NA,m) 
  for(j in 1:m)
  {
    data.obs <- getObs(data.true,Nv,alpha,beta)
    degree.obs <- colSums(data.obs)
    ratio.obs[j] <- sum(degree.obs^2)/sum(degree.obs) 
  }
  return(ratio.obs)
}

# get bootstrap CI
getCI<-function(boot_case,boot_index,boot_n,O.c)
{
  a <- sapply(1:boot_n,function(i) apply(boot_case[boot_index[i,],],2,mean))
  b <- sapply(1:boot_n,function(i) apply(boot_case[boot_index[i,],],2,var))
  n_col <- ncol(boot_case)
  res <- matrix(c(apply(a,1,mean)-rep(O.c,n_col),apply(a,1,quantile,probs=0.025)-rep(O.c,n_col),
                  apply(a,1,quantile,probs=0.975)-rep(O.c,n_col),apply(b,1,mean),
                  apply(b,1,quantile,probs=0.025),apply(b,1,quantile,probs=0.975)),nrow=n_col)
  return(res)
}

################ load data ##################
indir <-  "../../Data/RandomNetworks/"
load( paste0(indir,"randomNetwork.RData") )

################ simulation ##################
# number of trials
m <- 1e4
# number of nodes
Nv <- 1e4

# case 1: ER networks, average degree = 50, beta = 0.1, 0.2, 0.3
index <- type1
data.true <- sparseMatrix(i=index[,1],j=index[,2],x = 1, symmetric = T)
NumOfEdge <- nrow(index)
NumOfNoEdge <- choose(Nv,2)-NumOfEdge
beta <- c(0.1,0.2,0.3)
alpha <- beta*NumOfEdge/NumOfNoEdge
case1 <- sapply(1:length(beta),function(i) getRatio(alpha[i],beta[i],Nv,data.true,m)  )

# case 2: ER networks, average degree = 100, beta = 0.1, 0.2, 0.3
index <- type2
data.true <- sparseMatrix(i=index[,1],j=index[,2],x = 1, symmetric = T)
NumOfEdge <- nrow(index)
NumOfNoEdge <- choose(Nv,2)-NumOfEdge
alpha <- beta*NumOfEdge/NumOfNoEdge
case2 <- sapply(1:length(beta),function(i) getRatio(alpha[i],beta[i],Nv,data.true,m)  )
 
# case 3: scale-free networks, average degree = 50, beta = 0.1, 0.2, 0.3
index <- type3
data.true <- sparseMatrix(i=index[,1],j=index[,2],x = 1, symmetric = T)
NumOfEdge <- nrow(index)
NumOfNoEdge <- choose(Nv,2)-NumOfEdge
alpha <- beta*NumOfEdge/NumOfNoEdge
case3 <- sapply(1:length(beta),function(i) getRatio(alpha[i],beta[i],Nv,data.true,m)  )

# case 4: scale-free networks, average degree = 100, beta = 0.1, 0.2, 0.3
index <- type4
data.true <- sparseMatrix(i=index[,1],j=index[,2],x = 1, symmetric = T)
NumOfEdge <- nrow(index)
NumOfNoEdge <- choose(Nv,2)-NumOfEdge
alpha <- beta*NumOfEdge/NumOfNoEdge
case4 <- sapply(1:length(beta),function(i) getRatio(alpha[i],beta[i],Nv,data.true,m)  )

################ bootstrap CI ##################
case_all <- cbind(case1,case2,case3,case4)
degree_true <- degree_matrix
ratio_true <- colSums(degree_true^2)/colSums(degree_true)
boot_n <- 10^3
boot_index <- replicate(boot_n,sample(1:m,m,replace = T))
boot_index <- matrix(boot_index,nrow=boot_n,byrow = T)
res1 <- t(getCI(case_all[,1:3],boot_index,boot_n,ratio_true[1]))
res2 <- t(getCI(case_all[,4:6],boot_index,boot_n,ratio_true[2]))
res3 <- t(getCI(case_all[,7:9],boot_index,boot_n,ratio_true[3]))
res4 <- t(getCI(case_all[,10:12],boot_index,boot_n,ratio_true[4]))

################ Figure 3.1 ##################
data.long <- data.frame(beta=beta,degree=50,Low=res1[2,],means=res1[1,],
                        High=res1[3,],type="Homogeneous",measures="Bias") %>%
  add_row(beta=beta,degree=50,Low=res1[5,],means=res1[4,],
          High=res1[6,],type="Homogeneous",measures="Variance")%>%
  add_row(beta=beta,degree=100,Low=res2[2,],means=res2[1,],
          High=res2[3,],type="Homogeneous",measures="Bias")%>%
  add_row(beta=beta,degree=100,Low=res2[5,],means=res2[4,],
          High=res2[6,],type="Homogeneous",measures="Variance")%>%
  add_row(beta=beta,degree=50,Low=res3[2,],means=res3[1,],
          High=res3[3,],type="Inhomogeneous",measures="Bias") %>%
  add_row(beta=beta,degree=50,Low=res3[5,],means=res3[4,],
          High=res3[6,],type="Inhomogeneous",measures="Variance")%>%
  add_row(beta=beta,degree=100,Low=res4[2,],means=res4[1,],
          High=res4[3,],type="Inhomogeneous",measures="Bias") %>%
  add_row(beta=beta,degree=100,Low=res4[5,],means=res4[4,],
          High=res4[6,],type="Inhomogeneous",measures="Variance")%>% 
  mutate(degree=case_when(degree==50~ "average degree: 50",degree==100~ "average degree: 100"))

Fig3.1 <- data.long %>% ggplot(aes(x = beta, y = means,color=factor(type)))+
  geom_point(alpha = 0.7)+geom_errorbar(aes(ymin=Low, ymax=High), width=.005)+
  facet_grid(measures~degree, scales = "free")+ 
  theme_bw() +scale_x_continuous(breaks = c(0.1, .2, .3), name=expression(beta))+ ylab("")+
  theme(plot.title = element_text(hjust = 0.5,size=16),axis.line = element_line(colour = "black"),panel.grid.minor = element_blank(),panel.background = element_blank())+
  scale_color_discrete(name="Network",labels =c("Homogeneous","Inhomogeneous")  )






 
