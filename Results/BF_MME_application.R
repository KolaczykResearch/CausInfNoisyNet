library(Matrix)
library(xtable)
library(igraph)
library(ggplot2)
library(dplyr)
library("RColorBrewer")
library(wesanderson)

############ functions ############
# get the value of branching factor
getBF <- function(degree)
{
  mean(degree^2)/mean(degree)
}
# get the value of R0
getR0 <- function(bf,theta,gamma)
{
  theta/(theta+gamma)*(bf-1)
}
source("../../R/MME_BF.R")

############ hospital ############
# load data
indir <-  "../../Data/HospitalNetworks/"
load( paste0(indir,"Replicates.RData") )
 
# get MME of Branching factor
P <- adj_list[[1]]@Dim[1]
alpha0 <- 0.1
Y <- as.matrix(adj_list[[1]])
Ystar <-  as.matrix(adj_list[[2]])
Ystar2 <- as.matrix(adj_list[[3]])
estBF_123 <- MME_BF(alpha0,P, Ystar2, Ystar, Y) ## Table 5.2 (hospital)

BF_est<- unlist(estBF_123[c(6,12)])[c(2,1,3)] +1
BF_obs <- sapply(1:3, function(i) getBF(colSums(adj_list[[i]])))
BF_hosp <- data.frame(Low=BF_est[1],means=BF_est[2],High=BF_est[3],type="MME") %>%
  add_row(means=BF_obs,type=paste0("Observed ",1:3)) %>%  mutate(Data="Hosptial",type2=case_when(type=="MME"~1, type!="MME"~2) )
table_hosp <- c(estBF_123$Alpha,estBF_123$cfAlpha,estBF_123$Beta,estBF_123$cfBeta) # MME of alpha and beta

# get MME of R0
theta_seq=c(0.016,0.026)
gamma_seq = seq(1/8, 1/24.6, length.out = 4)
R0_hosp <- expand.grid(theta=theta_seq,gamma=gamma_seq) %>% mutate(mme=theta/(theta+gamma)*(BF_est[2]-1),tmp=(BF_est[2]-BF_est[1])*theta/(theta+gamma),
                                                                  Low=mme-tmp,High=mme+tmp,contact=1,Data="Hospital") %>%
  select(-tmp)
 
 

 
############ school ############
# load data
indir <-  "../../Data/SchoolNetworks/"
load( paste0(indir,"Replicates.RData") )

# get MME of Branching factor
P1_undiCol <- school1_undiCol_rep[[1]]@Dim[1]
P2_undiCol <- school2_undiCol_rep[[1]]@Dim[1]
P3_undiCol <- school3_undiCol_rep[[1]]@Dim[1]
P4_undiCol <- school4_undiCol_rep[[1]]@Dim[1]

school1_undiCol_BF123 <- MME_BF(alpha0,P1_undiCol, as.matrix(school1_undiCol_rep[[1]]), as.matrix(school1_undiCol_rep[[2]]), as.matrix(school1_undiCol_rep[[3]]))
school2_undiCol_BF123 <- MME_BF(alpha0,P2_undiCol, as.matrix(school2_undiCol_rep[[1]]), as.matrix(school2_undiCol_rep[[2]]), as.matrix(school2_undiCol_rep[[3]]))
school3_undiCol_BF123 <- MME_BF(alpha0,P3_undiCol, as.matrix(school3_undiCol_rep[[2]]), as.matrix(school3_undiCol_rep[[3]]), as.matrix(school3_undiCol_rep[[1]]))
school4_undiCol_BF123 <- MME_BF(alpha0,P4_undiCol, as.matrix(school4_undiCol_rep[[2]]), as.matrix(school4_undiCol_rep[[3]]), as.matrix(school4_undiCol_rep[[1]])) 

school1_undiCol_BFseq <- sapply(1:4,function(i) getBF(colSums(school1_undiCol_rep[[i]])) )
school2_undiCol_BFseq <- sapply(1:4,function(i) getBF(colSums(school2_undiCol_rep[[i]])) )
school3_undiCol_BFseq <- sapply(1:4,function(i) getBF(colSums(school3_undiCol_rep[[i]])) )
school4_undiCol_BFseq <- sapply(1:4,function(i) getBF(colSums(school4_undiCol_rep[[i]])) )


BF_undiCol_est<- matrix(unlist(c(school1_undiCol_BF123[c(6,12)],school2_undiCol_BF123[c(6,12)],school3_undiCol_BF123[c(6,12)],school4_undiCol_BF123[c(6,12)])),ncol = 4)[c(2,1,3),] +1

BF_school <- data.frame(Low=BF_undiCol_est[1,],means=BF_undiCol_est[2,],High=BF_undiCol_est[3,],
                        Data=paste0("School ",1:4),type="MME") %>% add_row(means=school1_undiCol_BFseq,
                                                                           Data="School 1",type=paste0("Observed ",1:4))%>% 
  add_row(means=school2_undiCol_BFseq, Data="School 2",type=paste0("Observed ",1:4)) %>%
  add_row(means=school3_undiCol_BFseq, Data="School 3",type=paste0("Observed ",1:4)) %>%
  add_row(means=school4_undiCol_BFseq, Data="School 4",type=paste0("Observed ",1:4)) %>%
  mutate(type2=case_when(type=="MME"~1, type!="MME"~2) )

table_school <- matrix(c(school1_undiCol_BF123$Alpha,school1_undiCol_BF123$cfAlpha,school1_undiCol_BF123$Beta,school1_undiCol_BF123$cfBeta,
                         school2_undiCol_BF123$Alpha,school2_undiCol_BF123$cfAlpha,school2_undiCol_BF123$Beta,school2_undiCol_BF123$cfBeta,
                         school3_undiCol_BF123$Alpha,school3_undiCol_BF123$cfAlpha,school3_undiCol_BF123$Beta,school3_undiCol_BF123$cfBeta,
                         school4_undiCol_BF123$Alpha,school4_undiCol_BF123$cfAlpha,school4_undiCol_BF123$Beta,school4_undiCol_BF123$cfBeta),byrow = TRUE,nrow = 4 )  # MME of alpha and beta

# get MME of R0
R0_school <- expand.grid(theta=theta_seq,gamma=gamma_seq) %>%
  merge(data.frame(bf_low=BF_undiCol_est[1,],bf=BF_undiCol_est[2,],type=seq(1,4)))%>% 
  mutate(mme=theta/(theta+gamma)*(bf-1),tmp=(bf-bf_low)*theta/(theta+gamma),
         Low=mme-tmp,High=mme+tmp,contact=2,Data=case_when(type==1 ~ "School 1",type==2 ~ "School 2",
                                                           type==3 ~ "School 3",type==4 ~ "School 4"))%>%
  select(-c(tmp,bf,bf_low,type)) 

############ Table 5.1 ############
table_all <- rbind(table_hosp,table_school)


############ Figure 5.2 ############
BF_res <- rbind(BF_hosp,BF_school)
Fig5.2 <- BF_res %>%   mutate(Data = as.numeric(factor(Data) ) + as.numeric(factor(type2))/10 - 2/10 ) %>%
  ggplot(aes(x = Data, y = means,color=factor(type),shape=factor(type2)))+geom_point()+ theme_bw() +
  geom_errorbar(aes(ymin=Low, ymax=High), width=.05, show.legend=FALSE)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5),
        panel.grid.minor = element_blank(),panel.background = element_blank())+
  guides(colour = guide_legend(override.aes = list(shape = c(16,rep(17,4) ) )))+
  scale_color_manual(values = c("Black",wes_palette(n=4, name="GrandBudapest2") ) )+
  labs(color  = "Type")+scale_shape(guide = FALSE) +xlab(" ")+ylab(" ")+geom_blank(aes(y = 0))  +
  scale_x_continuous(breaks = 1:5,labels =c("Hospital",paste0("School ",1:4))  )

 
############ Figure 5.3 ############
R0_res <- rbind(R0_hosp,R0_school)
Fig5.3 <- R0_res %>% mutate(gamma = as.numeric(as.character(gamma)) + as.numeric(factor(Data))/400 - .0075 ) %>%
  ggplot(aes(x = gamma, y = mme,color=factor(Data),shape=factor(contact)))+geom_point(alpha = 0.7)+geom_errorbar(aes(ymin=Low, ymax=High), width=.002)+
  facet_grid(cols = vars(theta),labeller = label_both)+scale_x_continuous(breaks = round(gamma_seq,3),name=expression(gamma))+
  ylab(" ") +   theme_bw() +theme(panel.border = element_rect(fill=NA,color="black", size=0.5),
                                  panel.grid.minor = element_blank(),panel.background = element_blank())+
  geom_hline(yintercept=1, linetype="dashed")+guides(colour = guide_legend(override.aes = list(shape = c(16,rep(17,4) ) )))+
  labs(color  = "Network")+scale_shape(guide = FALSE)+geom_blank(aes(y = 0))  




