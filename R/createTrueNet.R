library(readxl)
library(reshape2)
library(igraph)

# set a destination for output
outdir <- "../../Data/"

# read data
data_all <- read_excel( paste0(outdir,"data.xlsx")) # available at https://doi.org/10.1371/journal.pone.0200090.s002
data <- data_all[which(data_all$participant!="0_ID_NA"),]
data[data=="0_ID_NA"] <- NA

school1 <- data[which(data$school==1),]
school1_1 <- school1[which(school1$round==1),(2:8)]
school1_2 <- school1[which(school1$round==2),(2:8)]
school1_3 <- school1[which(school1$round==3),(2:8)]
school1_4 <- school1[which(school1$round==4),(2:8)]

school2 <- data[which(data$school==2),]
school2_1 <- school2[which(school2$round==1),(2:8)]
school2_2 <- school2[which(school2$round==2),(2:8)]
school2_3 <- school2[which(school2$round==3),(2:8)]
school2_4 <- school2[which(school2$round==4),(2:8)]


school3 <- data[which(data$school==3),]
school3_1 <- school3[which(school3$round==1),(2:8)]
school3_2 <- school3[which(school3$round==2),(2:8)]
school3_3 <- school3[which(school3$round==3),(2:8)]
school3_4 <- school3[which(school3$round==4),(2:8)]

school4 <- data[which(data$school==4),]
school4_1 <- school4[which(school4$round==1),(2:8)]
school4_2 <- school4[which(school4$round==2),(2:8)]
school4_3 <- school4[which(school4$round==3),(2:8)]
school4_4 <- school4[which(school4$round==4),(2:8)]

# create observed networks for 4 school
school1_edge <-as.matrix( melt(school1[,(2:8)], id.vars = "participant",na.rm = TRUE)[,-2])
school1_di<-graph_from_edgelist(school1_edge, directed = TRUE)
school1_undiCol<- as.undirected(school1_di, mode = "collapse")

school1_1_edge <-as.matrix( melt(school1_1, id.vars = "participant",na.rm = TRUE)[,-2])
school1_1_di<-graph_from_edgelist(school1_1_edge, directed = TRUE)
school1_1_undiCol<- as.undirected(school1_1_di, mode = "collapse")

school1_2_edge <-as.matrix( melt(school1_2, id.vars = "participant",na.rm = TRUE)[,-2])
school1_2_di<-graph_from_edgelist(school1_2_edge, directed = TRUE)
school1_2_undiCol<- as.undirected(school1_2_di, mode = "collapse")

school1_3_edge <-as.matrix( melt(school1_3, id.vars = "participant",na.rm = TRUE)[,-2])
school1_3_di<-graph_from_edgelist(school1_3_edge, directed = TRUE)
school1_3_undiCol<- as.undirected(school1_3_di, mode = "collapse")

school1_4_edge <-as.matrix( melt(school1_4, id.vars = "participant",na.rm = TRUE)[,-2])
school1_4_di<-graph_from_edgelist(school1_4_edge, directed = TRUE)
school1_4_undiCol<- as.undirected(school1_4_di, mode = "collapse")

 
school2_edge <-as.matrix( melt(school2[,(2:8)], id.vars = "participant",na.rm = TRUE)[,-2])
school2_di<-graph_from_edgelist(school2_edge, directed = TRUE)
school2_undiCol<- as.undirected(school2_di, mode = "collapse")

school2_1_edge <-as.matrix( melt(school2_1, id.vars = "participant",na.rm = TRUE)[,-2])
school2_1_di<-graph_from_edgelist(school2_1_edge, directed = TRUE)
school2_1_undiCol<- as.undirected(school2_1_di, mode = "collapse")

school2_2_edge <-as.matrix( melt(school2_2, id.vars = "participant",na.rm = TRUE)[,-2])
school2_2_di<-graph_from_edgelist(school2_2_edge, directed = TRUE)
school2_2_undiCol<- as.undirected(school2_2_di, mode = "collapse")

school2_3_edge <-as.matrix( melt(school2_3, id.vars = "participant",na.rm = TRUE)[,-2])
school2_3_di<-graph_from_edgelist(school2_3_edge, directed = TRUE)
school2_3_undiCol<- as.undirected(school2_3_di, mode = "collapse")

school2_4_edge <-as.matrix( melt(school2_4, id.vars = "participant",na.rm = TRUE)[,-2])
school2_4_di<-graph_from_edgelist(school2_4_edge, directed = TRUE)
school2_4_undiCol<- as.undirected(school2_4_di, mode = "collapse")

school3_edge <-as.matrix( melt(school3[,(2:8)], id.vars = "participant",na.rm = TRUE)[,-2])
school3_di<-graph_from_edgelist(school3_edge, directed = TRUE)
school3_undiCol<- as.undirected(school3_di, mode = "collapse")

school3_1_edge <-as.matrix( melt(school3_1, id.vars = "participant",na.rm = TRUE)[,-2])
school3_1_di<-graph_from_edgelist(school3_1_edge, directed = TRUE)
school3_1_undiCol<- as.undirected(school3_1_di, mode = "collapse")

school3_2_edge <-as.matrix( melt(school3_2, id.vars = "participant",na.rm = TRUE)[,-2])
school3_2_di<-graph_from_edgelist(school3_2_edge, directed = TRUE)
school3_2_undiCol<- as.undirected(school3_2_di, mode = "collapse")

school3_3_edge <-as.matrix( melt(school3_3, id.vars = "participant",na.rm = TRUE)[,-2])
school3_3_di<-graph_from_edgelist(school3_3_edge, directed = TRUE)
school3_3_undiCol<- as.undirected(school3_3_di, mode = "collapse")

school3_4_edge <-as.matrix( melt(school3_4, id.vars = "participant",na.rm = TRUE)[,-2])
school3_4_di<-graph_from_edgelist(school3_4_edge, directed = TRUE)
school3_4_undiCol<- as.undirected(school3_4_di, mode = "collapse")

school4_edge <-as.matrix( melt(school4[,(2:8)], id.vars = "participant",na.rm = TRUE)[,-2])
school4_di<-graph_from_edgelist(school4_edge, directed = TRUE)
school4_undiCol<- as.undirected(school4_di, mode = "collapse")

school4_1_edge <-as.matrix( melt(school4_1, id.vars = "participant",na.rm = TRUE)[,-2])
school4_1_di<-graph_from_edgelist(school4_1_edge, directed = TRUE)
school4_1_undiCol<- as.undirected(school4_1_di, mode = "collapse")

school4_2_edge <-as.matrix( melt(school4_2, id.vars = "participant",na.rm = TRUE)[,-2])
school4_2_di<-graph_from_edgelist(school4_2_edge, directed = TRUE)
school4_2_undiCol<- as.undirected(school4_2_di, mode = "collapse")

school4_3_edge <-as.matrix( melt(school4_3, id.vars = "participant",na.rm = TRUE)[,-2])
school4_3_di<-graph_from_edgelist(school4_3_edge, directed = TRUE)
school4_3_undiCol<- as.undirected(school4_3_di, mode = "collapse")

school4_4_edge <-as.matrix( melt(school4_4, id.vars = "participant",na.rm = TRUE)[,-2])
school4_4_di<-graph_from_edgelist(school4_4_edge, directed = TRUE)
school4_4_undiCol<- as.undirected(school4_4_di, mode = "collapse")


# construct true networks for 4 schools
library(plyr)
df=data.frame(a=school1_edge[,1],b=school1_edge[,2])
school1_edge_rep <- ddply(df,.(a,b),nrow)
school1_edge_true <- as.matrix( school1_edge_rep[school1_edge_rep[,3]>1,-3])
school1_di_true<-graph_from_edgelist(school1_edge_true, directed = TRUE)
school1_undiCol_true<- as.undirected(school1_di_true, mode = "collapse")

df=data.frame(a=school2_edge[,1],b=school2_edge[,2])
school2_edge_rep <- ddply(df,.(a,b),nrow)
school2_edge_true <- as.matrix( school2_edge_rep[school2_edge_rep[,3]>1,-3])
school2_di_true<-graph_from_edgelist(school2_edge_true, directed = TRUE)
school2_undiCol_true<- as.undirected(school2_di_true, mode = "collapse")

df=data.frame(a=school3_edge[,1],b=school3_edge[,2])
school3_edge_rep <- ddply(df,.(a,b),nrow)
school3_edge_true <- as.matrix( school3_edge_rep[school3_edge_rep[,3]>1,-3])
school3_di_true<-graph_from_edgelist(school3_edge_true, directed = TRUE)
school3_undiCol_true<- as.undirected(school3_di_true, mode = "collapse")

df=data.frame(a=school4_edge[,1],b=school4_edge[,2])
school4_edge_rep <- ddply(df,.(a,b),nrow)
school4_edge_true <- as.matrix( school4_edge_rep[school4_edge_rep[,3]>1,-3])
school4_di_true<-graph_from_edgelist(school4_edge_true, directed = TRUE)
school4_undiCol_true<- as.undirected(school4_di_true, mode = "collapse")

save(school1_undiCol_true,school2_undiCol_true,school3_undiCol_true,school4_undiCol_true,file ="dataTrue.RData")

 
