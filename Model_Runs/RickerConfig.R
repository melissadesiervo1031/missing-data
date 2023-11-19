# Make config file for missing data Ricker runs

# Should include (1) data file (2) parameter file (3) number of nodes for cluster (4) save file location
# (5) beginnning index (6) ending index (7) optional- model list

# Split it into chunks- every 2500 
n=75000/2500

# data file name
datFile=rep("Model_Runs/modelruns_ricker.R",n)

# par file name
parFile=rep("data/missingDatasets/pois_sim_randMiss_A.rds",n)

# cluster size
clsize=rep(8,n)

# save file location
saveFile=paste0(rep("Model_Runs/RickerA_resultTable",n),1:n,".rds")

# beginning index
index1=seq(1,75000,by=2500)

# ending index
index2=seq(2500,75000,by=2500)

configx=as.data.frame(matrix(data=c(1:n,datFile,parFile,clsize,saveFile,index1,index2),nrow=n,ncol=7,byrow=F))

colnames(configx)=c("ArrayTaskID","datFile","parFile","clsize","saveFile","index1","index2")
rownames(configx)=NULL

write.table(configx, file="/Users/amypatterson/Documents/Laramie_postdoc/Missing_data_TS/missing-data/Model_Runs/RickerConfig.txt",sep=" ",row.names = F)

