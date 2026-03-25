# Make config file for missing data Ricker runs B

# Should include (1) data file (2) parameter file (3) number of nodes for cluster (4) save file location
# (5) beginnning index (6) ending index (7) optional- model list

##########################################
# Simulated data, random missingness
##########################################
# Split it into chunks- every 2500 
# Ricker A
n=75000/2500
# data file name
datFile=rep("data/missingDatasets/pois_sim_randMiss_trim_A.rds",n)
# par file name
parFile=rep("data/missingDatasets/pois_sim_params.rds",n)
# cluster size
clsize=rep(5,n)
# save file location
saveFile=paste0(rep("Model_Runs/RickerA_randMiss_resultTable",n),1:n,".rds")
# beginning index
index1=seq(1,75000,by=2500)
# ending index
index2=seq(2500,75000,by=2500)
configx=as.data.frame(matrix(data=c(1:n,datFile,parFile,clsize,saveFile,index1,index2),nrow=n,ncol=7,byrow=F))
colnames(configx)=c("ArrayTaskID","datFile","parFile","clsize","saveFile","index1","index2")
rownames(configx)=NULL
configx$ArrayTaskID<-as.numeric(configx$ArrayTaskID)
configx$clsize<-as.numeric(configx$clsize)
configx$index1<-as.numeric(configx$index1)
configx$index2<-as.numeric(configx$index2)

write.table(configx, file="/Users/amypatterson/Documents/Laramie_postdoc/Missing_data_TS/missing-data/Model_Runs/poiss_Sim_mods/RickerConfigA_randMiss.txt",sep=" ",row.names = F, quote=F)


# Ricker B
n=75000/2500
# data file name
datFile=rep("data/missingDatasets/pois_sim_randMiss_trim_B.rds",n)
# par file name
parFile=rep("data/missingDatasets/pois_sim_params.rds",n)
# cluster size
clsize=rep(5,n)
# save file location
saveFile=paste0(rep("Model_Runs/RickerB_randMiss_resultTable",n),1:n,".rds")
# beginning index
index1=seq(1,75000,by=2500)
# ending index
index2=seq(2500,75000,by=2500)
configx=as.data.frame(matrix(data=c(1:n,datFile,parFile,clsize,saveFile,index1,index2),nrow=n,ncol=7,byrow=F))
colnames(configx)=c("ArrayTaskID","datFile","parFile","clsize","saveFile","index1","index2")
rownames(configx)=NULL
configx$ArrayTaskID<-as.numeric(configx$ArrayTaskID)
configx$clsize<-as.numeric(configx$clsize)
configx$index1<-as.numeric(configx$index1)
configx$index2<-as.numeric(configx$index2)

write.table(configx, file="/Users/amypatterson/Documents/Laramie_postdoc/Missing_data_TS/missing-data/Model_Runs/poiss_Sim_mods/RickerConfigB_randMiss.txt",sep=" ",row.names = F, quote=F)

##########################################
# Simulated data, MinMax missingness
##########################################
# Split it into chunks- every 2500 
# Ricker 
n=75000/2500
# data file name
datFile=rep("data/missingDatasets/pois_sim_MinMaxMiss_trim.rds",n)
# par file name
parFile=rep("data/missingDatasets/pois_sim_params.rds",n)
# cluster size
clsize=rep(5,n)
# save file location
saveFile=paste0(rep("Model_Runs/Ricker_MinMaxMiss_resultTable",n),1:n,".rds")
# beginning index
index1=seq(1,75000,by=2500)
# ending index
index2=seq(2500,75000,by=2500)
configx=as.data.frame(matrix(data=c(1:n,datFile,parFile,clsize,saveFile,index1,index2),nrow=n,ncol=7,byrow=F))
colnames(configx)=c("ArrayTaskID","datFile","parFile","clsize","saveFile","index1","index2")
rownames(configx)=NULL
configx$ArrayTaskID<-as.numeric(configx$ArrayTaskID)
configx$clsize<-as.numeric(configx$clsize)
configx$index1<-as.numeric(configx$index1)
configx$index2<-as.numeric(configx$index2)

write.table(configx, file="/Users/amypatterson/Documents/Laramie_postdoc/Missing_data_TS/missing-data/Model_Runs/poiss_Sim_mods/RickerConfig_MinMaxMiss.txt",sep=" ",row.names = F, quote=F)



##########################################

# Other config files:
# Make config file for missing data Ricker runs- we need to make each job have just 1 data set
# so that we can rerun the data set if the fatal cholesky error happens

# Make config file for missing data Ricker runs
# Should include (1) data file (2) parameter file (3) number of nodes for cluster (4) save file location
# (5) beginnning index (6) ending index (7) optional- model list (8) seed- so that we can change see for reruns
# Split it into 10 arrays of 7500 each
# Ricker A MI
for(i in 1:10){
  n=7500
  # data file name
  datFile=rep("data/missingDatasets/pois_sim_randMiss_A.rds",n)
  # par file name
  parFile=rep("data/missingDatasets/pois_sim_params.rds",n)
  # cluster size
  clsize=rep(5,n)
  # save file location
  saveFile=paste0(rep("Model_Runs/poiss_Sim_mods/RickerA_resultTable_",n),1:n,".csv")
  # beginning index
  index1=seq((i-1)*7500+1,i*7500,by=1)
  # ending index
  index2=index1
  seed=rep(1493,n)
  configx=as.data.frame(matrix(data=c(1:n,datFile,parFile,clsize,saveFile,index1,index2,seed),nrow=n,ncol=8,byrow=F))
  colnames(configx)=c("ArrayTaskID","datFile","parFile","clsize","saveFile","index1","index2","seed")
  rownames(configx)=NULL
  configx$ArrayTaskID<-as.numeric(configx$ArrayTaskID)
  configx$clsize<-as.numeric(configx$clsize)
  configx$index1<-as.numeric(configx$index1)
  configx$index2<-as.numeric(configx$index2)
  
  write.table(configx, file=paste0("/Users/amypatterson/Documents/Laramie_postdoc/Missing_data_TS/missing-data/Model_Runs/poiss_Sim_mods/RickerConfigA_",i,".txt"),sep=" ",row.names = F, quote=F)
  
}




# Ricker B MI
for(i in 1:10){
  n=7500
  # data file name
  datFile=rep("data/missingDatasets/pois_sim_randMiss_B.rds",n)
  # par file name
  parFile=rep("data/missingDatasets/pois_sim_params.rds",n)
  # cluster size
  clsize=rep(5,n)
  # save file location
  saveFile=paste0(rep("Model_Runs/poiss_Sim_mods/RickerB_resultTable_",n),1:n,".csv")
  # beginning index
  index1=seq((i-1)*7500+1,i*7500,by=1)
  # ending index
  index2=index1
  seed=rep(1493,n)
  configx=as.data.frame(matrix(data=c(1:n,datFile,parFile,clsize,saveFile,index1,index2,seed),nrow=n,ncol=8,byrow=F))
  colnames(configx)=c("ArrayTaskID","datFile","parFile","clsize","saveFile","index1","index2","seed")
  rownames(configx)=NULL
  configx$ArrayTaskID<-as.numeric(configx$ArrayTaskID)
  configx$clsize<-as.numeric(configx$clsize)
  configx$index1<-as.numeric(configx$index1)
  configx$index2<-as.numeric(configx$index2)
  
  write.table(configx, file=paste0("/Users/amypatterson/Documents/Laramie_postdoc/Missing_data_TS/missing-data/Model_Runs/poiss_Sim_mods/RickerConfigB_",i,".txt"),sep=" ",row.names = F, quote=F)
  
}

# Ricker MinMaxMiss MI
for(i in 1:2){
  n=8000
  # data file name
  datFile=rep("data/missingDatasets/pois_sim_minMaxMiss.rds",n)
  # par file name
  parFile=rep("data/missingDatasets/pois_sim_params.rds",n)
  # cluster size
  clsize=rep(5,n)
  # save file location
  saveFile=paste0(rep("Model_Runs/poiss_Sim_mods/RickerMinMaxMiss_MI_",n),1:n,".csv")
  # beginning index
  index1=seq((i-1)*n+1,i*n,by=1)
  # ending index
  index2=index1
  seed=rep(1493,n)
  configx=as.data.frame(matrix(data=c(1:n,datFile,parFile,clsize,saveFile,index1,index2,seed),nrow=n,ncol=8,byrow=F))
  colnames(configx)=c("ArrayTaskID","datFile","parFile","clsize","saveFile","index1","index2","seed")
  rownames(configx)=NULL
  configx$ArrayTaskID<-as.numeric(configx$ArrayTaskID)
  configx$clsize<-as.numeric(configx$clsize)
  configx$index1<-as.numeric(configx$index1)
  configx$index2<-as.numeric(configx$index2)
  
  write.table(configx, file=paste0("/Users/amypatterson/Documents/Laramie_postdoc/Missing_data_TS/missing-data/Model_Runs/poiss_Sim_mods/RickerConfigMinMax_MI_",i,".txt"),sep=" ",row.names = F, quote=F)
  
}
