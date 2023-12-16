# Make config file for missing data Ricker runs

# Should include (1) data file (2) parameter file (3) number of nodes for cluster (4) save file location
# (5) beginnning index (6) ending index (7) optional- model list

# Split it into chunks- every 2500 
#n=75000/2500
task1=21
n1start=50001
n1finish=52500

task2=26
n2start=62501
n2finish=65000

task3=27
n3start=65001
n3finish=67500

task4=28
n4start=67501
n4finish=70000

n=100

# data file name
datFile=rep("data/missingDatasets/pois_sim_randMiss_B.rds",n)

# par file name
parFile=rep("data/missingDatasets/pois_sim_params.rds",n)

# cluster size
clsize=rep(5,n)

# save file location
numbers1=paste0(task1,"_",1:25)
numbers2=paste0(task2,"_",1:25)
numbers3=paste0(task3,"_",1:25)
numbers4=paste0(task4,"_",1:25)
numbersAll=c(numbers1,numbers2,numbers3,numbers4)
saveFile=paste0(rep("Model_Runs/RickerB_resultTable",n),numbersAll,".rds")

# beginning index
beginInd1=seq(n1start,n1finish,by=100)
beginInd2=seq(n2start,n2finish,by=100)
beginInd3=seq(n3start,n3finish,by=100)
beginInd4=seq(n4start,n4finish,by=100)
index1=c(beginInd1,beginInd2,beginInd3,beginInd4)
#index1=seq(1,75000,by=2500)

# ending index
endInd1=seq(n1start+99,n1finish,by=100)
endInd2=seq(n2start+99,n2finish,by=100)
endInd3=seq(n3start+99,n3finish,by=100)
endInd4=seq(n4start+99,n4finish,by=100)
index2=c(endInd1,endInd2,endInd3,endInd4)
#index2=seq(2500,75000,by=2500)

configx=as.data.frame(matrix(data=c(1:n,datFile,parFile,clsize,saveFile,index1,index2),nrow=n,ncol=7,byrow=F))

colnames(configx)=c("ArrayTaskID","datFile","parFile","clsize","saveFile","index1","index2")
rownames(configx)=NULL

write.table(configx, file="/Users/amypatterson/Documents/Laramie_postdoc/Missing_data_TS/missing-data/Model_Runs/RickerConfig.txt",sep=" ",quote=F,row.names = F)

