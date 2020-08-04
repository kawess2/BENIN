#### This code is for combining BENIN output

# Set the working directory
home=="/Volumes/Seagate\ Backup\ Plus\ Drive/PhD/GRN_inference_prime/src"
setwd(home)
##########
### 
source("R/utile.R")
source("R/benin.R")
source("R/bootpracs.q")
source("R/bootfuns.q")
#### function for merging

network=5
sizenetwork=10
#### getting all file in the specific file
#locData
prior="locData"
pathfile=paste0("../final_res_dream4/",sizenetwork,"_",network,"/",prior)
fileforlocD<-list.files(pathfile,all.files=F, full.names=T,no.. = T)
#KO
prior="KO"
pathfile=paste0("../final_res_dream4/",sizenetwork,"_",network,"/",prior)
fileforKO<-list.files(pathfile,all.files=F, full.names=T,no.. = T)

result_KO_data<-read.table(fileforKO, header=F, stringsAsFactors = F, sep="\t")
colnames(result_KO_data)<-c("TF","TG","W")
method="average"
#nameoutputfile=paste0("../final_res_dream4/",sizenetwork,"_",network,"/",method,"globalres_benin_combined_",method,
  #                    "_",sizenetwork,"_",network,".txt")
for (namefile in fileforlocD)
{
 
  if (sizenetwork==10)
  {
    result_loc_data<-read.table(namefile, header=T, stringsAsFactors = F, sep="\t",row.names = 1)
    result_loc_data<-subset(result_loc_data,TF!=TG)
  }
  else
  {
    result_loc_data<-read.table(namefile, header=F, stringsAsFactors = F, sep="\t")
    colnames(result_loc_data)<-c("TF","TG","W")
  }
  merged_res<-merge(result_KO_data,result_loc_data,by=c("TF","TG"), all=T,suffixes = c(".KO",".locD"))
  #print(file)
  if (method=="max")
  {
    merged_res$W<-apply(merged_res[,3:4],1,FUN=function(x){max(x)})
  }
  else{
    if (method=="average")
    {
      merged_res$W<-apply(merged_res[,3:4],1,FUN=function(x){mean(x)})
    }
  }

  
  merged_res<-merged_res[,c("TF","TG","W")]
  merged_res<-merged_res[order(merged_res$W,decreasing = T),]
  if(sizenetwork==10)
  {
    file<-strsplit(namefile, "_")[[1]][6]
  }
  else
  {
    if(network==3)
    {
      file<-strsplit(namefile, "_")[[1]][10]
    }
    else
    {
      file<-strsplit(namefile, "_")[[1]][12]
    }
  }
  
  
  
  nameoutputfile=paste0("../final_res_dream4/",sizenetwork,"_",network,"/combined/",method,"/globalres_combined_benin_KO_locData_",method,
                        "_",file,"_",sizenetwork,"_",network,".txt")
  print(nameoutputfile)
  savedata(merged_res,nameoutputfile, colnames=F, rownames=F,sep="\t",quote=FALSE,append=FALSE)
}
