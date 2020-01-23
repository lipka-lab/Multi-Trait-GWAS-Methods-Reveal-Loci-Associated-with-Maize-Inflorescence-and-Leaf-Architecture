myGM<-read.table("GAPIT.Genotype.map55KLinkImputeFiltered.txt",header = T)
myGD<-read.table("GAPIT.Genotype.Numerical55KLinkImputeFiltered.txt",header = T)
source("Function to detrimine false and true postive.R")
Settings<-c("GWAS_1_NullTrait","GWAS_2","GWAS_3","GWAS_4","GWAS_5","GWAS_6",
            "GWAS_2_partially","GWAS_3_partially","GWAS_4_partially",
            "GWAS_5_partially","GWAS_6_partially")
S<-Settings[1]
h1<-0
h2<-0
h3<-0
summary<-as.data.frame(matrix(NA,nrow = 1,ncol = 7))
setwd(paste("~/Box/brians_paper/SimulatedPheno2.0/",S,sep=""))
home.dir<-getwd()
if(S=="GWAS_1_NullTrait"){
  summary<-as.data.frame(matrix(NA,nrow = 1,ncol = 5))
  setwd(paste(home.dir,"/output",sep=""))
  for(j in 698:1000){
    print(paste("Rep",j))
        for (A in c(0.05,0.1)){
          for (thres in c("FDR", "bonferoni")){
            destfile=paste("Simulated_Data__Rep",j,"_Herit_",h1,"_",h2,"_",h3,".assoc.txt",sep="")
            if (file.exists(destfile)) {
  GWAS.RESULTS<-read.table(paste("Simulated_Data__Rep",j,"_Herit_",h1,"_",h2,"_",h3,".assoc.txt",sep=""),header=T)
  FDR<-p.adjust(GWAS.RESULTS$p_lrt, method = "fdr" , n = nrow(GWAS.RESULTS))
  GWAS.RESULTS<-cbind(GWAS.RESULTS,FDR)
  summary.temp<-data.frame(Trait="MV",Rep=j,Alpha=A,Threshold=thres,simulation.success(c,GWAS.RESULTS,
                                                                                                  myGM=myGM,
                                                                                                  myGD=myGD,
                                                                                                  chr=1,
                                                                                                  pos=3,
                                                                                                  rs=2,
                                                                                                  pvalue=17,
                                                                                                  FDR=18,
                                                                                                  Threshold=thres,#"bonferoni",  #options: "FDR" "bonferoni" "N" 
                                                                                                  Alpha=A,
                                                                                                  bp=TRUE,LD.logic=TRUE,
                                                                                                  LD.measurement="R2", #options: D  D' r  R2 X2
                                                                                                  removeCasual=TRUE,
                                                                                                  NULL.trait=TRUE))
            colnames(summary)<-colnames(summary.temp)
            summary<-rbind(summary,summary.temp)
            
            colnames(summary)<-colnames(summary.temp)
            #summary<-rbind(summary,summary.temp)
          }
       } }
    #put in info for getting false detection rate for null settings
  }
} else {
  setwd(paste(home.dir,"/output",sep=""))
  
  casualmutations<-fread("Additive_selected_QTNs.txt")
  casualmutations<-casualmutations[-c(which(casualmutations$type=="trait_specific")),]
  for(j in 1:100){
    print(paste("Rep",j))

    c<- casualmutations[which(casualmutations$rep==j)]
    c<-c[1:3,]
    destfile=paste("Simulated_Data__Rep",j,"_Herit_",h1,"_",h2,"_",h3,".assoc.txt",sep="")
    if (file.exists(destfile)) {
     GWAS.RESULTS<-read.table(paste("Simulated_Data__Rep",j,"_Herit_",h1,"_",h2,"_",h3,".assoc.txt",sep=""),header=T)
            FDR<-p.adjust(GWAS.RESULTS$p_lrt, method = "fdr" , n = nrow(GWAS.RESULTS))
            GWAS.RESULTS<-cbind(GWAS.RESULTS,FDR)
        
            for (A in c(0.05,0.1)){
          for (thres in c("FDR", "bonferoni")){
            
           
            summary.temp<-data.frame(Trait="MV",Rep=j,Alpha=A,Threshold=thres,simulation.success(c,GWAS.RESULTS,
                                                                                                  myGM=myGM,
                                                                                                  myGD=myGD,
                                                                                                  chr=1,
                                                                                                  pos=3,
                                                                                                  rs=2,
                                                                                                  pvalue=17,
                                                                                                  FDR=18,
                                                                                                  Threshold=thres,#"bonferoni",  #options: "FDR" "bonferoni" "N" 
                                                                                                  Alpha=A,
                                                                                                  bp=TRUE,LD.logic=TRUE,
                                                                                                  LD.measurement="R2", #options: D  D' r  R2 X2
                                                                                                  removeCasual=TRUE,
                                                                                                  NULL.trait=FALSE))
            colnames(summary)<-colnames(summary.temp)
            summary<-rbind(summary,summary.temp)
            
            colnames(summary)<-colnames(summary.temp)
            #summary<-rbind(summary,summary.temp)
          }
        } }
  }
} #end of else if if(S=="Setting_1_NullTrait"){
write.csv(summary[-1,],"~/Box/brians_paper/Results/Setting1MVSimulationResults.csv",row.names = F)#end of for(S in Settings){
