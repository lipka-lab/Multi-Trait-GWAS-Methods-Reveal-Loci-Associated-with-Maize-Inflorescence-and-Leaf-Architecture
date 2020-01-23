pheno<-read.csv("~/Box/GermplasmMaterial/PhenoData/PanzeaMaizeArchitecturaltraits/BLUPS_282_GBSv4names.csv",header=T)
home.dir<-"~/Box/PhD Research/Univariate_and_PC_GWAS_in_282_GBSv4"
load("~/Box/GermplasmMaterial/GenomeData/GBSv4 282 GAPIT Input Files/GBSv4282kinshipPerChr.Rdata")
myCV<-read.csv("~/Box/GermplasmMaterial/GenomeData/GBSv4 282 GAPIT Input Files/covariatesforGAPIT.csv")
#replace NA with mean
myCV[277,7:8]<-NA
myCV$GDD_tassel<-as.numeric(as.character(myCV$GDD_tassel))
myCV$GDD_silk<-as.numeric(as.character(myCV$GDD_silk))
myCV[277,7]<-mean(na.omit(myCV$GDD_tassel))
myCV[277,8]<-mean(na.omit(myCV$GDD_silk))

load("~/Box/GermplasmMaterial/GenomeData/GBSv4 282 GAPIT Input Files/GBSv4282numericdata.Rdata")
myPC<-myCV[,1:6]

PC<-prcomp(na.omit(pheno[,c(4:ncol(pheno))]),center=T,scale=T)
PC.pheno<-data.frame(taxa<-na.omit(pheno[,c(1,4:ncol(pheno))])[1],PC$x)
pheno<-merge(pheno,PC.pheno,by=1,all=T)

#library(GAPIT3)
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
library('MASS')
source("http://www.zzlab.net/GAPIT/previous/gapit_functions20191108.txt") # to install GAPIT
source("http://zzlab.net/GAPIT/emma.txt") # to install themodified library
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

library(parallel)
#for(i in 1:10){
  LOCOGWAS<-function(chr){
    con <- file(paste("GAPITchr_",chr,".log",sep=""))
    sink(con, append=TRUE)
    sink(con, append=TRUE, type="message")
    
    setwd(home.dir)
  if(!file.exists("chr",chr,sep="")){
    dir.create(paste("chr",chr,sep="")) }
  myKI<-kinshiplist[[chr]]
  setwd(paste(home.dir,"/chr",chr,sep=""))
  
  myGAPIT <- GAPIT(
  Y = pheno[,c(1,4:ncol(pheno))],#[,c(1,3)],
  #then set it as wd
  file.path="~/Box/GermplasmMaterial/GenomeData/Maize282GBSvcf/Maize282GBSv4hmpfiles/",
  file.G ="LinkImpute_filtered_ZeaGBSv27_282_imputedV5_AGPv4_chr",
  file.from = chr,
  file.to = chr,
  #PCA.total=5,
  file.total = 1,
  file.Ext.G ="hmp.txt",
  CV=myCV, 
  K = myKI,
  SNP.fraction=1,
  group.from = 281,
  group.to = 281,
  Model.selection = T,
  SNP.MAF=0,
  SNP.FDR=0.05,
  model="MLM" ,
  PCA.View.output = F,
  Geno.View.output = F
)
  # Restore output to console
  setwd(home.dir)
  sink() 
  sink(type="message")
  }
mclapply(1:10,LOCOGWAS,mc.cores =detectCores()) 

setwd(home.dir)

#########################################################
#COMBINED CHR GWAS RESULT FILES
setwd(home.dir)
dir.create(paste("Complete GWAS Results",sep=""))
#read back in each chr gwas results for each traits
PC<-TRUE
for(t in 1:21){#colnames(pheno[4:ncol(pheno)])
  if(PC==FALSE) {
    n<-colnames(pheno)[t+3]
  } else { n<-paste("PC",t,sep="") }  
  print(n)
  for(i in 1:10){
    print(paste("reading in chr",i))
    setwd(paste(home.dir,"/Chr",i,sep=""))
    if(i==1){ 
      #combined GWAS result
      RESULTS<- read.csv(paste("GAPIT.MLM.",n,".GWAS.Results.csv",sep=""), head = TRUE)
    } else {
      GWAS.RESULTS<-read.csv(paste("GAPIT.MLM.",n,".GWAS.Results.csv",sep=""), head = TRUE)
      #GWAS.RESULTS$Chromosome<-i
      RESULTS<-rbind(GWAS.RESULTS,RESULTS)    
    } 
  }
  #remove FDR column
  RESULTS.new<-RESULTS[,-9]
  #calculate new FDR
  Conduct.FDR <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = RESULTS.new, 
                                                                    FDR.Rate = 0.05, FDR.Procedure = "BH")
  GWAS.Results.with.FDR<- Conduct.FDR$PWIP
  colnames(GWAS.Results.with.FDR)[9:10]<-c("effect","FRD Adjusted Pvalue")
  
  setwd(paste(home.dir,"/Complete GWAS Results",sep=""))
  
  #Make new Manhattan plots with the combined results
  GAPIT.Manhattan(GI.MP = GWAS.Results.with.FDR[,2:4], name.of.trait = n, 
                  DPP=50000, plot.type = "Genomewise",cutOff=0.00)
  #Mak a QQ plot of the results
  GAPIT.QQ(P.values = GWAS.Results.with.FDR$P.value, name.of.trait = n)
  #write new complete results
  write.table(GWAS.Results.with.FDR, paste("GAPIT.MLM.",n, ".GWAS.RESULTS.csv", sep = ""), 
              quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
  #print if any results passed threhold
  paste(n)
  print(length(which(RESULTS$FDR_Adjusted_P.values <=0.05)))
}

#prepping for METAL
PC<-TRUE
for(t in 1:21){#colnames(pheno[4:ncol(pheno)])
  if(PC==FALSE) {
    n<-colnames(pheno)[t+1]
  } else { n<-paste("PC",t,sep="") }  
  print(n)
  RESULTS<- fread(paste("GAPIT.MLM.",n,".GWAS.Results.csv",sep=""), head = TRUE)
  #CHR #column 2
  #POS #3
  #SNP #1
  #colnames(RESULTS)[1]<-"MARKER"
  #N	 #6
  colnames(RESULTS)[6]<-"N"
  #EFFECT_ALLELE
  RESULTS$EFFECT<-sign(RESULTS$effect)
  RESULTS$EFFECT[which(RESULTS$EFFECT==1)]<-"+"
  RESULTS$EFFECT[which(RESULTS$EFFECT==-1)]<-"-"
  #NON_EFFECT_ALLELE
  #RESULTS$NON_EFFECT_ALLELE<-"T"
  #EFFECT_ALLELE_FREQ	
  #RESULTS$EFFECT_ALLELE_ALLELE<-NA
  #BETA	#9
  colnames(RESULTS)[9]<-"BETA"
  #SE	
  #RESULTS$SE<-NA
  #r2	#7
  colnames(RESULTS)[7]<-"r2"
  #r2hat #8
  colnames(RESULTS)[8]<-"r2hat"
  #P_VAL #4
  colnames(RESULTS)[4]<-"PVALUE"
  colnames(RESULTS)[10]<-"FDR"
  write.table(RESULTS,paste("GAPIT.MLM.",n,".GWAS.RESULTS.txt",sep=""),row.names = F,quote = FALSE)
#print(length(which(RESULTS$`FRD Adjusted Pvalue` <=0.15)))
}

