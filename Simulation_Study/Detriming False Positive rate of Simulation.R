library(data.table)
library(genetics)
library(GAPIT3)
source("Function to detrimine false and true postive.R")
#setRepositories(ind = 1:2)
#devtools::install_github("samuelbfernandes/simplePHENOTYPES")
library(simplePHENOTYPES)

#INPUTS
#genotypic input
# Numeric genotype files (myGM, myGD, myPC, myKI)
load("SNP55kv2_282_GenoFiles_GapitInput.RData")
rm(myKI) #remove this version of myKI it is calculated using all 10 chromosome
rm(myGM)
rm(myGD)
#10 kinship files in a list ordered (1:10 coresponding to chr 1:10)
load("SNP55kv2_282_individual_chr_kinship_files_filter0.05.RData")
load("GBSv4_282_kinship_per_chromosome.RData")
KinNames<-read.csv("GAPIT.GBS.Kin.VanRaden.WGSnames.csv",header = F)
Names<-KinNames[,1]
for(i in 1:10){
  kinshiplist[[i]]$V1<-Names
}
save(kinshiplist,file="GBSv4_282_kinship_per_chromosome_WGSnames.RData")
genoforsimulation<-data.frame(SNP=myGM$SNP,allele=as.character(NA),chr=myGM$Chromosome,pos=myGM$Position,cm=NA,t(myGD[,2:ncol(myGD)]))
alleles<-fread("SNP55K_maize282_AGPv2.hmp.txt")
alleles<-alleles[,1:2]
test<-merge(alleles,genoforsimulation,by.x="rs",by.y="SNP")
genoforsimulation<-test[,-3]
colnames(genoforsimulation)[1:5]<-c("snp", "allele", "chr", "pos", "cm")
colnames(genoforsimulation)[6:ncol(genoforsimulation)]<-as.character(unlist(myGD$taxa))

gemma<-TRUE
Null.trait<-FALSE
RDATA<-FALSE
ntraits = 3
rep = 100
output_format = "wide"
geno_obj = genoforsimulation
add_QTN_num = 3
h2 = c(0.5, 0.1, 0.1)
add_effect = c(0.9,0.5,0.2)
vary_QTN = TRUE
architecture = "partially"
pleio_a = 3
trait_spec_a_QTN_num=c(2,2,2)
output_dir = "Setting_test"
to_r = TRUE
seed = 10
model = "A"
sim_method = "geometric"

if(RDATA){
  sim.dir<-getwd()
  setwd(paste(sim.dir,output_dir,sep="/"))
  load("Simulate.Rdata")
} else {
  #simulate pheno data
  simulation<-  create_phenotypes(
    geno_obj = genoforsimulation,
    add_QTN_num = add_QTN_num,
    #dom_QTN_num = dom_QTN_num,
    h2 = h2,
    add_effect = add_effect,
    pleio_a = pleio_a,
    trait_spec_a_QTN_num=trait_spec_a_QTN_num,
    #dom_effect = c(0.0001,0.0002,0.0003),
    ntraits =  ntraits,
    rep =rep,
    vary_QTN = vary_QTN,
    output_format = output_format,
    quiet=T,
    architecture = architecture,
    output_dir = output_dir,
    home_dir = getwd(),
    to_r = TRUE,
    seed = 10,
    model = "A",
    sim_method = "geometric")
  
  sim.dir<-getwd()
  setwd(paste(sim.dir,output_dir,sep="/"))
  save(simulation,file="Simulate.Rdata")
}
#Add in step to repeat simple phenotypes with same seed number but gemma output in a folder named gemma 
if(gemma){
  
  output_dir_gemma ="gemma"
  output_format_gemma = "gemma"
  
  simulation<-  create_phenotypes(
    geno_obj = genoforsimulation,
    add_QTN_num = add_QTN_num,
    #dom_QTN_num = dom_QTN_num,
    h2 = h2,
    add_effect = add_effect,
    pleio_a = pleio_a,
    trait_spec_a_QTN_num=trait_spec_a_QTN_num,
    #dom_effect = c(0.0001,0.0002,0.0003),
    ntraits =  ntraits,
    rep =rep,
    vary_QTN = vary_QTN,
    output_format = output_format_gemma,
    quiet=T,
    architecture = architecture,
    output_dir = output_dir_gemma,
    home_dir = getwd(),
    to_r = TRUE,
    seed = 10,
    model = "A",
    sim_method = "geometric")
}


#the following work if simplyphenotypes was just run and the data is saved in the environment
#need to arrange data as a list
phenolist<-vector(mode = "list", length = rep)




for(r in 1:rep){
  if(r==1){temp.rep<-simulation[which(simulation$Rep==r),c(1:ntraits,ntraits+1)]
  colnames(temp.rep)<-paste(colnames(temp.rep),"_rep_",r,sep="")
  PC<-prcomp(temp.rep[2:ncol(temp.rep)],center=T,scale=T)
  PC.temp<-data.frame(taxa=temp.rep$`<Taxa>`,PC$x)
  colnames(PC.temp)<-paste(colnames(PC.temp),"_rep_",r,sep="")
  traits<-data.frame(Trait=colnames(temp.rep[c(2:ntraits,ntraits+1)]) ,Rep=r)
  PC.traits<-data.frame(Trait=colnames(PC.temp[c(2:ntraits,ntraits+1)]) ,Rep=r)
  }
  else{
    temp<-simulation[which(simulation$Rep==r),c(1:ntraits,ntraits+1)]
    colnames(temp)<-paste(colnames(temp),"_rep_",r,sep="")
    PC<-prcomp(temp[2:ncol(temp)],center=T,scale=T)
    temp.PC<-data.frame(taxa=temp.rep$`<Taxa>`,PC$x)
    colnames(temp.PC)<-paste(colnames(temp.PC),"_rep_",r,sep="")
    temp.rep<-merge(temp.rep,temp,by=1)
    PC.temp<-merge(PC.temp,temp.PC,by=1) 
    traits.temp<-data.frame(Trait=colnames(temp[c(2:ntraits,ntraits+1)]) ,Rep=r)
    traits<-rbind(traits,traits.temp)
    PC.traits.temp<-data.frame(Trait=colnames(temp.PC[c(2:ntraits,ntraits+1)]) ,Rep=r)
    PC.traits<-rbind(PC.traits,PC.traits.temp)
  } 
  
  print(r)}
traits<-rbind(traits,PC.traits)


#make file for GAPIT with both raw traits and PCs of traits
phenoforGAPIT<-merge(temp.rep,PC.temp,by=1)

#pipelines options
bp<-TRUE  #record base pair (bp) distance between significant peaks and causal variants
LD.logic<-TRUE #record LD between significant peaks and causal variants 
Threshold<-"FDR" #options: "FDR" "bonferoni" "N" 
#- N will just test top N most signficant SNPS reguardless of significance
#Bascially ask the question - is GWAS at least correctling ordering SNPs
Alpha<-0.1#detrimine significance threshold # set to whole number if N threshold chosen
removeCasual<-FALSE #if TRUE function removes causal simulated SNPs from GWAS output
LD.measurement<-"R2" #options: D  D' r  R2 X2
NULL.trait<-TRUE #If NULL.trait is TRUE than will only report number of signficant SNPs
#Null traits have no casual mutations and hert 0

#geno info for casual mutations
if(Null.trait==FALSE){
  casualmutations<-fread(paste("geno_info_for_",rep,"_reps_and_",add_QTN_num,"_Add_QTN.txt",sep=""))
}
colnames(phenoforGAPIT)<-c("Taxa",colnames(phenoforGAPIT[,2:ncol(phenoforGAPIT)]))
#run LOCO GWAS pipeline on every single trait for each replicate
results.dir<-"Univariate GWAS and PC GWAS results"
dir.create(results.dir)
setwd(paste(getwd(),results.dir,sep="/"))
home.dir<-getwd()

for(i in 1:10){
  #step 4 run GWAS
  setwd(home.dir)
  dir.create(paste("Chr",i,sep=""))
  setwd(paste(home.dir,"/Chr",i,sep=""))
  myKI<-kinshiplist[[i]]
  colnames(myKI)[1]<-"taxa"
  GM.K<-myGM[which(myGM$Chromosome==i),]
  GD.K<-myGD[,c(1,(which(myGM$Chromosome==i)+1))]
  
  #Run GAPIT
  myGAPIT <- GAPIT( #file.fragment = file.fragment,#Geno.View.output=FALSE #PCA.total=3,
    Y=phenoforGAPIT[c(1,2:ncol(phenoforGAPIT))], CV=myPC, GD=GD.K, GM=GM.K, K=myKI, SNP.fraction=1,PCA.total = 5, 
    group.from = 280, group.to = 280, Model.selection = F, SNP.MAF=0.05, model="MLM",Geno.View.output = F)
  
}
#[,c(1,534:ncol(phenoforGAPIT))]
#combined chromosome results for each trait and get new FDR
#during this loop also calculate false/true positive nature of
#GWAS.RESULTS
#create file that includes each trait name and the corresponding
# rep number to match with casual mutation file
#trait_reps<-strsplit(traits, "_")
summary<-as.data.frame(matrix(NA,nrow = 1,ncol = 5))
setwd(home.dir)
dir.create(paste("Complete GWAS Results",sep=""))
#read back in each chr gwas results for each traits
for(j in 1:2){#nrow(traits)
  t<-traits$Trait[j]
  r<-traits$Rep[j]
  print(paste(t))
  for(i in 1:10){ print(paste("reading in chr",i))
    setwd(paste(home.dir,"/Chr",i,sep=""))
    if(i==1){#combine GWAS result
      RESULTS<- read.csv(paste("GAPIT.MLM.",t,".GWAS.Results.csv",sep=""), head = TRUE)
      print(unique(RESULTS$Chromosome))
    } else {
      GWAS.RESULTS<-read.csv(paste("GAPIT.MLM.",t,".GWAS.Results.csv",sep=""), head = TRUE)
      GWAS.RESULTS$Chromosome<-i
      RESULTS<-rbind(GWAS.RESULTS,RESULTS)    
      print(unique(RESULTS$Chromosome))
    } }
  #remove FDR column
  RESULTS.new<-RESULTS[,-9]
  #calculate new FDR
  Conduct.FDR <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = RESULTS.new,FDR.Rate = 0.05, FDR.Procedure = "BH")
  GWAS.Results.with.FDR<- Conduct.FDR$PWIP
  colnames(GWAS.Results.with.FDR)[9:10]<-c("effect","FRD Adjusted Pvalue")
  setwd(paste(home.dir,"/Complete GWAS Results",sep=""))
  write.table(GWAS.Results.with.FDR, paste("GAPIT.MLM.",t, ".GWAS.RESULTS.csv", sep = ""), 
              quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
  
  #run simulation.success on GWAS results for each trait and replicate
  c<- casualmutations[which(casualmutations$rep==j)]
  
  summary.temp<-data.frame(Trait=t,Rep=r,simulation.success(c,GWAS.RESULTS,myGM,Threshold,Alpha,bp,LD.logic,removeCasual,NULL.trait))
  colnames(summary)<-colnames(summary.temp)
  summary<-rbind(summary,summary.temp)
}


