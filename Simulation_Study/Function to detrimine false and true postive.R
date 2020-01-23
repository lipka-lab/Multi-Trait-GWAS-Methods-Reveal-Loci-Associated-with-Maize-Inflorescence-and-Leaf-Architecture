#Function for detriming if detected SNP is false or true positive
 #also removes casual mutation from GWAS output so they are not considered for detection
simulation.success<-function(c,GWAS.RESULTS,chr,pos,rs,pvalue,FDR,LD.measurement,
                             myGM,myGD,Threshold,Alpha,bp,LD.logic,
                             removeCasual,NULL.trait){
  if(NULL.trait){
    if(Threshold=="FDR"){sigSNPs<-GWAS.RESULTS[which(GWAS.RESULTS[,FDR]<Alpha),]}
    if(Threshold=="bonferoni"){sigSNPs<-GWAS.RESULTS[which(GWAS.RESULTS[,pvalue]<Alpha/nrow(GWAS.RESULTS)),]}
    if(Threshold=="N"){sigSNPs<-GWAS.RESULTS[1:Alpha,]}
    breif<-nrow(sigSNPs)
  } else{
  if(removeCasual){
    for(k in 1:nrow(c)){
      print(paste("Removing Marker",c$snp[k],"from GWAS RESULTS"))
    if(length(which(GWAS.RESULTS[,rs]==as.character(c$snp[k])))>0){
    GWAS.RESULTS<-GWAS.RESULTS[-which(GWAS.RESULTS[,rs]==as.character(c$snp[k])),]
    } else { print(paste("!!!!! Marker",c$snp[k],"not found in GWAS RESULTS !!!!!"))}
    print(nrow(GWAS.RESULTS))
    }
  } else {print("Leaving Casual Mutations in GWAS RESULTS")}
  if(Threshold=="FDR"){sigSNPs<-GWAS.RESULTS[which(GWAS.RESULTS[,FDR]<Alpha),]}
  if(Threshold=="bonferoni"){sigSNPs<-GWAS.RESULTS[which(GWAS.RESULTS[,pvalue]<Alpha/nrow(GWAS.RESULTS)),]}
  if(Threshold=="N"){sigSNPs<-GWAS.RESULTS[1:Alpha,]}
  #if(nrow(sigSNPs)>0){
  
  
  if(nrow(sigSNPs)>0){#if no signficant SNPs to test we just skip it all
   print(paste("Getting summary for",nrow(sigSNPs),"SNPs"))
    breif<-as.data.frame(matrix(NA,nrow = nrow(sigSNPs),ncol=3))
    colnames(breif)<-c("SNP","bp Distance",LD.measurement)
  for(i in 1:nrow(sigSNPs)){
    print(paste("SNP",i))
    if(bp){
      pos<-myGM[which(myGM$SNP==as.character(unlist(sigSNPs[i,rs]))),]
    if(length(which(c$chr==pos$Chromosome))>0){
    distance<-min(abs(as.numeric(as.character(pos$Position))-as.numeric(as.character(casualmutations$pos[which(casualmutations$chr==pos$Chromosome)]))))  
    } else {distance <- NA}
    breif$SNP[i]<- as.character(unlist(sigSNPs[i,rs]))
    breif$`bp Distance`[i]<-distance
    }
    if(LD.logic){
 if(length(which(c$chr==pos$Chromosome))>0){
      #get numeric SNP info for signficant snps and convert to A/T
      SigSNPsGeno<-myGD[,which(colnames(myGD)==pos$SNP)]
      #un-numericalize genotypes for function LD()
      SigSNPsGeno[which(SigSNPsGeno==0)]<-"A/A"
      SigSNPsGeno[which(SigSNPsGeno==1)]<-"A/T"
      SigSNPsGeno[which(SigSNPsGeno==2)]<-"T/T"
      #colnames(SigSNPsGeno)<-as.character(unlist(pos$SNP[i]))
      #convert info for casual mutations to A/T
      temp<-which(c$chr==pos$Chromosome)
      temp<-c[c(temp),]
tempSNP<-temp$snp[which.min(abs(as.numeric(as.character(pos$Position))-as.numeric(as.character(c$pos[which(c$chr==pos$Chromosome)]))))] 
      
      genos<-myGD[,which(colnames(myGD)==tempSNP)]
      #if(ncol(genos)>0){
       genos[which(genos==0)]<-"A/A"
       genos[which(genos==1)]<-"A/T"
       genos[which(genos==2)]<-"T/T"
       tempG<-cbind(SigSNPsGeno,genos)
       LD<-LD(makeGenotypes(tempG))
       
     
  if(LD.measurement=="D"){LD.temp<-unlist(LD[2])}
    if(LD.measurement=="D'"){LD.temp<-unlist(LD[3])}
    if(LD.measurement=="r"){LD.temp<-unlist(LD[4])}
    if(LD.measurement=="R2"){LD.temp<-unlist(LD[5])}
    if(LD.measurement=="X2"){LD.temp<-unlist(LD[7])}
          breif[i,3]<-LD.temp[3]
 #} else {print(paste("Marker",tempSNP,"not found in myGD"))}
      } 
      } else {LD <- NA}
    

      }
  return(breif)
  print(breif) } 
  else{print("No SNPs passes threshold")
  breif<-as.data.frame(matrix(NA,nrow = 1,ncol=3))
  colnames(breif)<-c("SNP","bp Distance",LD.measurement)
  breif$SNP<- NA
  breif$`bp Distance`<-NA
  breif[1,3]<-NA
  }
  return(breif)}
  }



  