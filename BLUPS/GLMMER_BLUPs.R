############
#October 1st
#Brian Rice
  library(bestNormalize)
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

#Step one : data processing 
data<-read.delim2("traitMatrix_maize282NAM_v15-130212.txt",header=T,stringsAsFactors = FALSE)
data<-data[-1,]
#Pipeline options
#Select Genotype file to make sure the correct names are used
SNPCHP55k<-TRUE
GBS<-FALSE
GBSnames<-read.csv("~/Box/GermplasmMaterial/Reference GBS taxanames.csv",header=T,stringsAsFactors = FALSE)
#tranform raw values?
tranform_raw<-TRUE
#get BLUPS using 282 and NAM?
NAM<-TRUE
#if checking normality and transforming raw values
normality_1<-TRUE
#if checking normality and only transforming non normal BLUPs from step 1 is desired
#FALSE will normalize all traits
normality_2<-TRUE


if(SNPCHP55k){
#genotype files for GWAS later and to confirm taxa names match
load("~SNP55kv2_282_GenoFiles_GapitInput.RDS")
#change names to original phenotype file
data[which(data[,1]=="Ia5125"),1]<-"IA5125"
data[which(data[,1]=="W22_R-r:std"),1]<-"W22_R.r.std"
data[which(data[,1]=="Ki14"),1]<-"KI14"
data[which(data[,1]=="Oh603"),1]<-"OH603"
data[which(data[,1]=="Va14"),1]<-"VA14"

#change GBSnames file names
GBSnames[which(GBSnames$TAXA=="Ia5125"),1]<-"IA5125"
GBSnames[which(GBSnames$TAXA=="W22_R-r:std"),1]<-"W22_R.r.std"
GBSnames[which(GBSnames$TAXA=="Ki14"),1]<-"KI14"
GBSnames[which(GBSnames$TAXA=="Oh603"),1]<-"OH603"
GBSnames[which(GBSnames$TAXA=="Va14"),1]<-"VA14"
}

#merge with GBS names to reduce to 281 individuals
  #remove individuals without GBS names
    #if column = -1 (GBS names) if column  = -2 non GBS names
if(NAM==FALSE){
  data282<-merge(GBSnames[,c(1,4)],data,by=1)#,all.y = T

if(SNPCHP55k){
  data282<-data282[-which(is.na(data282$GBS.V4names)==T),-2]
} else {
  data282<-data282[-which(is.na(data282$GBS.V4names)==T),-1]
}} else { data282<-data}
colnames(data282)[1]<-"TAXA"

#checking taxa names in phenotype file
#test<-merge(data282[,1:2],myGD[,1:2],by.x=1,by.y=1,all=T)

#change NaN into NA
for(i in 2:ncol(data282)){
  data282[which(data282[,i]=="NaN"),i]<-NA
  data282[,i]<-as.numeric(as.character((data282[,i])))
}

#remove environments with no data for particular phenotypes
  #remove environment if number of obervations is less than 3
cols_to_rm<-c()
for(y in 2:ncol(data282)){
  if(length(which(is.na(data282[,y])==F))<=3){
  #print(y)
  cols_to_rm<-c(cols_to_rm,y)
  }}
data282<-data282[,-cols_to_rm]

#subset to only the traits I want
#I have a file called trait_names which has the replication number removed from the colname title
  #this allows us to sort traits based on column name
c_names<-read.csv("trait_names.csv",header=F)  #check for normality
colnames(data282)<-c_names$V1
traits<-read.csv("BLUP Model Trait Pipeline Summaries.csv",header = T)

#step 2 : evalutate normality
  #a) normaility of each environment
#and
  #b) normality of residuals 
if(normality_1){
Normality<-as.data.frame(matrix(NA,1,3))
colnames(Normality)<-c("trait","w","p")
residualNormality<-as.data.frame(matrix(NA,1,2))
colnames(residualNormality)<-c("trait","p")

for(t in traits$trait){
  #subset data282 into a single trait across replications
  traitdata<-data282[,c(1,which(colnames(data282)==t))]
  
  #section to check normality of the residuals of a linear model
  dataforLMER<-NULL
  for (i in 2:ncol(traitdata)){
    E<-rep(i-1,nrow(traitdata))
    P<-as.numeric(as.character(traitdata[,i]))
    EP<-cbind.data.frame(traitdata[,1],P,E)
    dataforLMER<-rbind.data.frame(dataforLMER,EP)
  }  
  colnames(dataforLMER)[1]<-"taxa"
 myModel<-lmer(P~(1|taxa)+(1|E),data=dataforLMER,REML=TRUE,
                control = lmerControl(optimizer ="Nelder_Mead"))

  s.res<-shapiro.test(resid(myModel)[sample (c(1:5707), size=5000, replace =F)])
  pdf(paste(t,"Panzea NAM and 282 LMER Residual Distributions.pdf"))
  
  hist(resid(myModel),breaks=1000,main=t,
       xlab=paste("S.W. ",s.res$p.value,sep=""),
       cex.main=0.8,
       cex.lab=0.8)
  dev.off()
  #Draw the qq-plot of the normally distributed data using pch=19 to produce solid circles.
  #qqnorm(dataforLMER$P,main="QQ plot",pch=19)
  #Add a line where x = y to help assess how closely the scatter fits the line.
  #qqline(dataforLMER$P)
  #library(nortest)
  #ad.test(myModel)$p.value
  residualNormality<-rbind(residualNormality,c(t,s.res$p.value))
  
  print(paste("Trait ",t," has shapiro wilks pvalue of ",s.res$p.value,sep=""))
  
    pdf(paste(t,"Panzea 282 Raw Data Distributions.pdf"))
  par(mfrow=c(2, 5))
  for(i in 2:ncol(traitdata)){
    #shapiro wilks
    s<-shapiro.test(traitdata[sample(c(1:5707), size=5000, replace =F),i])
    #put hist for each rep onto a single pdf
    #add text for pvalue of s.w.
    hist(traitdata[,i],main=paste(t," Rep",i-1,sep=""),xlab=paste("S.W. ",s$p.value,sep=""),
         cex.main=0.8,
         cex.lab=0.8)
    Normality<-rbind(Normality,c(t,i-1,s$p.value))
  } 
  dev.off()
}
Normality<-na.omit(Normality)
residualNormality<-na.omit(residualNormality)
write.csv(Normality,paste(t," Shapiro.Wilk for each environment.csv"))
write.csv(residualNormality,paste(t," Shapiro.Wilk for LMER residuals.csv"))

}

#Step 3 compute BLUPS using model based on normaility assesment of step 2
library(glmm)
library(lme4)
BLUPS<-data.frame(taxa=data282$TAXA)
for(j in 1:nrow(traits)){ #c(13,17,19,23)
  print(paste(traits$trait[j]))
  #put all the environments into 3 columns (taxa,environment and value)
  #j<-22
  traitdata<-data282[,c(1,which(colnames(data282)==traits$trait[j]))]
  dataforGLM<-NULL
  for (i in 2:ncol(traitdata)){
    E<-as.factor(rep(i-1,nrow(traitdata)))
    P<-as.numeric(as.character(traitdata[,i]))
    print(length(which(is.na(P)==TRUE)))
    EP<-cbind.data.frame(traitdata[,1],P,E)
    dataforGLM<-rbind.data.frame(dataforGLM,EP)
  }  
  dataforGLM<-na.omit(dataforGLM)
  colnames(dataforGLM)[1]<-"taxa"
  print(paste("Number of environments = ",length(unique(dataforGLM$E))))
  #dataforGLM$E<-as.factor(dataforGLM$E)
  if(paste(traits$trait[j])=="EarWeight"){dataforGLM$P[which(dataforGLM$P==0)]<-NA
   print("Ear weight has a phenotype of 0 which makes no sense so I changed it to NA")
}
  if(paste(traits$trait[j])=="EarDiameter"){dataforGLM$P[which(dataforGLM$P==0)]<-NA
  print("Ear diameter has a phenotype of 0 which makes no sense so I changed it to NA")
  }
  
  #if(paste(traits$Link[j])=="poisson"){dataforGLM$P[which(dataforGLM$P==0)]<-NA}
  ####
  #Run model depending on Link function determined in item traits
  if(paste(traits$Link[j])=="gaussian"){
    print(paste(traits$trait[j]," Gaussian",sep=""))
    myModel<-lmer(P~(1|taxa)+(1|E),data=dataforGLM,REML=TRUE,
                  control = lmerControl(optimizer ="Nelder_Mead"))
  } 
  if(paste(traits$Link[j])=="poisson"){
    m1<-glmer(as.integer(P)~(1|taxa)+(1|E),family="poisson",data=dataforGLM)
    print(paste(traits$trait[j]," poisson i.e. negative binomial",sep=""))
   myModel<- glmer.nb(as.integer(P)~(1|taxa)+(1|E),data=dataforGLM,
                      #na.action=na.exclude,
                      control=glmerControl(optCtrl=list(maxfun=100000)))
   print(overdisp_fun(myModel))
  }
  if(paste(traits$Link[j])=="Gamma"){
    print(paste(traits$trait[j]," Gamma",sep=""))
  myModel<-glmer(P~(1|taxa)+(1|E),family="Gamma",data=dataforGLM,
            control=glmerControl(optCtrl=list(maxfun=100000)))
  }
  if(paste(traits$Link[j])=="transform"){
    print(paste(traits$trait[j]," transformed",sep=""))
    #tranform using 
    dataforGLM<-NULL
    for (i in 2:ncol(traitdata)){
      E<-rep(i-1,nrow(traitdata))
      P1<-as.numeric(as.character(traitdata[,i]))
      P<-predict(orderNorm(P1))
     print(shapiro.test(P[sample(c(1:5707), size=300, replace =F)])$p.value)
      EP<-cbind.data.frame(traitdata[,1],P,E)
      dataforGLM<-rbind.data.frame(dataforGLM,EP)
    }  
    colnames(dataforGLM)<-c("taxa","P","E")
    if(paste(traits$Model.Notes[j])=="boundary (singular) fit"){
    print("Removing Environment from model")
      myModel<-lmer(P~(1|taxa),data=dataforGLM,REML=TRUE,
                    control = lmerControl(optimizer ="Nelder_Mead"))
    } else {
      myModel<-lmer(P~(1|taxa)+(1|E),data=dataforGLM,REML=TRUE,
                  #na.action=na.exclude,
                  control = lmerControl(optimizer ="Nelder_Mead"))
    }
    
  s.res<-shapiro.test(resid(myModel)[sample (c(1:5707), size=300, replace =F)])
    pdf(paste(traits$trait[j],"Panzea NAM and 282 LMER Residual Distributions Post Transformation.csv"))
    hist(resid(myModel),breaks=1000,main=traits$trait[j],
         xlab=paste("S.W. ",s.res$p.value,sep=""),
         cex.main=0.8,
         cex.lab=0.8)
    dev.off()
    }
  temp<-coef(myModel)$taxa
  colnames(temp)<-traits$trait[j]
  BLUPS<-merge(BLUPS,temp,by.x=1,by.y="row.names",all=T)
  
}
print("Its ok to see warnings()")
  print(warnings()[1])
Normality_BLUPS<-as.data.frame(matrix(NA,1,3))
colnames(Normality_BLUPS)<-c("trait","p")

#step 4 compute normailty of model BLUPS
  # and transform them using normal qunatile transformation
for(w in 2:ncol(BLUPS)){
  print(colnames(BLUPS)[w])
  s<-(shapiro.test(BLUPS[sample(c(1:5706), size=1000, replace =F),w])$p)
  print(s)
  hist(as.numeric(unlist(BLUPS[w])),breaks=100,main=colnames(BLUPS)[w],xlab=s)
  Normality_BLUPS<-rbind(Normality_BLUPS,c(colnames(BLUPS)[w],s))
  
  }
Normality_BLUPS<-na.omit(Normality_BLUPS)

#normal qunatile transformation
if(normality_2){ #if TRUE will only transform non normal data
  BLUPS_postQNORM<-as.data.frame(matrix(NA,1,2))
  colnames(BLUPS_postQNORM)<-c("trait","p")
  normBLUPS<-as.data.frame(matrix(NA,nrow(BLUPS),ncol(BLUPS)))
  normBLUPS[,1]<-BLUPS$taxa
for(w in 2:24){
  if(traits$Tranform.BLUPs=="Yes"){
    normBLUPS[,w]<- predict(orderNorm(BLUPS[,w]))
    print(paste("Normalized trait",Normality_BLUPS$trait[w],sep=" "))
  } else { normBLUPS[,w]<-BLUPS[,w] 
  print(paste("Did not normalized trait",Normality_BLUPS$trait[w-1],sep=" "))
  }}} else {
    BLUPS_postQNORM<-as.data.frame(matrix(NA,1,2))
    colnames(BLUPS_postQNORM)<-c("trait","p")
    normBLUPS<-as.data.frame(matrix(NA,281,24))
    normBLUPS[,1]<-BLUPS$taxa
    print(paste("Normalized all traits"))
    
    for(w in 2:24){
        normBLUPS[,w]<- predict(orderNorm(BLUPS[,w]))
      }
  }
#You will get warnings - basically the function is letting you know it does garuntee normailty

#write this version of BLUPS to csv
colnames(normBLUPS)<-c(colnames(BLUPS))
write.csv(normBLUPS,"enter you working directory",row.names=FALSE)

#some taxa still have missing data due to missing data in all environments
#you can check to what extent with this
for(i in 2:ncol(normBLUPS)){
  print(length(which(is.na(normBLUPS[,i])==T)))
}
BLUPS282<-merge(GBSnames[,c(1,4)],normBLUPS,by=1)
BLUPS282<-BLUPS282[,-1] #-2 will remove gbs names, -1 will remove non gbs names
for(i in 2:ncol(BLUPS282)){
  #print(i-1)
  print(length(which(is.na(BLUPS282[,i])==T)))
}
