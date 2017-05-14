#Code for Stevens and Latimer 2015
#Jens Stevens stevensjt@gmail.com 11/21/14
##The purpose of this script is to build statistical models of each transition probability that account for treatment effects of canopy closure, fire and snowpack (where these treatments have a meaningful impact on a given demographic transition). The script generates a table "bst" for each species, which contains the parameter estimates on a logit scale for each transition probability, along with associated variance for that parameter estimate. These tables are then saved and used in the demographic analyses in the CanopyProjections.R and FireProjections.R scripts.

##To start: Ready to go. Run this as per the workflow. Then move on to the Projections scripts and delete this line.

library(plyr) #for ddply(). verison 1.8.4
library(MuMIn) #for AICc(). version 1.15.6
library(rethinking) #for sample.naive.posterior(). version 1.55. Not on CRAN.
library(lme4) #for glmer(). version 1.1-12
library(ez) #for ezPredict(). version 4.4-0
#library(boot)
#library(MASS)

stdErr <- function(x) {#Function to generate standard error
  sd(na.omit(x))/ sqrt(length(na.omit(x)))} 


####1. Initializations and variable definition####
spp <- "sj" #Define species as Cytisus scoparius (cs; run first) or Spartium junceum (sj)
d <- read.csv("data/Individual.csv")
d <- d[d$Species==spp,] #Subset data to species in question
d$Canopy <- factor(d$Canopy,levels=c("O","T","C","D"))
Canopy <- c("O","T","C","D")
Scenarios <- c("O.S10","O.S40","O.S70",
               "T.UB.S10","T.UB.S40","T.UB.S70",
               "T.B.S10","T.B.S40","T.B.S70",
               "C.S10","C.S40","C.S70",
               "D.S10","D.S40","D.S70") #Canopy.Burned/UnBurned/Snow (inches). For the MS, ran rougly approximate simulations in cm (25,100,175), but kept the 10-40-70 convention in the code.

if(spp=="cs"){ #Define template table to hold best-model output parameters ("bst")
  bst = read.csv("data/csParamTableBST_Template.csv")
} else{
  bst = read.csv("data/sjParamTableBST_Template.csv")}
m.out.l <- m.out <- plots <- list()

####2.Data generation function####
#Sets the response variable as a given transition stage. Generates data representing the appropriate population to estimate the given transition. Returns "md" (="model data").
GenerateData=function(P,d){
  if(P=="Pg0"){
    d$y=ifelse(is.na(d$F1),0,1) # Response variable is germination success in Fall 0
  }
  if(P=="Pg1"){
    d=d[is.na(d$F1),]
    d$y=ifelse(is.na(d$Sp1),0,1) # Response variable is germination success in Spring 1
    d$SnowTot=d$S1Tot
  }
  if(P=="Pg2"){
    d=d[is.na(d$Sp1),]
    d$y=ifelse(is.na(d$Sp2),0,1) # Response variable is germination success in Spring 1
    d$SnowTot=d$S2Tot
  }
  if(P=="PgB"){
    d=d[is.na(d$Sp2),]
    d$y=ifelse(is.na(d$Sp3),0,1) # Response variable is germination success in Spring 1
    d$SnowTot=d$S3Tot
  }
  if(P=="Pss1"){#1st Summer survival
    CandidateRows1=which(d[,"Sp1"]==1 & is.na(d[,"F1"]))
    CandidateRows2=which(d[,"Sp2"]==1 & is.na(d[,"F2"]))
    d[CandidateRows1,"SnowTot"]=d[CandidateRows1,"S1Tot"]
    d[CandidateRows2,"SnowTot"]=d[CandidateRows2,"S2Tot"]
    d[CandidateRows1,"y"]=d[CandidateRows1,"F2"]
    d[CandidateRows2,"y"]=d[CandidateRows2,"F3"]    
    d=d[sort(c(CandidateRows1,CandidateRows2)),]
  }
  if(P=="Pws0"){#Winter survival for fall germinants (relevant for Spanish broom only)
    CandidateRows1=which(d[,"F1"]==1 )
    d[CandidateRows1,"SnowTot"]=d[CandidateRows1,"S1Tot"]
    d[CandidateRows1,"y"]=d[CandidateRows1,"Sp1"]
    d=d[sort(c(CandidateRows1)),]
  }  
  if(P=="Pws1"){#1st Winter survival
    CandidateRows1=which(d[,"F1"]==1 )
    CandidateRows2=which(d[,"F2"]==1 & is.na(d[,"F1"]))
    CandidateRows3=which(d[,"F3"]==1 & is.na(d[,"F2"]))
    d[CandidateRows1,"SnowTot"]=d[CandidateRows1,"S1Tot"]
    d[CandidateRows2,"SnowTot"]=d[CandidateRows2,"S2Tot"]
    d[CandidateRows3,"SnowTot"]=d[CandidateRows3,"S3Tot"]
    d[CandidateRows1,"y"]=d[CandidateRows1,"Sp1"]
    d[CandidateRows2,"y"]=d[CandidateRows2,"Sp2"]    
    d[CandidateRows3,"y"]=d[CandidateRows3,"Sp3"]    
    d=d[sort(c(CandidateRows2,CandidateRows3)),]
  }
  if(P=="Pss2"){#2nd Summer survival
    CandidateRows1=which(d[,"Sp1"]==1 & d[,"Sp2"]==1)
    d[CandidateRows1,"SnowTot"]=d[CandidateRows1,"S2Tot"]
    d[CandidateRows1,"y"]=d[CandidateRows1,"F3"]
    d=d[sort(c(CandidateRows1)),]
  }
  if(P=="Pws2"){#2nd Winter survival
    CandidateRows1=which(d[,"F1"]==1 & d[,"F2"]==1)#Alive in 1st and 2nd fall
    CandidateRows2=which(is.na(d[,"F1"])& d[,"F2"]==1 & d[,"F3"]==1)#Alive in 2nd and 3rd fall
    d[CandidateRows1,"SnowTot"]=d[CandidateRows1,"S2Tot"]
    d[CandidateRows2,"SnowTot"]=d[CandidateRows2,"S3Tot"]
    d[CandidateRows1,"y"]=d[CandidateRows1,"Sp2"]
    d[CandidateRows2,"y"]=d[CandidateRows2,"Sp3"]    
    d=d[sort(c(CandidateRows1,CandidateRows2)),]
  }
  if(P=="Pm3"){#Maturity in 3rd growing season (Maybe do just for cs?)
    CandidateRows1=which(d[,"Sp1"]==1 & d[,"Sp3"]==1)#Alive in 1st and 3rd spring
    d[CandidateRows1,"SnowTot"]=d[CandidateRows1,"S3Tot"] #All flowering individuals were in their 3rd year.
    d[CandidateRows1,"y"]=d[CandidateRows1,"Su3.Fl"]
    d=d[sort(CandidateRows1),]
    d$y=ifelse(is.na(d$y),0,1) #Convert NA's to 0's ("not mature")
  }
  if(P=="S1"){#Seed set in 3rd growing season
    CandidateRows1=which(!(is.na(d[,"Su3.Fl"]))) #Which individuals were flowering?
    d[CandidateRows1,"SnowTot"]=d[CandidateRows1,"S3Tot"] 
    d[CandidateRows1,"y"]=d[CandidateRows1,"Su3.Fl"]
    d=d[sort(CandidateRows1),]
    d$y=d$Su3.Sd
  }
  return(d)
}

####3.Statistical multilevel models####
####3a. Function definitions####
#might be some glitches in model structure on Canopy Effects, some wonkiness related to choosing to do multilevel models in Snow Effects
CanopyEffects=function(P,md,subset=""){
  #It is justified to not subset the data for canopy effects, because if a snow effect subsequently becomes significant, it overrides the main canopy effect anyway.
  #glm is the proper model here because the variation in canopy is *at* the block level
  md.glm=ddply(md,.(Block),summarize,Canopy=unique(Canopy),y=sum(y),n=length(Block),p=y/n)
  m0=glm(p~1, data=md.glm, weights=n, family="binomial")
  m1=glm(p~Canopy, data=md.glm, weights=n, family="binomial") 
  m1b=glm(p~0+Canopy,data=md.glm,weights=n,family="binomial") #Coding the model this way (0+) makes predictions easier; this is the same model as m1.
  daicc=AICc(m0)-AICc(m1)
  d.pred=data.frame(Canopy=c("O","T","C","D"))
  d.pred$pred=predict(m1,d.pred)#Generate predicted values for each canopy level
  post=sample.naive.posterior(m1b,n=10000)#Sample posterior density of canopy model fit
  post0=sample.naive.posterior(m0,n=10000)#Sample posterior density of null model fit
  d.pred$pred.sample=apply(post,2,mean)#Take mean of posterior samples of canopy model parameter estimates (should be ~ the same as predicted values above)
  d.pred$sd.sample=apply(post,2,sd)#Generate a standard deviation of posterior samples of parameter estimates.
  global.mean=c(coef(m0)[[1]],apply(post0,2,sd)[[1]])#Find mean and standard deviation of population mean parameter estimate (from null model). This global.mean value is used when there are no significant treatment effects at any level.
  return(list(daicc=daicc,summary=summary(m1),d.pred=d.pred,global.mean=global.mean))
}
FireEffects=function(P,md){
  md=md[md$Canopy=="T",] #Fire effects only relevant in thinned-canopy units.
  #glm is the proper model here because the variation in fire is *at* the block level
  md.glm=ddply(md,.(Block),summarize,Fire=unique(Fire),y=sum(y),n=length(Block),p=y/n)
  m0=glm(p~1, data=md.glm, weights=n, family="binomial")
  m1=glm(p~Fire, data=md.glm, weights=n, family="binomial")
  daicc=AICc(m0)-AICc(m1)
  d.pred=data.frame(Fire=c(0,1))
  d.pred$pred=predict(m1,d.pred)
  post=sample.naive.posterior(m1,n=10000)
  post[,2]=post[,1]+post[,2] #This line was not in the code for the MS; would reduce the standard deviations for the bold values in Table 2, column T-B very slightly (not affecting any results).
  post0=sample.naive.posterior(m0,n=10000)
  d.pred$pred.sample=apply(post,2,mean)
  d.pred$sd.sample=apply(post,2,sd)
  global.mean=c(coef(m0)[[1]],apply(post0,2,sd)[[1]])
  return(list(daicc=daicc,summary=summary(m1),d.pred=d.pred,global.mean=global.mean))
}
SnowEffects=function(P,md,subset="",fam="binomial"){
  C=subset[1]
  if(!""%in%subset){ #Subset the data if needed (i.e. into a specific canopy group). Should always run.
    md=md[which(md$Canopy==subset[which(subset%in%md$Canopy)] &
        md$Fire%in%subset[which(subset%in%md$Fire)]  ) ,]
  }
  #Run models
  if(length(unique(md$Block))>1 & var(md$y)>0 & #If there are multiple blocks and variance in the data, and if it is not one of the special cases below, run mixed model (glmer)
  !(spp=="sj"&C=="O"&(P=="Pg2"|P=="Pws1")) &
  !(spp=="sj"&C=="D"&(P=="PgB"|P=="Pss1")) 
  #Special Case: The above exemptions !() skip variables that crash glmer (model does not converge).  These models are run as simple glm's. 
  ){
    m0=glmer(y~1 + (1|Block),data=md,family=fam)
    m1=glmer(y~SnowTot + (1|Block),data=md,family=fam)
    #m2=glmer(y~Snow + (1|Block),data=md,family=fam) #This is a deprecated model that treated "Snow" as a categorical variable (increased/ambient/decreaesd); in the final version we treated snow as a continuous variable (SnowTot= cumulative winter snowfall).
    #For centimeters, c(10,40,70)=c(25.4, 101.6, 177.8)
    d.pred=data.frame(SnowTot=c(25,100,175))
    #Predict using snow model
    post=ezPredict(m1,d.pred,iterations=10000) #returns predicted value and variance of that value
    d.pred=post$cells
    names(d.pred)=c("SnowTot","pred.sample","var.sample")
    d.pred$sd.sample=sqrt(d.pred$var.sample)
    global.mean=fixef(m0)[[1]]
    daicc=c(AICc(m0)-AICc(m1))
  } else if(var(md$y)==0 | is.na(var(md$y))){ #If there's no variance in the response
    print("no variance")
    d.pred=data.frame(SnowTot=c(25,100,175))
    d.pred[,c("pred.sample")]=ifelse(mean(md$y)==0,-20,20) #prob is either low or high (inv.logit of -20 or 20) 
    d.pred[,"sd.sample"]=d.pred[,"var.sample"]=0
    return(list(daicc=0,d.pred=d.pred))
  } else { #If there's no blocking structure, use a glm
    m0=glm(y~1,data=md,family=fam)
    m1=glm(y~SnowTot,data=md,family=fam)
    d.pred=data.frame(SnowTot=c(25,100,175))
    global.mean=coef(m0)[[1]]
    daicc=c(AICc(m0)-AICc(m1))
    print("glm for the model below (not glmer)")
  }
  ##Below are four "Special Cases" where small sample size and/or poor model convergence leads to unreliable results and requires overriding/setting daicc for the snow effect.
  if( 
      (spp=="sj" & C=="O" & P=="Pg2") | #Special Case: Only three germinants (of 275 seeds) during the second growing season, all in the same block and same snow treatment. Mixed model does not converge, and glm has daicc>2 because all 3 germinants were in the deepest snow plot, so parameter estimates are very large but unreliable. Thus we treat it as an outlier and remove the snow effect.
      (spp=="sj" & C=="D" & P=="Pss1") | #Special Case: Very small sample size (n=14 seedlings alive in their first spring in dense shade, only 2 survived [14.2%]). Mixed model does converge; parameter estimate for snow was very large because oen of the two surviving seedlings had its preceeding winter in the dry winter of 2012-2013, while the rest of the seedlings had their preceeding winter in the wetter winter of 2011-2012, so this individual heavily weights the parameter estimate towards increased survival under reduced snow. Thus we treat it as an outlier and remove the snow effect.
      (spp=="cs" & C=="O" & P=="Pws2") |#Special Case: Very small sample size (n=10); mixed model converges (barely) and daicc >2, but parameter estimates are extremely large because of the small sample size (the only two brooms that died in their 2nd winter experienced more snow relative to those that survived, heavily biased by one broom that germinated in the first fall, so its 2nd winter was the wetter winter of 2012-13). Thus we treat it as an outlier and remove the snow effect.
      (spp=="cs" & C=="T" & P=="Pm3") #Special Case: in both unburned and burned treatments, there was only one block each that had any mature individuals, and the sample size was low. For the unburned treatment, the only two 3rd-year individuals that were alive in the reduced-snow plot both flowered, while in the burned treatment, the only 3rd-year individuals that flowered were both in the increased-snow plot. These coefficients have a large-magnitude effect size because of the small sample sizes, although the model dAICC's were small for both unburned and burned treatments, and the snow effect coefficients are non-significant. Thus we treat it as an outlier and remove the snow effect.
  ){daicc=0}
      
  
  print(c(P,C))
  print(paste("daicc=",daicc))
  #print(summary(m1))
  print(d.pred)
  return(list(daicc=daicc,summary=summary(m1),d.pred=d.pred,global.mean=global.mean))
}  
CanopyEffectsFecundity=function(P,md){#Fecundity statistical model requires poisson
  #glm is the proper model here because the variation in canopy is *at* the block level
  m0=glm(y~1, data=md, family="poisson")
  m1=glm(y~Canopy, data=md, family="poisson")
  m1b=glm(y~0+Canopy, data=md, family="poisson") #Coding the model this way (0+) makes predictions easier; this is the same model as m1.
  daicc=AICc(m0)-AICc(m1)
  d.pred=data.frame(Canopy=c("O","T"))
  d.pred$pred=predict(m1,d.pred)
  post=sample.naive.posterior(m1b,n=10000)
  post0=sample.naive.posterior(m0,n=10000)
  d.pred$pred.sample=apply(post,2,mean)
  d.pred$sd.sample=apply(post,2,sd)
  global.mean=c(coef(m0)[[1]],apply(post0,2,sd)[[1]])
  return(list(daicc=daicc,summary=summary(m1),d.pred=d.pred,global.mean=global.mean))
}
FireEffectsFecundity=function(P,md){#Fecundity statistical model requires poisson
  md=md[md$Block<9,]
  m0=glm(y~1, data=md, family="poisson")
  m1=glm(y~Fire, data=md, family="poisson")
  m1b=glm(y~0+Fire, data=md, family="poisson") 
  daicc=AICc(m0)-AICc(m1)
  d.pred=data.frame(Fire=c(0,1))
  d.pred$pred=predict(m1,d.pred)
  post=sample.naive.posterior(m1,n=10000)
  post[,2]=post[,1]+post[,2] #This line was not in the code for the MS, but does not affect any results.
  post0=sample.naive.posterior(m0,n=10000)
  d.pred$pred.sample=apply(post,2,mean)
  d.pred$sd.sample=apply(post,2,sd)
  global.mean=c(coef(m0)[[1]],apply(post0,2,sd)[[1]])
  return(list(daicc=daicc,summary=summary(m1),d.pred=d.pred,global.mean=global.mean))
}


####3b. Function Calls####
#Special calls first, where Scotch and Spanish differ and no snow and/or fire effects
####3b.1 Special Case Function Calls (Pg0 and S)####
#Special case for initial germination, which has no snow effect
P="Pg0"
md=GenerateData(P=P,d=d)
m.out.l[["Pg0.CC"]]=CanopyEffects(P=P,md=md,subset="")
m.out.l[["Pg0.F"]]=FireEffects(P=P,md=md)
#Special case for fecundity, which requires a poisson model and has no snow/fire effects
if(spp=="cs"){#Only have fecundity data for Scotch broom
P="S1"
md=GenerateData(P=P,d=d)
m.out.l[["S1.CC"]]=CanopyEffectsFecundity(P=P,md=md)
#md$Fire[c(3,5)]=1; md$Block[c(3,5)]=7 #** #Deprecated?
md[md$ID==1103,"Snow"]="A"; md[md$ID==1103,"SnowTot"]=25.40 #**
m.out.l[["S1.F"]]=FireEffectsFecundity(P=P,md=md)
md$y[10]=md$y[12];md$y[12]=0 #**
md=md[-c(11),] #** remove outliers with huge seed set
m.out.l[["S1.Snow.O"]]=SnowEffects(P=P,md=md,subset=c("O",0),fam="poisson")
m.out.l[["S1.Snow.T.UB"]]=SnowEffects(P=P,md=md,subset=c("T",0),fam="poisson")
#m.out.l[["S1.Snow.T.B"]]=SnowEffects(P=P,md=md,subset=c("T",1),fam="poisson")
#This one doesn't work because there is only one plot and therefore only one snow level represented in the burned plot. Make NA in table for paper.
}#End Scotch broom-only fecundity section.

#Fill out parameter table
SpecialPars <- c("Pg0","S1")
if(spp=="sj"){SpecialPars = "Pg0"}

for(P in SpecialPars){
if(m.out.l[[paste(P,".CC",sep="")]]$daicc>2){#If canopy affects parameters
  d.pred=m.out.l[[paste(P,".CC",sep="")]]$d.pred
  for(C in as.character(d.pred$Canopy)){
    bst[grep(C,bst$Scenario),grep(P,names(bst))]=d.pred[d.pred$Canopy==C,c(2,4)]
  }
  }else {
    print("no canopy effect")
    bst[,grep(P,names(bst))]=
      data.frame(matrix(m.out.l[[paste(P,".CC",sep="")]]$global.mean,nrow=1))
  }
if(m.out.l[[paste(P,".F",sep="")]]$daicc>2){#If fire affects parameters (Pg0 only)
  d.pred=m.out.l[[paste(P,".F",sep="")]]$d.pred
  bst[c(4:6),grep(P,names(bst))]=d.pred[1,c(2,4)] #Unburned
  bst[c(7:9),grep(P,names(bst))]=d.pred[2,c(2,4)] #Burned
}#End parameter extraction from fire model
}#End loop for filling parameter table with special parameters

####3b.2 Function Calls (standard binomial probabilities)####
Probs <- c("Pg1","Pg2","PgB","Pss1","Pws1","Pss2","Pws2","Pm3")
if(spp=="sj"){Probs = c("Pg1","Pg2","PgB","Pws0","Pss1","Pws1","Pss2","Pws2","Pm3")} 

for(P in Probs){# Takes a while
md=GenerateData(P=P,d=d)
m.out.l[[paste(P,".CC",sep="")]]=CanopyEffects(P=P,md=md,subset="")
if(m.out.l[[paste(P,".CC",sep="")]]$daicc>2){#If canopy affects parameters
  d.pred=m.out.l[[paste(P,".CC",sep="")]]$d.pred
  for(C in as.character(d.pred$Canopy)){
  bst[grep(C,bst$Scenario),grep(P,names(bst))]=d.pred[d.pred$Canopy==C,c(2,4)]
  }}else {
    print("no canopy effect for transition probability listed below")
    bst[,grep(P,names(bst))]=
      data.frame(matrix(m.out.l[[paste(P,".CC",sep="")]]$global.mean,nrow=1))
  }
m.out.l[[paste(P,".F",sep="")]]=FireEffects(P=P,md=md)
if(m.out.l[[paste(P,".F",sep="")]]$daicc>2){#If fire affects parameters
  d.pred=m.out.l[[paste(P,".F",sep="")]]$d.pred
  bst[c(4:6),grep(P,names(bst))]=d.pred[1,c(2,4)] #Burned
  bst[c(7:9),grep(P,names(bst))]=d.pred[2,c(2,4)] #Unburned
  }
for(C in Canopy){ #C="O"
  m.out=SnowEffects(P=P,md=md,subset=c(C,0))
  if(m.out$daicc[1]>2){#Override canopy effect with snow effect
    bst[grep(C,bst$Scenario)[1:3],grep(P,names(bst))]= #[1:3] skips burned plots, handled later
      m.out$d.pred[,c("pred.sample","sd.sample")] }
  if(C=="T"){#If plot was thinned, also need to extract parameters for fire
  m.out=SnowEffects(P=P,md=md,subset=c(C,1))
  if(m.out$daicc[1]>2){#Override canopy effect with snow effect
    bst[grep("T.B",bst$Scenario),grep(P,names(bst))]=m.out$d.pred[,c("pred.sample","sd.sample")]}
  }
}
} #End of statistical model loop.

####4. Tidy and print####
#Finally, edit parameters where there were 0 observed successes (setting the logit value to -20):
if(spp=="cs"){ 
  bst[c(13:15),grep("Pg0",names(bst))]=as.data.frame(matrix(c(-20,0),nrow=1))  
  bst[c(10:15),grep("S1",names(bst))]=as.data.frame(matrix(c(-20,0),nrow=1))
}
if(spp=="sj"){ 
  cs.pars=read.csv("data/csParamTableBST.csv")
 bst[,c("S1","S1sd")]=cs.pars[,c("S1","S1sd")] #Assign Spanish broom the 1st-year seed production data that Scotch broom had.
}

bst[,c(3:ncol(bst))]=round(bst[,c(3:ncol(bst))],3)
bst[,c(3:ncol(bst))][bst[,c(3:ncol(bst))]>100]=0 #Responses with 0 probability get SD=0
#write.csv(bst,file=paste(spp,"dataParamTableBST.csv",sep=""),row.names=F) #Option to write to file. CAUTION this overwrites existing file