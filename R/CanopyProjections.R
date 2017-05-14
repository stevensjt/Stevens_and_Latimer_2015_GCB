#Code for Stevens and Latimer 2015
#Jens Stevens stevensjt@gmail.com 11/21/14
##This code generates transition matrices from many samples of parameter distributions, which it reads in from csParamTableBST.csv and sjParamTableBST.csv. Need to manually set spp to "cs" (Scotch broom) and then "sj" (Spanish broom) and run #s 1 and 2, in that order.

library(boot) #for inv.logit(). version 1.3-18
library(popbio) #for lambda(). version 2.4.3
library(reshape) #for melt(). version 0.8.6
library(ggplot2) #for ggplot(). version 2.2.1
library(gridExtra) #for grid.arrange(). version 2.2.1
library(grid) #for textGrob(). version 3.3.2

spp="sj"

####1. Data Setup####
ParameterList=c("Pg0","Pg1","Pg2","PgB","Pd1","Pd2","PdB",
                "Pws0","Pss1","Pws1","Pss2","Pws2","Py1s0","Py1s1","Py1sB",
                "Py2","Pspa","Psa","Pm3","Ppa","Pmpa","Pa","S1" )
Scenarios <- c("O.S10","O.S40","O.S70",
               "T.UB.S10","T.UB.S40","T.UB.S70",
               "T.B.S10","T.B.S40","T.B.S70",
               "C.S10","C.S40","C.S70",
               "D.S10","D.S40","D.S70") #Canopy.Burned/UnBurned/Snow (inches). For the MS, ran rougly approximate simulations in cm (25,100,175), but kept the 10-40-70 convention in the code.

#Initialize matrices
#pm is parameter matrix, one row for each parameter, one column for each scenario. Parameters will be on logit scale.
pm <- data.frame(row.names=ParameterList) 
pm[,c(1:length(Scenarios))]  <- NA
names(pm) <- Scenarios
#sdm is a matrix of standard deviations that correspond to the parameters, on logit scale.
sdm <- pm
#pred.pm is a matrix of predicted parameters, the inverse logit of a sample from the parameter distribution with mean drawn from pm and standard deviation drawn from sdm
pred.pm <- pm

#Fill in parameter matrix and standard deviation matrix
if(spp=="cs"){ #**Run "cs" first, "sj" needs maturity parameters from "cs"
  ParamTable <- read.csv("data/csParamTableBST.csv")
  csParamHolder <- ParamTable} else {
    ParamTable <- read.csv("data/sjParamTableBST.csv")
    ParamTable[,grep("Pm3",names(ParamTable))] <- 
      csParamHolder[,grep("Pm3",names(ParamTable))]
}    

for(i in ParameterList){ #Fill out parameters
  if(length(grep(i,names(ParamTable)[c(3:ncol(ParamTable))]))>0){ #if parameter was modeled:
  mui=grep(i,names(ParamTable))[1] #set the column index for mu Beta
  sdi=grep(i,names(ParamTable))[2] #set the column index for standard deviation
  pm[i,]=as.numeric(as.character(ParamTable[,mui])) #get parameter for Scenario i, on logit scale
  sdm[i,]=as.numeric(as.character(ParamTable[,sdi])) #get stdev for Scenario i
  }
}
#The list pred.pm.l will hold j different predicted parameter matrices for all Scenarios.
pred.pm.l=list() 

####2. Sample parameter distributions to generate parameter tables####
for(j in c(1:1000)){#Fill out many parameter tables; will take about 1:15 to run 1000.

####2a. Sample parameter distributions j times, filling out the predicted parameter matrix####
for(i in ParameterList){ #Sample modeled parameters from distribution
  #might eventually want to speed up w/ apply
  if(length(grep(i,names(ParamTable)[c(3:ncol(ParamTable))]))>0){ #if parameter was modeled:
  pred.pm[i,]=round(inv.logit(rnorm(n=15,mean=as.numeric(pm[i,]),sd=as.numeric(sdm[i,]))),4)
  if(i=="S1"){#Treat seeds differently because parameter comes from Poisson process
  pred.pm[i,]=round(exp(rnorm(n=15,mean=as.numeric(pm[i,]),sd=as.numeric(sdm[i,]))),4)}
  }
}

####2b. Fill in derived & a-priori probabilities####
pred.pm["Pd1",]=(1-pred.pm["Pg0",]-pred.pm["Pg1",]) #Prob dormancy year 1
pred.pm["Pd2",]=(1-pred.pm["Pg2",]) #Prob dormancy year 2
#PdB: if there had been a fire, PdB is the ratio of seed bank germination rate in Sp3 to germination rate in Sp2, otherwise it is 68% (published estimate for Scotch broom from Haubensak et al 2004).
pred.pm["PdB",]=0.680
pred.pm["PdB",c(4:6)]= #columns 4-6 being the burned columns, overwrite the published estimate.
  round(max(pred.pm["PgB",c(4:6)]/pred.pm["Pg2",c(4:6)]),4)
#Py1s0 (prob of going from new seed to seedling after 1 year) depends on species
#For scotch, small initial cohort survives winter with Pws1, not enough data to estimate Pws0, so set Pws0=Pws1
if(spp=="cs"){pred.pm["Pws0",]=pred.pm["Pws1",]} #Set Scotch Pws0 equal to Pws1
pred.pm["Py1s0",]=round((pred.pm["Pg0",]*pred.pm["Pws1",]*pred.pm["Pss1",])+(pred.pm["Pg1",]*pred.pm["Pss1",]),4)
#For spanish, larger initial cohort gets its own Pws0
if(spp=="sj"){ pred.pm["Py1s0",]=round((pred.pm["Pg0",]*pred.pm["Pws0",]*pred.pm["Pss1",])+(pred.pm["Pg1",]*pred.pm["Pss1",]),4) 
}
pred.pm["Py1s1",]=round((pred.pm["Pg2",]*pred.pm["Pss1",]),4)
pred.pm["Py1sB",]=round((pred.pm["PgB",]*pred.pm["Pss1",]),4)
pred.pm["Py2",]=round((pred.pm["Pws1",]*pred.pm["Pss2",]),4)
#A-priori probability of going straight to reproductive adult from Yr2:
pred.pm["Pspa",]=round(pred.pm["Pws2",]*pred.pm["Pss2",]*(1-pred.pm["Pm3",]),4)#Prob of transitioning from a seedling in year 2 to a pre-adult
if(spp=="sj"){pred.pm["Pspa",]=round(pred.pm["Pws2",]*pred.pm["Pss2",])} #No loss to maturity for sj
pred.pm["Psa",]=round(pred.pm["Pws2",]*pred.pm["Pss2",]*pred.pm["Pm3",],4)#Prob of transitioning from a seedling in year 2 to an adult
if(spp=="sj"){pred.pm["Psa",]=0} #no 3rd year maturity for sj
pred.pm["Ppa",]=0.8 #A-priori prob of pre-adults surviving (this is in the ballpark of highest survival from Yr1 to Yr2)
pred.pm["Pmpa",]=pred.pm["Pm3",]
pred.pm["Pmpa",c(10:15)]=pred.pm["Pmpa",c(4:9)]#Set Closed and Dense = Thinned (lowest Pm3 value)
pred.pm["Pa",]=0.8 #A-priori prob of an adult surviving.
pred.pm["S1",c(10:15)]=pred.pm["S1",c(4:9)]#Set Closed and Dense = Thinned (lowest Pm3 value)

#transfer this predicted parameter matrix into the list at element j:
pred.pm.l[[j]]=pred.pm
}
#End of the Sample Parameter Distributions for-loop is the end of 2b

####2c. Calculate the mean predicted matrix for each scenario ####
temp.mat.l=lapply(pred.pm.l,as.matrix)
if(spp=="cs"){
    cs.pm.bst=apply(simplify2array(temp.mat.l), c(1,2), mean)} else {
    sj.pm.bst=apply(simplify2array(temp.mat.l), c(1,2), mean)}
#Export mean predicted matrices
#**If writing to file (below), run the above twice, once for "cs" and once for "sj", in that order.**#
#write.csv(cs.pm.bst,file="data/cs.pm.bst.csv",row.names=T) #Write best projection matrix to file (optional)
#write.csv(sj.pm.bst,file="data/sj.pm.bst.csv",row.names=T) #Write best projection matrix to file (optional)
#Once the pm.bst files are generated, don't need to write to file anymore. They are only used for the Sensitivity and Elasticity analyses and can be imported from the saved files.

####3a. Make a transition matrix for each simulation####
AnnMatStages=c("Seed1","SeedB","Yr1","Yr2","PreA","A")
Scenarios=c("O.S10","O.S40","O.S70","T.UB.S10","T.UB.S40","T.UB.S70","C.S10","C.S40","C.S70","D.S10","D.S40","D.S70")#Revise the scenarios to be unburned only after reading in the data.
#tm.template is the empty transition matrix which will get filled in. 
tm.template=data.frame(row.names=AnnMatStages);tm.template[,c(1:length(AnnMatStages))]=NA
names(tm.template)=AnnMatStages
tm.l=list();pp=list();lambdas=list();scenarios=list();
lambda.mat=matrix(nrow=length(pred.pm.l),ncol=length(Scenarios),dimnames=list(NULL,c(Scenarios))  )
ci=1
#ci=which(Scenarios=="O.S40") #To select scenario manually
for(j in c(1:length(pred.pm.l))){#Will take 1:09 with 1000 samples of j
for(ci in c(1:length(Scenarios))){ 
tm=tm.template
c=Scenarios[ci]
#Fill out transition matrix!
tm["SeedB","Seed1"]=pred.pm.l[[j]]["Pd2",c] #Seed from initial dormancy to seed bank
tm["Yr1","Seed1"]=pred.pm.l[[j]]["Py1s1",c] #Seed from initial dormancy to seedling
tm["SeedB","SeedB"]=pred.pm.l[[j]]["PdB",c] #Seed from bank staying viable
tm["Yr1","SeedB"]=pred.pm.l[[j]]["Py1sB",c] #Seed from bank to seedling
tm["Yr2","Yr1"]=pred.pm.l[[j]]["Py2",c] #Seedling in yr 1 to seedling in yr 2
tm["PreA","Yr2"]=pred.pm.l[[j]]["Pspa",c] #Seedling in yr 2 to a pre-adult
tm["A","Yr2"]=pred.pm.l[[j]]["Psa",c] #Seedling in yr 2 to an adult
tm["PreA","PreA"]=pred.pm.l[[j]]["Ppa",c] #Staying as a pre-adult
tm["A","PreA"]=pred.pm.l[[j]]["Pmpa",c] #Pre-adult to adult
tm["Seed1","A"]=#Adult to dormant seed = prob. of initial dormancy * N seeds produced
  pred.pm.l[[j]]["Pd1",c]*pred.pm.l[[j]]["S1",c] 
tm["Yr1","A"]=#Adult to seedling = prob of seed becoming a 1st year sdlg (Py1s0) * N seeds produced
  pred.pm.l[[j]]["Py1s0",c]*pred.pm.l[[j]]["S1",c] 
tm["A","A"]=pred.pm.l[[j]]["Pa",c] #Staying as an adult
tm[is.na(tm)]=0
tm.l[[j]]=as.matrix(tm)
lambda.mat[j,c]=lambda(tm.l[[j]])
}
}

####3b: Sensitivity and elasticity analyses (can skip if just making plots)####
ScenarioNamesMetric=c("O.25","O.100","O.175",
                      "T.UB.25","T.UB.100","T.UB.175",
                      "T.B.25","T.B.100","T.B.175",
                      "C.25","C.100","C.175",
                      "D.25","D.100","D.175")
matrix.names=c("cs.bst","sj.bst")
pred.pm.l=list();tm.l=list();tm.all.l=list();sens=list();elas=list()
pred.pm.l[["sj.bst"]]=read.csv("data/sj.pm.bst.csv",header=T,row.names=1)#Read in the average parameters
pred.pm.l[["cs.bst"]]=read.csv("data/cs.pm.bst.csv",header=T,row.names=1)
for(j in matrix.names){#This will print sensitivity and elasticity scores for all selected matrices
  spp.name=ifelse(j=="cs.bst","Scotch", "Spanish")
  for(ci in c(1:length(Scenarios))){ 
  tm=tm.template
  c=Scenarios[ci];cm=ScenarioNamesMetric[ci]
  #Fill out transition matrix!
  tm["SeedB","Seed1"]=pred.pm.l[[j]]["Pd2",c] #Seed from initial dormancy to seed bank
  tm["Yr1","Seed1"]=pred.pm.l[[j]]["Py1s1",c] #Seed from initial dormancy to seedling
  tm["SeedB","SeedB"]=pred.pm.l[[j]]["PdB",c] #Seed from bank staying viable
  tm["Yr1","SeedB"]=pred.pm.l[[j]]["Py1sB",c] #Seed from bank to seedling
  tm["Yr2","Yr1"]=pred.pm.l[[j]]["Py2",c] #Seedling in yr 1 to seedling in yr 2
  tm["PreA","Yr2"]=pred.pm.l[[j]]["Pspa",c] #Seedling in yr 2 to a pre-adult
  tm["A","Yr2"]=pred.pm.l[[j]]["Psa",c] #Seedling in yr 2 to an adult
  tm["PreA","PreA"]=pred.pm.l[[j]]["Ppa",c] #Staying as a pre-adult
  tm["A","PreA"]=pred.pm.l[[j]]["Pmpa",c] #Pre-adult to adult
  tm["Seed1","A"]=#Adult to dormant seed = prob. of initial dormancy * N seeds produced
    pred.pm.l[[j]]["Pd1",c]*pred.pm.l[[j]]["S1",c] 
  tm["Yr1","A"]=#Adult to seedling = prob of seed becoming a 1st year sdlg (Py1s0) * N seeds produced
    pred.pm.l[[j]]["Py1s0",c]*pred.pm.l[[j]]["S1",c] 
  tm["A","A"]=pred.pm.l[[j]]["Pa",c] #Staying as an adult
  tm[is.na(tm)]=0
  tm.l[[c]]=as.matrix(tm)
  #pp[[j]]=pop.projection(tm.l[[j]],n=c(rep(0,times=5),10),iterations=100) #Start with 10 adults
  #lambda.mat[j,c]=lambda(tm.l[[j]])
  sens[[c]]=sensitivity(tm.l[[c]])
  sens[[c]][tm.l[[c]]==0]=0 #Set unrealized transitions to 0 (justified in Caswell book pp 215-217)
  #image2(sens[[c]], mar=c(1,4,5,1), box.offset=.1) #Optional: Run to view sensitivity matrices
  #title(paste("Sensitivity",spp.name,cm), line=2.5) #Optional: Run to view sensitivity matrices
  #dev.copy2pdf(file=paste("Sensitivity",j,c,"pdf",sep=".")) #Optional: Write to file
  elas[[c]]=elasticity(tm.l[[c]])
  #image2(elas[[c]], mar=c(1,4,5,1), box.offset=.1) #Optional: Run to view elasticity matrices
  #title(paste("Elasticity",spp.name,cm), line=2.5) #Optional: Run to view elasticity matrices
  #dev.copy2pdf(file=paste("Elasticity",j,c,"pdf",sep="."))  #Optional: Write to file
}
tm.all.l[[j]]=tm.l
}

####4. Make some plots####
####4a. Set up data for visualization####
lambdas.gg=melt(lambda.mat)
names(lambdas.gg)=c("Sim","scenario","lambda")
lambdas.gg[grep("O.",lambdas.gg$scenario),"Canopy"]="Open"
lambdas.gg[grep("T.",lambdas.gg$scenario),"Canopy"]="Thinned"
lambdas.gg[grep("C.",lambdas.gg$scenario),"Canopy"]="Closed"
lambdas.gg[grep("D.",lambdas.gg$scenario),"Canopy"]="Dense"
lambdas.gg$Canopy=factor(lambdas.gg$Canopy,levels=c("Open","Thinned","Closed","Dense"))
lambdas.gg[grep("S10",lambdas.gg$scenario),"Snow"]="25 cm"
lambdas.gg[grep("S40",lambdas.gg$scenario),"Snow"]="100 cm"
lambdas.gg[grep("S70",lambdas.gg$scenario),"Snow"]="175 cm"
lambdas.gg$Snow=factor(lambdas.gg$Snow,levels=c("25 cm","100 cm", "175 cm"))

if(spp=="cs"){scotch.lambdas=lambdas.gg
  
  }else{spanish.lambdas=lambdas.gg
  spanish.lambdas[which(spanish.lambdas$scenario=="C.S10"&spanish.lambdas$lambda<1),"lambda"]=
    mean(spanish.lambdas$lambda)} #Gets rid of one outlier

##STOP: Run lines 1-212 for spp="cs" and then again for spp="sj"##
####4b. Make and arrange plots####
p.cs <- ggplot(scotch.lambdas)
p.sj <- ggplot(spanish.lambdas)
p.scotch <- p.cs+
  geom_line(stat="density",aes(x=lambda,y=..scaled..,color=Snow),trim=T,show.legend=F)+
  #ylim(c(0,5))+
  xlim(c(0.5,2.6))+
  facet_wrap(~Canopy,ncol=1)+labs(title="Scotch Broom", y="scaled density")+
  geom_vline(xintercept=1)+
  theme_bw()+   theme(plot.title = element_text(hjust = 0.5))
p.spanish <- p.sj+
  geom_line(stat="density",aes(x=lambda,y=..scaled.., color=Snow),trim=T)+
  #ylim(c(0,5))+
  xlim(c(0.5,2.6))+
  facet_wrap(~Canopy,ncol=1)+labs(title="Spanish Broom", y="scaled density")+
  geom_vline(xintercept=1)+  
  theme_bw()+  theme(plot.title = element_text(hjust = 0.5))

p.scotch;p.spanish  

a2 <- ggplot_gtable(ggplot_build(p.scotch))
a2$layout$clip[a2$layout$name == "panel"] <- "off"
b2 <- ggplot_gtable(ggplot_build(p.spanish))
b2$layout$clip[b2$layout$name == "panel"] <- "off"
grid.arrange(arrangeGrob(a2,b2,
            top = textGrob("(b)", x=0.42, hjust = -0.1),
            left = textGrob("(a)", x=0, y=0.98, hjust = -1.7),
            padding = unit(0.1, "line"),
            ncol=2,widths=c(0.4,0.6)),nrow=1 )
#For the final version, I moved the legend inside one of the faceted plots, and changed the lines to be dashed lines, but I can't find the code I used to do that. It's possible it was done in Illustrator.
#7x8 is a good dimension for the paper
#dev.copy2eps(file="Figures/Fig3Export.eps", width=7, height=8,fonts=c("Helvetica","sans"))