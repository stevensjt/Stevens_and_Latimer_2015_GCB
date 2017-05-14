#Code for Stevens and Latimer 2015
#Jens Stevens stevensjt@gmail.com 11/21/14
##This code generates transition matrices from many samples of parameter distributions, which it reads in from csParamTableBST.csv and sjParamTableBST.csv (specific to fire simulations). This code calculates lambda (with confidence intervals) for fire return intervals from 1-30.

require(boot) #for inv.logit(). version 1.3-18
require(popbio) #for lambda(). version 2.4.3
require(reshape) #for melt(). version 0.8.6
require(ggplot2) #for ggplot(). version 2.2.1
nthroot=function(x,n){x^(1/n)}
matexp<-function(A,n){
  if(n==1) A else {B<-A; for(i in (2:n)){A<-A%*%B}}; A
} 

#NEED TO RUN CANOPYPROJECTIONS.R
#Set species, and then run sections 1+2
####1. Data Setup: Run from CanopyProjections.R####
####2. Sample parameter distributions: Run from CanopyProjections.R####
#(This generates parameter tables)

####3. Make a transition matrix for each fireXsnow scenario####
####This is where we start to play around with fire
####This part has been modified for fire-effects only
ScenariosF=c("T.B.S10","T.B.S40","T.B.S70","T.UB.S10","T.UB.S40","T.UB.S70")
AnnMatStages=c("Seed1","SeedB","Yr1","Yr2","PreA","A")
#tm.template is the empty transition matrix which will get filled in. 
tm.template=data.frame(row.names=AnnMatStages);tm.template[,c(1:length(AnnMatStages))]=NA
names(tm.template)=AnnMatStages
tm.b.l=tm.ub.l=list();pp=list();lambdas=list();scenarios=list()
fri.mat=data.frame(row.names=c(1:30))
ci=1

for(j in c(1:length(pred.pm.l))){
for(ci in c(1:3)){#For each snow scenario  
tm.b=tm.ub=tm.template
cb=ScenariosF[ci];cub=ScenariosF[ci+3]
##Fill out burned transition matrix
tm.b["Yr1","Seed1"]=pred.pm.l[[j]]["Py1s1",cb]
tm.b["SeedB","SeedB"]=pred.pm.l[[j]]["PdB",cb]
tm.b["Yr1","SeedB"]=pred.pm.l[[j]]["Py1sB",cb]
tm.b["Yr2","Yr1"]=0 #Complete seedling mortality (Agee 96)
tm.b["PreA","Yr2"]=0 #Complete seedling mortality
tm.b["A","Yr2"]=0 #Complete seedling mortality
tm.b["PreA","PreA"]=0.13 #13% survival of pre-adults (Tveten 99)
tm.b["A","PreA"]=0 #No pre-adults transition to adults
tm.b["Seed1","A"]=pred.pm.l[[j]]["Pd1",cb]*pred.pm.l[[j]]["S1",cb]
tm.b["Yr1","A"]=pred.pm.l[[j]]["Py1s0",cb]*pred.pm.l[[j]]["S1",cb]
tm.b["A","A"]=0.13 #13% survival of pre-adults (Tveten 99)
tm.b[is.na(tm.b)]=0
tm.b.l[[j]]=as.matrix(tm.b)
#View(tm.b.l[[j]])
##Fill out unburned transition matrix
tm.ub["Yr1","Seed1"]=pred.pm.l[[j]]["Py1s1",cub]
tm.ub["SeedB","SeedB"]=pred.pm.l[[j]]["PdB",cub]
tm.ub["Yr1","SeedB"]=pred.pm.l[[j]]["Py1sB",cub]
tm.ub["Yr2","Yr1"]=pred.pm.l[[j]]["Py2",cub]
tm.ub["PreA","Yr2"]=pred.pm.l[[j]]["Pspa",cub]
tm.ub["A","Yr2"]=pred.pm.l[[j]]["Psa",cub]
tm.ub["PreA","PreA"]=pred.pm.l[[j]]["Ppa",cub]
tm.ub["A","PreA"]=pred.pm.l[[j]]["Pmpa",cub]
tm.ub["Seed1","A"]=pred.pm.l[[j]]["Pd1",cub]*pred.pm.l[[j]]["S1",cub]
tm.ub["Yr1","A"]=pred.pm.l[[j]]["Py1s0",cub]*pred.pm.l[[j]]["S1",cub]
tm.ub["A","A"]=pred.pm.l[[j]]["Pa",cub]
tm.ub[is.na(tm.ub)]=0
tm.ub.l[[j]]=as.matrix(tm.ub)
#View(tm.ub.l[[j]])

fri.mat[1,cb]=lambda(tm.b.l[[j]]) #Despite differences in matrices between scenarios, lambda always end up the same here, because of the seedling bottleneck (populations die in the same amount of time despite more establishment with less snow)
for(n in 2:30){ #Calculate lambda for range of fire return intervals
  mat.tmp=matexp(tm.ub.l[[j]],n-1)%*%tm.b.l[[j]]
  fri.mat[n,cb]=nthroot(x=as.double(eigen(mat.tmp)$values[1]),n=n)
}
}
lambdas[[j]]=as.matrix(fri.mat)
scenarios[[j]]=cb
}


####4. Make some plots####
lambdas.mu=apply(simplify2array(lambdas), c(1,2), mean)
lambdas.sd=apply(simplify2array(lambdas), c(1,2), sd)
lambdas.mu.gg=melt(lambdas.mu)
names(lambdas.mu.gg)=c("FRI","Snowfall","lambda")
lambdas.mu.gg$sd=melt(lambdas.sd)$value
levels(lambdas.mu.gg$Snowfall)=c("25 cm","100 cm", "175 cm")
if(spp=="cs"){lambdas.mu.gg.cs=lambdas.mu.gg}else{lambdas.mu.gg.sj=lambdas.mu.gg}#Storage
p=ggplot(lambdas.mu.gg, aes(x=FRI))
p2=p+
  geom_ribbon(aes(ymin=lambda-sd,ymax=lambda+sd,fill=Snowfall),alpha=0.2)+
  geom_line(aes(y=lambda,color=Snowfall))+
  ylab(expression(lambda))+  geom_hline(yintercept=1)+  
  theme_bw()+ theme(plot.title = element_text(hjust = 0.5))

#Customize plots for species
if(spp=="cs"){p2.cs=p2+ggtitle("Scotch broom")+ 
  ylim(c(0.2,1.8))+ theme(legend.position="none")+
  annotate("text",label="(a)",x=Inf,y=Inf,hjust=0.8,vjust=-0.5)}else
{p2.sj=p2+ggtitle("Spanish broom")+annotate("text",label="(b)",x=Inf,y=Inf,hjust=0.8,vjust=-0.5)}  

####5. Arrange p2.cs and p2.sj using grid arrange (after doing #1-4 for each species)####
a <- ggplot_gtable(ggplot_build(p2.cs))
a$layout$clip[a$layout$name == "panel"] <- "off"
b <- ggplot_gtable(ggplot_build(p2.sj))
b$layout$clip[b$layout$name == "panel"] <- "off"
grid.arrange(a,b,ncol=2,widths=c(1,1.33),clip=F) 

#dev.copy2pdf(file="Figures/Fig4.pdf", width=8, height=4,fonts=c("Helvetica","sans")) #Optional: Write to file
#For the final version, I moved the legend inside one of the faceted plots, and changed the lines to be dashed lines, but I can't find the code I used to do that. It's possible it was done in Illustrator.
