#Code for Stevens and Latimer 2015
#Jens Stevens stevensjt@gmail.com 11/21/14
##This code analyzes soil temperature and soil moisture from Block 6 at Blodgett Forest.

library(reshape) #for melt.data.frame(). version 0.8.6
library(ggplot2) #for ggplot(). version 2.2.1
library(gridExtra) #for grid.arrange(). version 2.2.1

#setDefaults('as.Date.character', format = '%d/%M/%Y')
d<- read.csv("data/BL_B6_SoilProbes.csv",  stringsAsFactors=FALSE)
names(d)=c("DT","SM+","T+","SM-","T-","SMA","TA")

#Set to the data to one daytime and one nighttime reading
day=subset(d,subset=grepl("3:00 PM",d$DT))
night=subset(d,subset=grepl("3:00 AM",d$DT))
day$DT=as.Date(sapply(strsplit(day$DT," 3",fixed = TRUE),"[[",1),format="%m/%d/%Y") #Only take the first part [[1 of the string
night$DT=as.Date(sapply(strsplit(night$DT," 3",fixed = TRUE),"[[",1),format="%m/%d/%Y")

#d=day
d.SM.am=melt.data.frame(data=night,measure.vars=c("SM+","SM-","SMA"))
d.SM.pm=melt.data.frame(data=day,measure.vars=c("SM+","SM-","SMA"))
levels(d.SM.pm[,5])=c("Increased","Decreased","Ambient")

d.T.am=melt.data.frame(data=night,measure.vars=c("T+","T-","TA"))
d.T.pm=melt.data.frame(data=day,measure.vars=c("T+","T-","TA"))

####Block 6####
p.SM.am=ggplot(d.SM.am, aes(DT, value)) + geom_line(aes(linetype=variable)) +
  scale_x_date(labels = date_format("%b-%Y"),limits=c(as.Date("2011-12-05"),as.Date("2012-11-05"))) + 
  xlab("") + ylab("Nighttime Soil Moisture (%VWC)")+
  theme_bw()
p.SM.pm=ggplot(d.SM.pm, aes(DT, value)) + geom_line(aes(linetype=variable)) +
  scale_x_date(labels = date_format("%b-%Y"),limits=c(as.Date("2011-12-05"),as.Date("2012-11-05"))) + 
  xlab("") + ylab(expression(atop("Daytime Soil Moisture", "(%VWC)")))+
  annotate("text",x=as.Date("2011-12-05"),y=Inf,label="(c)",vjust=1)+
  guides(linetype=guide_legend(title="Snowpack Treatment"))+
  theme_bw()
p.T.am=ggplot(d.T.am, aes(DT, value)) + geom_line(aes(linetype=variable)) +
  scale_x_date(labels = date_format("%b-%Y"),limits=c(as.Date("2011-12-05"),as.Date("2012-08-05"))) + 
  xlab("") + ylab(expression(atop("Nighttime Temperature ", ( degree~C))))+ 
  geom_hline(yintercept=0)+guides(linetype=F)+
  annotate("text",x=as.Date("2011-12-05"),y=Inf,label="(a)",vjust=1)+
  theme_bw()
p.T.pm=ggplot(d.T.pm, aes(DT, value)) + geom_line(aes(linetype=variable)) +
  scale_x_date(labels = date_format("%b-%Y"),limits=c(as.Date("2011-12-05"),as.Date("2012-08-05"))) + 
  xlab("") +ylab(expression(atop("Daytime Temperature ", ( degree~C))))+ 
  geom_hline(yintercept=0)+guides(linetype=F)+
  annotate("text",x=as.Date("2011-12-05"),y=Inf,label="(b)",vjust=1)+
  theme_bw()

#p.SM.am;p.SM.pm;p.T.am;p.T.pm


####Plots####
a2 <- ggplot_gtable(ggplot_build(p.T.am))
b2 <- ggplot_gtable(ggplot_build(p.T.pm))
c2 <- ggplot_gtable(ggplot_build(p.SM.pm))
grid.arrange(arrangeGrob(a2,b2,nrow=1),c2,nrow=2,clip=F)

dev.copy2eps(file="AppendixSoilProbes_FigS1.eps", width=10, height=6,fonts=c("Helvetica","sans"))
