#move to microarthropod directory
setwd("~/Desktop/rDirectory/microarthropod")

##reading in the txt files as a data frame 
#site data 
metadata=read.table("~/Desktop/Microarthropod/statSheets/metadata_dwALL.txt",header=T,sep="\t")
View(metadata)
#I separated the taxa counts in each sample, called taxadata here 
taxadata=read.table("~/Desktop/Microarthropod/statSheets/taxaData_ALL.txt",header=T,sep="\t")
View(taxadata)
#temperature data - environmental data rather than sample-specific 
temps=read.table("~/Desktop/Microarthropod/statSheets/tempData.txt",header=T,sep="\t")
View(temps)
#dry weight data, using n=20 per substrate 
bryodw=read.table("~/Desktop/Microarthropod/statSheets/bryophyteDW20rep.txt",header=T,sep="\t")
View(bryodw)

#load libraries
library(vegan)
library(emmeans)
library(multcompView)
library(indicspecies)
library(betapart)
library(psych)
library(ggplot2)
library(ggthemes)
library(ggnewscale)
library(gridExtra)
library(knitr)
library(patchwork)

#adding a treatment x sample type factor to data frame
treatment=c(rep("corridor bryophyte",32),rep("corridor lichen",24),rep("shelterwood bryophyte",32,),rep("shelterwood lichen",24))
metadata=data.frame(metadata,treatment)

##Making accumulation curves - To check sampling effort 
#separating data frame for accumulation curves 
df1=data.frame(taxadata[1:8,])
df2=data.frame(taxadata[9:16,])
df3=data.frame(taxadata[17:24,])
df4=data.frame(taxadata[25:32,])
df5=data.frame(taxadata[33:40,])
df6=data.frame(taxadata[41:48,])
df7=data.frame(taxadata[49:56,])
df8=data.frame(taxadata[57:64,])
df9=data.frame(taxadata[65:72,])
df10=data.frame(taxadata[73:81,])
df11=data.frame(taxadata[82:88,])
df12=data.frame(taxadata[89:96,])
df13=data.frame(taxadata[97:104,])
df14=data.frame(taxadata[105:112,])

##plotting the curves 
#basic graphs
specaccum(df1)->AC
plot(AC, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
specaccum(df2)->BC
plot(BC, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
specaccum(df3)->NC
plot(NC, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
specaccum(df4)->PC
plot(PC, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
specaccum(df5)->LC
plot(LC, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
specaccum(df6)->SC
plot(SC, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
specaccum(df7)->GC
plot(GC, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
specaccum(df8)->AS
plot(AS, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
specaccum(df9)->BS
plot(BS, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
specaccum(df10)->NS
plot(NS, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
specaccum(df11)->PS
plot(PS, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
specaccum(df12)->LS
plot(LS, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
specaccum(df13)->SS
plot(SS, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
specaccum(df14)->GS
plot(GS, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

##combo graphs, grouped by treatment / epiphyte type 
#corridor cut bryophytes
plot(AC, ci.type="poly", col="skyblue3",lwd=2, ci.lty=0, ci.col="lightblue", main="Corridor Cut Bryophytes",xlab="Replicates",ylab="Taxa")
plot(BC, ci.type="poly", col="mediumpurple3", lwd=2, ci.lty=0, ci.col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.25), add=TRUE)
plot(NC, ci.type="poly", col="plum3", lwd=2, ci.lty=0, ci.col=rgb(red=0.3, green=0.5, blue=1.0, alpha=0.25), add=TRUE)
plot(PC, ci.type="poly", col="cyan4", lwd=2, ci.lty=0, ci.col=rgb(red=0.2, green=0.3, blue=0.6, alpha=0.25), add=TRUE)
legend("bottomright",legend=c("Anomodon","Brachythecium","Neckera","Porella"),col=c("skyblue3","mediumpurple3","plum3","cyan4"),cex=.6,pch=20)
#shelterwood bryophytes
plot(AS, ci.type="poly", col="skyblue3",lwd=2, ci.lty=0, ci.col="lightblue", ylim=c(0,40),main="Shelterwood Bryophytes",xlab="Replicates",ylab="Taxa")
plot(BS, ci.type="poly", col="mediumpurple3", lwd=2, ci.lty=0, ci.col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.25), add=TRUE)
plot(PS, ci.type="poly", col="cyan4", lwd=2, ci.lty=0, ci.col=rgb(red=0.2, green=0.3, blue=0.6, alpha=0.25), add=TRUE)
plot(NS, ci.type="poly", col="plum3", lwd=2, ci.lty=0, ci.col=rgb(red=0.3, green=0.5, blue=1.0, alpha=0.25), add=TRUE)
legend("bottomright",legend=c("Anomodon","Brachythecium","Neckera","Porella"),col=c("skyblue3","mediumpurple3","plum3","cyan4"),cex=.6,pch=20)
#corridor cut lichens
plot(LC, ci.type="poly", col="goldenrod", lwd=2, ci.lty=0, ci.col=rgb(red=245, green=191, blue=66, alpha=95,maxColorValue=255),ylim=c(0,40),main="Corridor Cut Lichens",xlab="Replicate",ylab="Taxa")
plot(SC, ci.type="poly", col="coral3", lwd=2, ci.lty=0, ci.col=rgb(red=199, green=91, blue=28, alpha=95,maxColorValue=255),add=TRUE)
plot(GC, ci.type="poly", col="khaki2", lwd=2, ci.lty=0, ci.col=rgb(red=227, green=227, blue=129, alpha=95,maxColorValue=255),add=TRUE)
legend("bottomright",legend=c("Lobaria pulmonaria", "Lobaria quercizans","Phaeophyscia"),col=c("goldenrod","coral3","khaki2"),cex=.6,pch=20)
#shelterwood lichens
plot(LS, ci.type="poly", col="goldenrod", lwd=2, ci.lty=0, ci.col=rgb(red=245, green=191, blue=66, alpha=95,maxColorValue=255),ylim=c(0,40),main="Shelterwood Lichens",xlab="Replicate",ylab="Taxa")
plot(SS, ci.type="poly", col="coral3", lwd=2, ci.lty=0, ci.col=rgb(red=199, green=91, blue=28, alpha=95,maxColorValue=255),add=TRUE)
plot(GS, ci.type="poly", col="khaki2", lwd=2, ci.lty=0, ci.col=rgb(red=227, green=227, blue=129, alpha=95,maxColorValue=255),add=TRUE)
legend("bottomright",legend=c("Lobaria pulmonaria", "Lobaria quercizans","Phaeophyscia"),col=c("goldenrod","coral3","khaki2"),cex=.6,pch=20)
#Comparing corridor and shelterwood treatments for all sample types
plot(AC, ci.type="poly", col="darkolivegreen3",lwd=2, ci.lty=0, ci.col=rgb(red=83, green=148, blue=67, alpha=95,maxColorValue=255),main="Anomodon rugelii", xlab="Replicate",ylab="Taxa")
plot(AS, ci.type="poly", col="deepskyblue3", lwd=2, ci.lty=0, ci.col=rgb(red=74, green=146, blue=161, alpha=95,maxColorValue = 255), add=TRUE)
legend("bottomright",legend=c("Corridor Cut","Shelterwood"),col=c("darkolivegreen3","deepskyblue3"),cex=.6,pch=20)
plot(BC, ci.type="poly", col="darkolivegreen3",lwd=2, ci.lty=0, ci.col=rgb(red=83, green=148, blue=67, alpha=95,maxColorValue=255),main="Brachythecium oxycladon",ylim=c(0,40), xlab="Replicate",ylab="Taxa")
plot(BS, ci.type="poly", col="deepskyblue3", lwd=2, ci.lty=0, ci.col=rgb(red=74, green=146, blue=161, alpha=95,maxColorValue = 255), add=TRUE)
legend("bottomright",legend=c("Corridor Cut","Shelterwood"),col=c("darkolivegreen3","deepskyblue3"),cex=.6,pch=20)
plot(NC, ci.type="poly", col="darkolivegreen3",lwd=2, ci.lty=0, ci.col=rgb(red=83, green=148, blue=67, alpha=95,maxColorValue=255),main="Neckera pennata",ylim=c(0,40), xlab="Replicate",ylab="Taxa")
plot(NS, ci.type="poly", col="deepskyblue3", lwd=2, ci.lty=0, ci.col=rgb(red=74, green=146, blue=161, alpha=95,maxColorValue = 255), add=TRUE)
legend("bottomright",legend=c("Corridor Cut","Shelterwood"),col=c("darkolivegreen3","deepskyblue3"),cex=.6,pch=20)
plot(PC, ci.type="poly", col="darkolivegreen3",lwd=2, ci.lty=0, ci.col=rgb(red=83, green=148, blue=67, alpha=95,maxColorValue=255),main="Porella platyphylla",ylim=c(0,40), xlab="Replicate",ylab="Taxa")
plot(PS, ci.type="poly", col="deepskyblue3", lwd=2, ci.lty=0, ci.col=rgb(red=74, green=146, blue=161, alpha=95,maxColorValue = 255), add=TRUE)
legend("bottomright",legend=c("Corridor Cut","Shelterwood"),col=c("darkolivegreen3","deepskyblue3"),cex=.6,pch=20)
plot(LC, ci.type="poly", col="darkolivegreen3",lwd=2, ci.lty=0, ci.col=rgb(red=83, green=148, blue=67, alpha=95,maxColorValue=255),main="Lobaria pulmonaria",ylim=c(0,40), xlab="Replicate",ylab="Taxa")
plot(LS, ci.type="poly", col="deepskyblue3", lwd=2, ci.lty=0, ci.col=rgb(red=74, green=146, blue=161, alpha=95,maxColorValue = 255), add=TRUE)
legend("bottomright",legend=c("Corridor Cut","Shelterwood"),col=c("darkolivegreen3","deepskyblue3"),cex=.6,pch=20)
plot(SC, ci.type="poly", col="darkolivegreen3",lwd=2, ci.lty=0, ci.col=rgb(red=83, green=148, blue=67, alpha=95,maxColorValue=255),main="Lobaria quercizans",ylim=c(0,40), xlab="Replicate",ylab="Taxa")
plot(SS, ci.type="poly", col="deepskyblue3", lwd=2, ci.lty=0, ci.col=rgb(red=74, green=146, blue=161, alpha=95,maxColorValue = 255), add=TRUE)
legend("bottomright",legend=c("Corridor Cut","Shelterwood"),col=c("darkolivegreen3","deepskyblue3"),cex=.6,pch=20)
plot(GC, ci.type="poly", col="darkolivegreen3",lwd=2, ci.lty=0, ci.col=rgb(red=83, green=148, blue=67, alpha=95,maxColorValue=255),main="Phaeophyscia sp.",ylim=c(0,40), xlab="Replicate",ylab="Taxa")
plot(GS, ci.type="poly", col="deepskyblue3", lwd=2, ci.lty=0, ci.col=rgb(red=74, green=146, blue=161, alpha=95,maxColorValue = 255), add=TRUE)
legend("bottomright",legend=c("Corridor Cut","Shelterwood"),col=c("darkolivegreen3","deepskyblue3"),cex=.6,pch=20)

####looking at the environmental data
##t-test to compare mean temperatures between treatments 
#one for average daily highs, one for average daily lows 
t.test(maxTemp~treatment,data=temps)
t.test(minTemp~treatment,data=temps)

##dry weights
#sep dry weights by species 
AC20=data.frame(bryodw[1:20,])
BC20=data.frame(bryodw[21:40,])
NC20=data.frame(bryodw[41:60,])
PC20=data.frame(bryodw[61:80,])
AS20=data.frame(bryodw[81:100,])
BS20=data.frame(bryodw[101:120,])
NS20=data.frame(bryodw[121:140,])
PS20=data.frame(bryodw[141:160,])

skirt20=rbind(A20,B20)
bole20=rbind(N20,P20)

A20=rbind(AC20,AS20)
B20=rbind(BC20,BS20)
N20=rbind(NC20,NS20)
P20=rbind(PC20,PS20)

skirtCC20=rbind(AC20,BC20)
skirtSW20=rbind(AS20,BS20)
boleCC20=rbind(NC20,PC20)
boleSW20=rbind(NS20,PS20)

#descriptive stats 
describe(skirtCC20$dryWeight)
describe(skirtSW20$dryWeight)
describe(boleCC20$dryWeight)
describe(boleSW20$dryWeight)
describe(AC20$dryWeight)
describe(AS20$dryWeight)
describe(BC20$dryWeight)
describe(BS20$dryWeight)
describe(NC20$dryWeight)
describe(NS20$dryWeight)
describe(PC20$dryWeight)
describe(PS20$dryWeight)

#kw rank sums for each species 
kruskal.test(dryWeight ~ standType, data = A20)
kruskal.test(dryWeight ~ standType, data = B20)
kruskal.test(dryWeight ~ standType, data = N20)
kruskal.test(dryWeight ~ standType, data = P20)
kruskal.test(dryWeight ~ standType, data = skirt20)
kruskal.test(dryWeight ~ standType, data = bole20)

#violin plots for weights
mytheme=theme_tufte()+theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),text = element_text(size = 14, family = "Tahoma"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10),)

ggplot(P20, aes(x=standType, y=dryWeight)) + 
geom_violin(fill="#edffec",trim="FALSE")+
labs(title="Porella platyphylla",x="",y="Dry Weight (g)")+
mytheme+
stat_summary(fun=median, geom="point", size=2, color="#7eca9c")

ggplot(N20, aes(x=standType, y=dryWeight))+ 
geom_violin(fill="#edffec",trim="FALSE")+
labs(title="Neckera pennata",x="",y="Dry Weight (g)")+
mytheme+
stat_summary(fun=median, geom="point", size=2, color="#7eca9c")

ggplot(B20, aes(x=standType, y=dryWeight)) + 
geom_violin(fill="#edffec",trim="FALSE")+
labs(title="Brachythecium oxycladon",x="",y="Dry Weight (g)")+
mytheme+ 
stat_summary(fun=median, geom="point", size=2, color="#7eca9c")

ggplot(A20, aes(x=standType, y=dryWeight))+ 
geom_violin(fill="#edffec",trim="FALSE")+
labs(title="Anomodon rugelii",x="",y="Dry Weight (g)")+
mytheme+ 
stat_summary(fun=median, geom="point", size=2, color="#7eca9c")

#boxplots 
ggplot(A20, aes(x=standType, y=dryWeight))+ 
geom_boxplot(fill="#caf7e3")+
labs(title="Porella platyphylla",x="",y="Dry Weight (g)")+
mytheme

ggplot(B20, aes(x=standType, y=dryWeight))+
geom_boxplot(fill="#caf7e3")+
labs(title="Porella platyphylla",x="",y="Dry Weight (g)")+
mytheme

ggplot(N20, aes(x=standType, y=dryWeight))+ 
geom_boxplot(fill="#caf7e3")+
labs(title="Porella platyphylla",x="",y="Dry Weight (g)")+
mytheme

ggplot(P20, aes(x=standType, y=dryWeight))+ 
geom_boxplot(fill="#caf7e3")+
labs(title="Porella platyphylla",x="",y="Dry Weight (g)")+
mytheme

#calculating richness, adding column to metadata
rowSums(taxadata>0)->"richness"
metadata=data.frame(metadata,richness)
View(metadata)
#abundance
rowSums(taxadata)->"abundance"
metadata=data.frame(metadata,abundance)
View(metadata)
#diversity
diversity(taxadata,index="shannon")->"diversity"
metadata=data.frame(metadata,diversity)
View(metadata)

#new df, separating the epiphyte categories 
epiphytesep=rep(c("Basal","Bole","Lichen","Bkirt","Bole","Lichen"),times=c(16,16,24,16,16,24))
newFactor=rep(c("corridorSkirt","corridorBole","corridorLichen","shelterwoodSkirt","shelterwoodBole","shelterwoodLichen"),times=c(16,16,24,16,16,24))
metadataMossType=cbind(metadata,epiphytesep,newFactor)
View(metadataMossType)

##univariate ANOVA
#abundance 
abundAnova=aov(abundance~as.factor(standType)*as.factor(epiphytesep),data=metadataMossType)
summary(abundAnova)
#richness
richAnova=aov(richness~as.factor(standType)*as.factor(epiphytesep),data=metadataMossType)
summary(richAnova)
#diversity
divAnova=aov(diversity~as.factor(standType)*as.factor(epiphytesep),data=metadataMossType)
summary(divAnova)

##interaction plots
interaction.plot(x.factor=metadataMossType$standType,trace.factor=metadataMossType$epiphytesep,response  = metadataMossType$abundance,fun=mean,type="b",fixed=TRUE,leg.bty="o")
interaction.plot(x.factor=metadataMossType$standType,trace.factor=metadataMossType$epiphytesep,response  = metadataMossType$richness,fun=mean,type="b",fixed=TRUE,leg.bty="o")
interaction.plot(x.factor=metadataMossType$standType,trace.factor=metadataMossType$epiphytesep,response  = metadataMossType$diversity,fun=mean,type="b",fixed=TRUE,leg.bty="o")

##interaction was present in abundance and richness 
#Abundance
tukeyAbund <- TukeyHSD(abundAnova)
abundCLD <- multcompLetters4(abundAnova, tukeyAbund)
#view groups 
print(abundCLD)
#view significance values 
View(tukeyAbund)
#Richness
tukeyRich <- TukeyHSD(richAnova)
richCLD <- multcompLetters4(richAnova, tukeyRich)
print(richCLD)
View(tukeyRich)

##interaction was not present in diversity
#stand type 
standDiv=emmeans(divAnova,~standType)
pairs(standDiv,adjust="tukey")
#epiphyte type
epDiv = emmeans(divAnova,~ epiphytesep)
pairs(epDiv,adjust="tukey")
#tukey letters 
tukeyDiv <- TukeyHSD(divAnova)
divCLD <- multcompLetters4(divAnova, tukeyDiv)
print(divCLD)

#pairwise comparisons to compare variables between each specific treatment block 
pairwise.t.test(metadataMossType$richness,metadataMossType$partition)
pairwise.t.test(metadataMossType$abundance,metadataMossType$partition)
pairwise.t.test(metadataMossType$diversity,metadataMossType$partition)

##descriptive stats, for the magnitude of the significant differences in variables of interest for substrate
#subset the df into the six factors
corSkirt=data.frame(metadataMossType[1:16,])
corBole=data.frame(metadataMossType[17:32,])
corLichn=data.frame(metadataMossType[33:56,])
shelShirt=data.frame(metadataMossType[57:72,])
shelBole=data.frame(metadataMossType[73:88,])
shelLichn=data.frame(metadataMossType[89:112,])
#another subset for treatment differences 
corridorStats=data.frame(metadataMossType[1:56,])
shelterwoodStats=data.frame(metadataMossType[57:112,])
#another for epiphyte type
skirtStats=data.frame(rbind(corSkirt,shelShirt))
boleStats=data.frame(rbind(corBole,shelBole))
lichenStats=data.frame(rbind(corLichn,shelLichn))

#richness
mean(corSkirt$richness)
mean(corBole$richness)
mean(corLichn$richness)
mean(shelShirt$richness)
mean(shelBole$richness)
mean(shelLichn$richness)
#abundance
mean(corSkirt$abundance)
mean(corBole$abundance)
mean(corLichn$abundance)
mean(shelShirt$abundance)
mean(shelBole$abundance)
mean(shelLichn$abundance)
#diversity
mean(corSkirt$diversity)
mean(corBole$diversity)
mean(corLichn$diversity)
mean(shelShirt$diversity)
mean(shelBole$diversity)
mean(shelLichn$diversity)

#richness 
mean(corridorStats$richness)
mean(shelterwoodStats$richness)
#abundance
mean(corridorStats$abundance)
mean(shelterwoodStats$abundance)
#diversity
mean(corridorStats$diversity)
mean(shelterwoodStats$diversity)

#richness
mean(skirtStats$richness)
mean(boleStats$richness)
mean(lichenStats$richness)
#abundance
mean(skirtStats$abundance)
mean(boleStats$abundance)
mean(lichenStats$abundance)
#diversity
mean(skirtStats$diversity)
mean(boleStats$diversity)
mean(lichenStats$diversity)

#plot
abund_wrap<-ggplot(metadataMossType, aes(x=reorder(epiphytesep,-abundance), y=log10(abundance),fill=newFactor,color=newFactor))+
introdataviz::geom_split_violin(alpha=.6,trim = FALSE)+
 stat_summary(fun.data = "mean_se", geom = "pointrange",size=0.25, show.legend = F, 
               position = position_dodge(.4)) +
labs(x = element_blank(),y=(bquote("(B) "~log[10]~Abundance/25~cm^2~{})))+
scale_color_manual(values=c("#5D8852","#6D92B1","#477567","#6EA161","#87B4DA","#70A08F"),guide="none")+
scale_fill_manual(values=c("#7FB069","#90C2E7","#477567","#92CE7B","#9FDDF7","#70A08F"),guide="none")+
theme_tufte()+
theme(plot.title = element_text(size = 13, family = "Times", face = "bold"),
  text = element_text(size = 12, family = "Times"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect( size = 1, linetype = "solid"))+
theme(legend.position ="none")

rich_wrap<-ggplot(metadataMossType, aes(x=reorder(epiphytesep,-abundance), y=richness,fill=newFactor,color=newFactor))+
introdataviz::geom_split_violin(alpha=.6,trim = FALSE)+
 stat_summary(fun.data = "mean_se", geom = "pointrange",size=0.25, show.legend = F, 
               position = position_dodge(.4)) +
labs(x = element_blank(),y=(bquote("(A) "~Richness/25~cm^2~{})))+
scale_color_manual(values=c("#5D8852","#6D92B1","#477567","#6EA161","#87B4DA","#70A08F"),guide="none")+
scale_fill_manual(values=c("#7FB069","#90C2E7","#477567","#92CE7B","#9FDDF7","#70A08F"),guide="none")+
theme_tufte()+
theme(plot.title = element_text(size = 13, family = "Times", face = "bold"),
  text = element_text(size = 12, family = "Times"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect( size = 1, linetype = "solid"))+
theme(legend.position ="none")

div_wrap<-ggplot(metadataMossType, aes(x=reorder(epiphytesep,-abundance), y=diversity,fill=newFactor,color=newFactor))+
introdataviz::geom_split_violin(alpha=.6,trim = FALSE)+
 stat_summary(fun.data = "mean_se", geom = "pointrange", size=0.25,show.legend = F, 
               position = position_dodge(.4)) +
labs(x = element_blank(),y=(bquote("(C) "~Diversity/25~cm^2~{})))+
scale_color_manual(values=c("#5D8852","#6D92B1","#477567","#6EA161","#87B4DA","#70A08F"),guide="none")+
scale_fill_manual(values=c("#7FB069","#90C2E7","#477567","#92CE7B","#9FDDF7","#70A08F"),guide="none")+
theme_tufte()+
theme(plot.title = element_text(size = 13, family = "Times", face = "bold"),
  text = element_text(size = 12, family = "Times"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect( size = 1, linetype = "solid"))+
theme(legend.position ="none")

wrapped <- rich_wrap + abund_wrap +div_wrap 

##multivariate analyses of community composition
#PerMANOVA 
adonis2(taxadata~epiphytesep*standType,data=metadataMossType,distance="bray")
#PERMDISP 
dist=vegdist(taxadata)
#by epiphyte type 
permutest(betadisper(dist,metadataMossType$epiphytesep),permutations=99,pairwise=TRUE)
#by stand type
permutest(betadisper(dist,metadataMossType$standType),permutations=99,pairwise=TRUE)
#combined trmt x substrate
permutest(betadisper(dist,metadataMossType$newFactor),permutations=99,pairwise=TRUE)

#multilevel pattern analysis
indicators=multipatt(taxadata,metadataMossType$newFactor,func = "r.g", control = how(nperm=999))
summary(indicators)

#NMDS
NMDS=metaMDS(taxadata,distance="binomial",trymax=999)
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDSord= data.frame(MDS1 = MDS1, MDS2 = MDS2)

#plot NMDS 
ordPlot<-ggplot(NMDSord, aes(x=MDS1, y=MDS2))+
geom_point(aes(color=metadataMossType$newFactor, fill=metadataMossType$newFactor,shape=factor(metadataMossType$newFactor)),size=2.5)+
scale_fill_manual(values=alpha(c("#477567","#5D8852","#6D92B1","#70A08F","#6EA161","#87B4DA"),.6),
  name ="Sample",
  breaks=c("corridorSkirt","corridorBole","corridorLichen","shelterwoodSkirt","shelterwoodBole","shelterwoodLichen"),
  labels=c("Corridor Basal","Corridor Bole","Corridor Lichen","Shelterwood Basal","Shelterwood Bole","Shelterwood Lichen"))+
scale_color_manual(values=c("#477567","#5D8852","#6D92B1","#70A08F","#6EA161","#87B4DA"),
  name ="Sample",
  breaks=c("corridorSkirt","corridorBole","corridorLichen","shelterwoodSkirt","shelterwoodBole","shelterwoodLichen"),
  labels=c("Corridor Basal","Corridor Bole","Corridor Lichen","Shelterwood Basal","Shelterwood Bole","Shelterwood Lichen"))+
new_scale("color")+
scale_shape_manual(values=c(21,21,21,24,24,24),
  name="Sample",breaks=c("corridorSkirt","corridorBole","corridorLichen","shelterwoodSkirt","shelterwoodBole","shelterwoodLichen"),
  labels=c("Corridor Basal","Corridor Bole","Corridor Lichen","Shelterwood Basal","Shelterwood Bole","Shelterwood Lichen"))+
stat_ellipse(size=1,aes(x = MDS1,y=MDS2,color=metadataMossType$epiphytesep,linetype=factor(metadataMossType$standType)))+
scale_linetype_discrete(name="Treatment", labels=c("Corridor","Shelterwood"))+
scale_color_manual(values=c("#477567","#7FB069","#90C2E7"),name="Epiphyte", breaks=c("Basal","Bole","Lichen"),labels=c("Basal","Bole","Lichen"))+
labs(x="NMDS 1", y="NMDS 2")+
theme_bw()+
theme(plot.title = element_text(size = 15, family = "Times", face = "bold"),text = element_text(size = 12, family = "Times"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+ 
theme(legend.justification=c(1,1), legend.position=c(1,1))+
theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))+
xlim(NA,11)+theme(panel.background = element_rect(fill = "transparent",colour = "transparent"))+ 
theme(plot.background = element_rect(fill = "transparent"))

ordPlot
ggsave("textOrd2green.png",plot=last_plot(),bg="transparent")


##nestedness / turnover analysis
#new df - sum taxadata per factor
sumCC_Skirt = NA
for (i in 1:70){sumCC_Skirt = c(sumCC_Skirt, sum(taxadata[1:16, i]))}
sumCC_Skirt = sumCC_Skirt[2:71]

sumCC_Bole = NA
for (i in 1:70){sumCC_Bole = c(sumCC_Bole, sum(taxadata[17:32, i]))}
sumCC_Bole = sumCC_Bole[2:71]

sumCC_Lic = NA
for (i in 1:70){sumCC_Lic = c(sumCC_Lic, sum(taxadata[33:56, i]))}
sumCC_Lic = sumCC_Lic[2:71]

sumSW_Skirt = NA
for (i in 1:70){sumSW_Skirt = c(sumSW_Skirt, sum(taxadata[57:72, i]))}
sumSW_Skirt = sumSW_Skirt[2:71]

sumSW_Bole = NA
for (i in 1:70){sumSW_Bole = c(sumSW_Bole, sum(taxadata[73:88, i]))}
sumSW_Bole = sumSW_Bole[2:71]

sumSW_Lic = NA
for (i in 1:70){sumSW_Lic = c(sumSW_Lic, sum(taxadata[89:112, i]))}
sumSW_Lic = sumSW_Lic[2:71]
#making the data frame, and converting it to binary rather than abundances 
nestData = data.frame(rbind(sumCC_Skirt, sumCC_Bole, sumCC_Lich, sumSW_Skirt, sumSW_Bole, sumSW_Lich)); names(nestData) = names(taxadata)
nestData[nestData>0]<-1 
#analysis
beta.multi(nestData,index.family = "sorensen")

#quick plot 
plot(nestedtemp(nestData), kind="incid", weighted=TRUE, names=TRUE, col = c("#edf6f9", "#83c5be"))

#seeing the outcomes from taxa that have >1 occurences 
taxadata.1=read.table("~/Desktop/Microarthropod Project/Stats Sheets/taxadata>1.txt",header=T, sep ="\t")

CC_Skirt = NA
for (i in 1:54){CC_Skirt = c(CC_Skirt, sum(taxadata.1[1:16, i]))}
CC_Skirt = CC_Skirt[2:55]

CC_Bole = NA
for (i in 1:54){CC_Bole = c(CC_Bole, sum(taxadata.1[17:32, i]))}
CC_Bole = CC_Bole[2:55]

CC_Lic = NA
for (i in 1:54){CC_Lic = c(CC_Lic, sum(taxadata.1[33:56, i]))}
CC_Lic = CC_Lic[2:55]

SW_Skirt = NA
for (i in 1:54){SW_Skirt = c(SW_Skirt, sum(taxadata.1[57:72, i]))}
SW_Skirt = SW_Skirt[2:55]

SW_Bole = NA
for (i in 1:54){SW_Bole = c(SW_Bole, sum(taxadata.1[73:88, i]))}
SW_Bole = SW_Bole[2:55]

SW_Lic = NA
for (i in 1:54){SW_Lic = c(SW_Lic, sum(taxadata.1[89:112, i]))}
SW_Lic = SW_Lic[2:55]

nestData.1=data.frame(rbind(CC_Skirt,CC_Bole,CC_Lic,SW_Skirt,SW_Bole,SW_Lic));names(nestData.1)=names(taxadata.1)
plot(nestedtemp(nestData.1), kind="incid",  names=TRUE, col = c("#B4D6D3", "#456063"))

