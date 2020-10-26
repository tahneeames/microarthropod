#to add variables (diversity, richness, abundance)
library(vegan)
View(sitetaxa)
sitedata=read.table("~/Desktop/STANDDATA_tabdelim.txt",header=T,sep="\t")
View(sitedata)
sitedata=data.frame(sitedata,richness)
abundance=rowSums(sitetaxa)
sitedata=data.frame(sitedata,abundance)
diversity(sitetaxa,index="shannon")
shandiv=diversity(sitetaxa,index="shannon")
sitedata=data.frame(sitedata,shandiv)

#t tests of difference between treatment means 
#richnes
t.test(richness~standType,data=sitedata)

#output
	Welch Two Sample t-test
data:  richness by standType
t = 2.3172, df = 32.441, p-value = 0.02696
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.2869856 4.4402871
sample estimates:
mean in group Corridor Cut  mean in group Shelterwood 
                  8.000000                   5.636364 t

#abundance
.test(abundance~standType,data=sitedata)

#output

	Welch Two Sample t-test

data:  abundance by standType
t = 0.28205, df = 38.999, p-value = 0.7794
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -11.32410  14.99396
sample estimates:
mean in group Corridor Cut  mean in group Shelterwood 
                  28.78947                   26.95455 

#diversity
t.test(shandiv~standType,data=sitedata)

#output

	Welch Two Sample t-test

data:  shandiv by standType
t = 2.6199, df = 38.669, p-value = 0.0125
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.08923743 0.69444155
sample estimates:
mean in group Corridor Cut  mean in group Shelterwood 
                  1.648693                   1.256854 


#adonis plot code
adonis2(sitetaxa~epiphyteType*standType,data=sitedata)

#output
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = sitetaxa ~ epiphyteType * standType, data = sitedata)
                       Df SumOfSqs      R2      F Pr(>F)    
epiphyteType            1   0.9477 0.07767 3.8306  0.001 ***
standType               1   1.5780 0.12932 6.3779  0.001 ***
epiphyteType:standType  1   0.5223 0.04281 2.1112  0.032 *  
Residual               37   9.1544 0.75021                  
Total                  40  12.2025 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#to do the kruskal-wallis multiple comparison test
ibrary(pgirmess)
kruskalmc(sitedata$shandiv~sitedata$Factor)

#the output of the tests 
Multiple comparison test after Kruskal-Wallis 
p.value: 0.05 
Comparisons
                     obs.dif critical.dif difference
CorXBryo-CorXLich 11.4555556     14.52107      FALSE
CorXBryo-ShXBryo  17.1555556     14.52107       TRUE
CorXBryo-ShXLich  12.2222222     13.93609      FALSE
CorXLich-ShXBryo   5.7000000     14.13377      FALSE
CorXLich-ShXLich   0.7666667     13.53206      FALSE
ShXBryo-ShXLich    4.9333333     13.53206      FALSE

#code
kruskalmc(sitedata$richness~sitedata$Factor)

#output
Multiple comparison test after Kruskal-Wallis 
p.value: 0.05 
Comparisons
                    obs.dif critical.dif difference
CorXBryo-CorXLich 12.038889     14.52107      FALSE
CorXBryo-ShXBryo  14.538889     14.52107       TRUE
CorXBryo-ShXLich  13.347222     13.93609      FALSE
CorXLich-ShXBryo   2.500000     14.13377      FALSE
CorXLich-ShXLich   1.308333     13.53206      FALSE
ShXBryo-ShXLich    1.191667     13.53206      FALSE

#code
kruskalmc(sitedata$abundance~sitedata$Factor)

#output
Multiple comparison test after Kruskal-Wallis 
p.value: 0.05 
Comparisons
                    obs.dif critical.dif difference
CorXBryo-CorXLich 12.694444     14.52107      FALSE
CorXBryo-ShXBryo   6.544444     14.52107      FALSE
CorXBryo-ShXLich  12.819444     13.93609      FALSE
CorXLich-ShXBryo   6.150000     14.13377      FALSE
CorXLich-ShXLich   0.125000     13.53206      FALSE
ShXBryo-ShXLich    6.275000     13.53206      FALSE

#indicator species analysis by Factor
library(indicspecies)
inv = multipatt(sitetaxa, sitedata$Factor, func = "r.g", control = how(nperm=999))
summary(inv)

#output
 Multilevel pattern analysis
 ---------------------------

 Association function: r.g
 Significance level (alpha): 0.05

 Total number of species: 45
 Selected number of species: 10 
 Number of species associated to 1 group: 7 
 Number of species associated to 2 groups: 3 
 Number of species associated to 3 groups: 0 

 List of species associated to each combination: 

 Group CorXBryo  #sps.  5 
               stat p.value    
Isotomidae    0.693   0.001 ***
Ologamastidae 0.532   0.003 ** 
Thripidae     0.527   0.001 ***
Oppiidae      0.508   0.004 ** 
Bdelidae      0.434   0.034 *  

 Group ShXBryo  #sps.  1 
          stat p.value   
Acaridae 0.535   0.006 **

 Group ShXLich  #sps.  1 
            stat p.value  
Liacaridae 0.429   0.022 *

 Group CorXBryo+CorXLich  #sps.  1 
                stat p.value    
Hypogastruidae 0.636   0.001 ***

 Group CorXLich+ShXLich  #sps.  1 
            stat p.value  
Carbodidade 0.45   0.011 *

 Group ShXBryo+ShXLich  #sps.  1 
              stat p.value  
Oribatulidae 0.419    0.05 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

#indicator species analysis by epiphyte type 
 inv = multipatt(sitetaxa, sitedata$epiphyteType, func = "r.g", control = how(nperm=999))
summary(inv)

#output
 Multilevel pattern analysis
 ---------------------------

 Association function: r.g
 Significance level (alpha): 0.05

 Total number of species: 45
 Selected number of species: 8 
 Number of species associated to 1 group: 8 

 List of species associated to each combination: 

 Group Bryophyte  #sps.  6 
               stat p.value   
Isotomidae    0.388   0.004 **
Acaridae      0.384   0.008 **
Ctenacaridae  0.343   0.041 * 
Bdelidae      0.342   0.035 * 
Thripidae     0.326   0.002 **
Pyroglyphidae 0.320   0.030 * 

 Group Lichen  #sps.  2 
             stat p.value   
Carbodidade 0.449   0.002 **
Liacaridae  0.316   0.039 * 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

#indicator species analysis by treatment type
inv = multipatt(sitetaxa, sitedata$standType, func = "r.g", control = how(nperm=999))
summary(inv)

#output
 Multilevel pattern analysis
 ---------------------------

 Association function: r.g
 Significance level (alpha): 0.05

 Total number of species: 45
 Selected number of species: 9 
 Number of species associated to 1 group: 9 

 List of species associated to each combination: 

 Group Corridor Cut  #sps.  7 
                 stat p.value    
Hypogastruidae  0.633   0.001 ***
Isotomidae      0.388   0.003 ** 
Cymbaeremaeidae 0.343   0.016 *  
Ascidae         0.333   0.045 *  
Sminthuridae    0.332   0.037 *  
Thripidae       0.312   0.010 ** 
Oppiidae        0.304   0.043 *  

 Group Shelterwood  #sps.  2 
              stat p.value   
Oribatulidae 0.415   0.003 **
Acaridae     0.329   0.028 * 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

#to make the plots
library(ggplot2)
library(ggthemes)

#diversity boxplot
ggplot(sitedata, aes(x=epiphyteType, y=shandiv,fill=species))
+geom_boxplot(fill="grey97")+
+     geom_dotplot(binaxis='y', stackdir='center',dotsize = .65)
+facet_grid(~standType)
+labs(title="Arthropod Diversity Between Epiphytes by Stand Type",x="Epiphyte", y = "Shannon Diversity")
+ scale_fill_manual(values=c("turquoise4", "palegreen4", "darkolivegreen1", "slategray1")) 
+ theme_tufte()
+theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),text = element_text(size = 14, family = "Tahoma"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10)) 

#to save the plots (plus transparency)
ggsave("arth_div.png",plot=last_plot(),bg="transparent")

#richness boxplot
ggplot(sitedata, aes(x=epiphyteType, y=richness,fill=species))
+geom_boxplot(fill="grey97")+
+     geom_dotplot(binaxis='y', stackdir='center',dotsize = .65)
+facet_grid(~standType)
+labs(title="Arthropod Richness Between Epiphytes by Stand Type",x="Epiphyte", y = "Richness")
+ scale_fill_manual(values=c("turquoise4", "palegreen4", "darkolivegreen1", "slategray1")) 
+ theme_tufte()
+theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),text = element_text(size = 14, family = "Tahoma"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10)) 
ggsave("arth_rich.png",plot=last_plot(),bg="transparent")

#abundance boxplot
ggplot(sitedata, aes(x=epiphyteType, y=abundance,fill=species))
+geom_boxplot(fill="grey97")+
+     geom_dotplot(binaxis='y', stackdir='center',dotsize = .65)
+facet_grid(~standType)
+labs(title="Arthropod Abundance Between Epiphytes by Stand Type",x="Epiphyte", y = "Abundance")
+ scale_fill_manual(values=c("turquoise4", "palegreen4", "darkolivegreen1", "slategray1")) 
+ theme_tufte()
+theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),text = element_text(size = 14, family = "Tahoma"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10)) 
ggsave("arth_abund.png",plot=last_plot(),bg="transparent")


#NMDS Plots

#Factor NMDS
NMDS=metaMDS(sitetaxa,distance="binomial",trymax=999)
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Epiphyte = sitedata$Factor)

#NMDS Plot
ggplot(NMDS, aes(x=MDS1, y=MDS2, col=Epiphyte)) 
+geom_point() +
stat_ellipse()
+labs(title="NMDS of Bray Distance Between Epiphytes by Treatment", x="NMDS 1", y="NMDS 2")
+scale_color_discrete(name="Epiphyte",labels=c( "Corridor Cut Bryophyte", "Corridor Cut Lichen", "Shelterwood Bryophyte", "Shelterwood Lichen")) 
+theme_bw() 
+theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),text = element_text(size = 12, family = "Tahoma"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
+theme(legend.justification=c(0,0), legend.position=c(0,0))
+theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))

ggplot(NMDS, aes(x=MDS1, y=MDS2, col=Epiphyte)) 
+geom_point() 
+stat_ellipse()
+labs(title="NMDS of Bray Distance Between Epiphytes by Treatment", x="NMDS 1", y="NMDS 2")
+scale_color_discrete(name="Epiphyte",labels=c( "Corridor Cut Bryophyte", "Corridor Cut Lichen", "Shelterwood Bryophyte", "Shelterwood Lichen")) 
+theme_bw() 
+theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),text = element_text(size = 12, family = "Tahoma"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
+ theme(legend.justification=c(1,0), legend.position=c(1,0))
+theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))
+xlim(NA,10)

ggplot(NMDS, aes(x=MDS1, y=MDS2, col=Epiphyte)) 
+geom_point() 
+stat_ellipse()
+labs(title="NMDS of Bray Distance Between Epiphytes by Treatment", x="NMDS 1", y="NMDS 2")
+scale_color_discrete(name="Epiphyte",labels=c( "Corridor Cut Bryophyte", "Corridor Cut Lichen", "Shelterwood Bryophyte", "Shelterwood Lichen")) 
+theme_bw() 
+theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),text = element_text(size = 12, family = "Tahoma"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())
+ theme(legend.justification=c(1,1), legend.position=c(1,1))
+theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))+xlim(NA,10)+theme(panel.background = element_rect(fill = "transparent",colour = "transparent"))
+ theme(plot.background = element_rect(fill = "transparent"))
ggsave("NMDS_epiphyte_topleg.png",plot=last_plot(),bg="transparent")

#Treatment NMDS
NMDS=metaMDS(sitetaxa,distance="binomial",trymax=999)
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Epiphyte = sitedata$standType)

#NMDS Plot
ggplot(NMDS, aes(x=MDS1, y=MDS2, col=Treatment)) 
+geom_point() +stat_ellipse()
+labs(title="NMDS of Bray Distance Between Treatments", x="NMDS 1", y="NMDS 2")
+scale_color_manual(name="Treatment",labels=c( "Corridor Cut", "Shelterwood"),values=c("Corridor Cut"="cornflowerblue","Shelterwood"="darkslategray3")) 
+theme_bw() +theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),text = element_text(size = 12, family = "Tahoma"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())
+ theme(legend.justification=c(1,1), legend.position=c(1,1))
+theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))
+xlim(NA,10)+theme(panel.background = element_rect(fill = "transparent",colour = "transparent"))
+ theme(plot.background = element_rect(fill = "transparent"))
ggsave("NMDS_treatment.png",plot=last_plot(),bg="transparent")

#Epiphyte type NMDS
NMDS=metaMDS(sitetaxa,distance="binomial",trymax=999)
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Epiphyte = sitedata$epiphyteType)

#NMDS Plot
ggplot(NMDS, aes(x=MDS1, y=MDS2, col=Epiphyte))
+geom_point()+stat_ellipse()
+labs(title="NMDS of Bray Distance Between Epiphytes", x="NMDS 1", y="NMDS 2")
+scale_color_manual(name="Epiphyte Type",labels=c( "Bryophyte", "Lichen"),values=c("Bryophyte"="darkslateblue","Lichen"="lightsteelblue2"))
+theme_bw()
+theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),text = element_text(size = 12, family = "Tahoma"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())
+ theme(legend.justification=c(1,1), legend.position=c(1,1))
+theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))
+xlim(NA,10)+theme(panel.background = element_rect(fill = "transparent",colour = "transparent"))
+ theme(plot.background = element_rect(fill = "transparent"))
ggsave("NMDS_bryo_lich.png", plot=last_plot(),bg="transparent")


