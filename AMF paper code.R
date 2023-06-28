amf<-phyloseq(otu_table(abun.2, taxa_are_rows = TRUE), sample_data(sam), tax_table(tax.1))

nach<-subset_samples(amf, Site == 'Nachusa')

nach.1<-prune_taxa(taxa_sums(nach) > 2, nach)

#venn diagrams 
shared<-read.csv('C:/Users/JBell/Dropbox/NSF-ROL-MOR/PAPERS-AND-PROJECTS/Jennifer/AMF/shared.taxa.csv')
names(shared)[1]<-"Nachusa"

ggVennDiagram(shared, label_size = 8, set_size = 0) + scale_fill_distiller(palette = "RdBu") +
  theme(legend.position="none")

#community compostion
bray.all<-vegdist(t(nach.abun))
perm.all<-adonis2(bray.all~Bison+Burn+Age.cat+Type+Year, data=n.high.sam)
perm.all  

bray.high<-vegdist(high.nach.abun, method ='bray') 
perm.high<-adonis2(bray.high~Bison+Burn+Type+Year, data=n.high.sam)
perm.high  
  
bray.low<-vegdist(n.low.abun, method ='bray')
perm.low<-adonis2(bray.low~Bison+Burn+Age.cat+Type+Year, data=n.high.sam)
perm.low

no.ag.bray<-vegdist(t(no.ag.abun), method ='bray') 
perm.no.ag<-adonis2(no.ag.bray~Bison+Burn+Type+Year, data=no.ag.sam)
perm.no.ag  

#soil properties
no.ag.sam<-subset(n.high.sam, Type != 'ag')

aov.a<-aov(n.high.sam$C.N~n.high.sam$Bison) #repeat with other soil variables 
summary(aov.a)
boxplot(n.high.sam$percN~n.high.sam$Bison)

aov.b<-aov(n.high.sam$C.N~n.high.sam$Burn) #repeat with other soil variables
summary(aov.b)
boxplot(n.high.sam$percN~n.high.sam$Burn)

aov.c<-aov(n.high.sam$C.N~n.high.sam$Age.cat) #repeat with other soil variables
summary(aov.c)
boxplot(n.high.sam$percN~n.high.sam$Age.cat)

#diversity 
#all
all.div<-estimate_richness(nach.1, split=TRUE)

all.div<-read.csv('C:/Users/JBell/Dropbox/NSF-ROL-MOR/PAPERS-AND-PROJECTS/Jennifer/AMF/all.diversity.csv')
all.div$Year<-as.factor(all.div$Year)
all.div$Type<-as.factor(all.div$Type)
all.div$Bison<-as.factor(all.div$Bison)
all.div$Age.cat<-as.factor(all.div$Age.cat)
all.div$Burn<-as.factor(all.div$Burn)

aov7<-aov(all.div$InvSimpson~all.div$Age.cat)
summary(aov7)


ag.div<-subset(all.div, Type == 'ag')

aov8<-aov(ag.div$Observed~ag.div$Bison)
summary(aov8)

#no ag
no.ag.div<-subset(all.div, Type != 'ag')
aov9<-aov(no.ag.div$InvSimpson~no.ag.div$Age.cat) #repeat with other diversity measures 
summary(aov9)

#high
high.div<-estimate_richness(high.amf, split=TRUE)
high.div$Year<-as.factor(high.div$Year)
high.div$Type<-as.factor(high.div$Type)
high.div$Bison<-as.factor(high.div$Bison)
high.div$Age.cat<-as.factor(high.div$Age.cat)
high.div$Burn<-as.factor(high.div$Burn)

aov10<-aov(high.div$Simpson~high.div$Burn)
summary(aov10)

boxplot(high.div$Shannon~high.div$Burn)
TukeyHSD(aov10)

aov11<-aov(high.div$InvSimpson~high.div$Year)
summary(aov11)
boxplot(high.div$Observed~high.div$Year)

aov12<-aov(high.div$Fisher~high.div$Age.cat)
summary(aov12)

aov16<-aov(high.div$Fisher~high.div$Bison)
summary(aov16)

no.ag.div<-subset(high.div, Type !='ag')

aov17<-aov(no.ag.div$Simpson~no.ag.div$Burn)
summary(aov17)
TukeyHSD(aov17)

boxplot(no.ag.div$Shannon~no.ag.div$Burn)

#low
low.div<-estimate_richness(nach.low, split=TRUE)
low.abun<-data.frame(otu_table(nach.low))

low.div$Year<-as.factor(low.div$Year)
low.div$Type<-as.factor(low.div$Type)
low.div$Bison<-as.factor(low.div$Bison)
low.div$Age.cat<-as.factor(low.div$Age.cat)
low.div$Burn<-as.factor(low.div$Burn)


aov5<-aov(low.div$InvSimpson~low.div$Type) #repeat with other diversity measures 
summary(aov5)
boxplot(low.div$InvSimpson~low.div$Type)
TukeyHSD(aov5)

aov6<-aov(low.div$InvSimpson~low.div$Year) #repeat with other diversity measures 
summary(aov6)
boxplot(low.div$InvSimpson~low.div$Year)

aov7<-aov(low.div$InvSimpson~low.div$Age.cat) #repeat with other diversity measures 
summary(aov7)
boxplot(low.div$InvSimpson~low.div$Age.cat)
TukeyHSD(aov7)


#DeSeq2
NB_sig <- differential_abundance(physeq, grouping_column = "Burn",output_norm=NULL,pvalue.threshold=0.05,
                                 lfc.threshold=0)


diff.bison <- differential_abundance(physeq, grouping_column = "Bison",output_norm=NULL,pvalue.threshold=0.05,
                                     lfc.threshold=0)

diff.year <- differential_abundance(physeq, grouping_column = "Year",output_norm=NULL,pvalue.threshold=0.05,
                                    lfc.threshold=0,
                                    filename= 'C:/Users/JBell/Dropbox/NSF-ROL-MOR/PAPERS-AND-PROJECTS/Jennifer/AMF/year.deseq.csv')


diff.type <- differential_abundance(high.amf, grouping_column = "Type",output_norm=NULL,pvalue.threshold=0.05,
                                    lfc.threshold=0)

low.year<- differential_abundance(no.ag.low, grouping_column = "Year",output_norm=NULL,pvalue.threshold=0.05,
                                  lfc.threshold=0,
                                  filename= 'C:/Users/JBell/Dropbox/NSF-ROL-MOR/PAPERS-AND-PROJECTS/Jennifer/AMF/low.bison.deseq.csv')

diff.age<- differential_abundance(physeq, grouping_column = "Age.cat",output_norm=NULL,pvalue.threshold=0.05,
                                  lfc.threshold=0)


#remake the deseq graphs from above but prettier 

#bison
bison<-read.csv('C:/Users/JBell/Dropbox/NSF-ROL-MOR/PAPERS-AND-PROJECTS/Jennifer/AMF/High/abun.deseq.bison.csv')

ggplot(bison, aes(x= Bison, y = KJ960134, fill = Bison)) +
  geom_boxplot()+
  theme_classic() +
  scale_fill_manual(values = c('#D55E00', '#0072b2')) +
  theme(legend.position = "none", axis.text=element_text(size=20), axis.title = element_blank(), plot.title = element_text(size=20))+
  ggtitle("KJ960134")

aov12<-aov(bison$KJ960134~bison$Bison)
summary(aov12)
#p = 0.00255

ggplot(bison, aes(x= Bison, y = HG976190, fill = Bison)) +
  geom_boxplot()+
  theme_classic() +
  scale_fill_manual(values = c('#D55E00', '#0072b2')) +
  theme(legend.position = "none", axis.text=element_text(size=20), axis.title = element_blank(), plot.title = element_text(size=20))+
  ggtitle("HG976190")



#burn
burn<-read.csv('C:/Users/JBell/Dropbox/NSF-ROL-MOR/PAPERS-AND-PROJECTS/Jennifer/AMF/High/abun.deseq.burn.csv')

names(burn)

ggplot(burn, aes(x= Burn, y = LN620499, fill = Burn)) +
  geom_boxplot()+
  theme_classic() +
  scale_fill_manual(values = c('#009e73', '#56b4e9', '#999999')) +
  theme(legend.position = "none", axis.text=element_text(size=20), axis.title = element_blank(), plot.title = element_text(size=20))+
  ggtitle("LN620499")

aov14<-aov(burn$JX488680~burn$Burn)
summary(aov14)

ggplot(burn, aes(x= Burn, y = LT216954, fill = Burn)) +
  geom_boxplot()+
  theme_classic() +
  scale_fill_manual(values = c('#009e73', '#56b4e9', '#999999')) +
  theme(legend.position = "none", axis.text=element_text(size=20), axis.title = element_blank(), plot.title = element_text(size=20))+
  ggtitle("LT216954")


#site type
site<-read.csv('C:/Users/JBell/Dropbox/NSF-ROL-MOR/PAPERS-AND-PROJECTS/Jennifer/AMF/High/abun.deseq.type.csv')
names(site)

ggplot(site, aes(x= Type, y = KF386309, fill = Type)) +
  geom_boxplot()+
  theme_classic() +
  scale_fill_manual(values = c('#882255', '#999933', '#332288')) +
  theme(legend.position = "none", axis.text=element_text(size=20), axis.title = element_blank(), plot.title = element_text(size=20))+
  ggtitle("KF386309")

ggplot(site, aes(x= Type, y = KJ960140, fill = Type)) +
  geom_boxplot()+
  theme_classic() +
  scale_fill_manual(values = c('#882255', '#999933', '#332288')) +
  theme(legend.position = "none", axis.text=element_text(size=20), axis.title = element_blank(), plot.title = element_text(size=20))+
  ggtitle("KJ960140")


age<-read.csv('C:/Users/JBell/Dropbox/NSF-ROL-MOR/PAPERS-AND-PROJECTS/Jennifer/AMF/High/abun.deseq.age.csv')
names(age)
age$Age<-as.factor(age$Age)

ggplot(age, aes(x= Age, y = FR716551, fill = Age)) +
  geom_boxplot()+
  theme_classic() +
  scale_fill_manual(values = c('#44aa99', '#661100', '#6699cc', '#888888')) +
  theme(legend.position = "none", axis.text=element_text(size=20), axis.title = element_blank(), plot.title = element_text(size=20))+
  ggtitle("FR716551")

ggplot(age, aes(x= Age, y = LT158531, fill = Age)) +
  geom_boxplot()+
  theme_classic() +
  scale_fill_manual(values = c('#44aa99', '#661100', '#6699cc', '#888888')) +
  theme(legend.position = "none", axis.text=element_text(size=20), axis.title = element_blank(), plot.title = element_text(size=20))+
  ggtitle("LT158531")


#environmental correlation 
bison.cor <- taxa.env.correlation(physeq,"Bison")
plot_taxa_env(bison.cor)

age.cor <- taxa.env.correlation(physeq,"Age.cat")
plot_taxa_env(age.cor)


no.cor<-taxa.env.correlation(physeq,"Site")
plot_taxa_env(no.cor)

no.ag.sam<-data.frame(sample_data(no.ag))

no.fall<-subset_samples(no.ag, Burn != 'fall')

burn.cor<-taxa.env.correlation(no.fall,"Burn")
plot_taxa_env(burn.cor)
