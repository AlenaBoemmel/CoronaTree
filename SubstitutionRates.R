install.packages("devtools")
library(devtools)
install_github( 'emvolz/treedater')
library(treedater)
library(ape)
library(magrittr)
library(lubridate)
library(skyspline)
library(ggplot2)
library(ggsci)


# read the phylogenetic tree
coronaTree_all <- read.tree("data/corona_germany_trimmed.fasta.raxml.bestTree.txt")

# root the tree by the oldest sequence from China
id1 <- grep("MN908947",coronaTree_all$tip.label)
coronaTree_all$tip.label[id1]
coronaTree_all <- root(coronaTree_all,outgroup = coronaTree_all$tip.label[id1],resolve.root = T)
coronaTree_all <- ladderize(coronaTree_all,right=F)

#dates <- sampleYearsFromLabels(coronaTree_all$tip.label,delimiter='|',index=3)
## read the manually imputed dates of sampling 
dates.man <- read.delim("data/headers_dates_corona_all.txt",header=FALSE,row.names = 1,as.is=TRUE)
# add characters to be able to use the function sampleYearsFromLabels
dates.man2 <- paste(rownames(dates.man),dates.man[,1],sep="\\") 
dates <- sampleYearsFromLabels(dates.man2,delimiter='\\',index=2)
# assign the correct names matching coronaTree_all$tip.labels
names(dates) <- rownames(dates.man)

hist(dates)

coronaTree_strict <- dater(coronaTree_all, dates, 29000 ,clock='strict')

plot(coronaTree_strict , no.mar=T, cex = .2 )

plot(coronaTree_strict,
     edge.color = ifelse(coronaTree_strict$edge.length>0.05,"red3","royalblue"),
     edge.width = 2,
     show.tip.label = F)
title("Random colours and branch thickness")


rootToTipRegressionPlot(coronaTree_strict)

outliers <- outlierTips(coronaTree_strict, alpha = 0.20) 

coronaTree_sel <- ape::drop.tip( coronaTree_all, rownames(outliers[outliers$q < 0.20,]) )

coronaTree_sel_strict <- dater(coronaTree_sel, dates, 29000,clock='strict')

rootToTipRegressionPlot(coronaTree_sel_strict)

pdf("plots/coronaTree_sel_strict_dater.pdf",width=7,height=10)
plot(coronaTree_sel_strict , no.mar=T, cex = .2 )
dev.off()

ggtree(coronaTree_all)


coronaTree_uncor <- dater(coronaTree_all, dates, 29000,clock='uncorrelated')

coronaTree_uncor

saveRDS(coronaTree_uncor, "data/coronaTree_uncor.rds")

rootToTipRegressionPlot(coronaTree_uncor)

outliers <- outlierTips(coronaTree_uncor, alpha = 0.20) 
dim(outliers)
hist(outliers$rates)


coronaTree_add <- dater(coronaTree_all, dates, 29000,clock='additive')

coronaTree_add

#### merge with metadata ###
metadat.df <- read.csv2("data/metadata.all.manual.csv",header=TRUE,row.names = 1, as.is=TRUE)
rownames(metadat.df) <- metadat.df$header
## define location as city or country
metadat.df$location <- metadat.df$division
metadat.df[is.na(metadat.df$location),"location"] <- metadat.df[is.na(metadat.df$location),"country"]
metadat.df[metadat.df$location=="unknown","location"] <- "Germany"
metadat.df[is.na(metadat.df$division),"division"] <- "unknown"
metadat.df$label <- metadat.df$header
metadat.df$division <- as.factor(metadat.df$division)
metadat.df$location <- as.factor(metadat.df$location)
metadat.df$country <- factor(metadat.df$country, levels=names(sort(table(metadat.df$country),
                                                                       decreasing=T)))
metadat.tree <- metadat.df#[,c("header","country","division","date")]

coronaTree_strict_ph <- as.phylo(coronaTree_strict$intree)
## take the edge length from the strict tree
coronaTree_strict_ph$edge.length <- coronaTree_strict$edge.length
coronaTree_strict_ph <- as.treedata(coronaTree_strict_ph)

coronaTree_strict_meta <- full_join(coronaTree_strict_ph, metadat.tree, by = 'label')

id1 <- grep("MN908947",coronaTree_strict_meta@phylo$tip.label)
coronaTree_strict_meta@phylo <- root(coronaTree_strict_meta@phylo,outgroup = coronaTree_strict_meta@phylo$tip.label[id1],resolve.root = T)

coronaTree_strict.tab <- as_tibble(coronaTree_strict_meta)
rootnode(coronaTree.tab)


ggtree(coronaTree_strict_meta) + 
  layout_dendrogram() +
  geom_tippoint(aes(colour=country),alpha=0.7, size=3)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  ggtitle("Phylogenetic tree by country") +
  scale_size_manual(values = 1) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(colour="grey30", size = 6))
ggsave("plots/tree_strict_bycountry.pdf",width=10,height=7)


ggtree(coronaTree_strict_meta, aes(color=branch.length)) +
  layout_dendrogram() +
  scale_color_continuous(low='darkgreen', high='red') +
  theme(legend.position="right")





#########################
library("seqinr")
corona_all <- read.fasta(file = "data/corona_germany_trimmed.fasta")

all.sixmer <- sapply(corona_all,count, wordsize=6)
var.all.sixmer <- apply(all.sixmer,1,var)
mat.sixmer <- all.sixmer[var.all.sixmer!=0,]
mat.sixmer <- as.matrix(t(mat.sixmer))
PCA.sixmer <- prcomp(mat.sixmer, scale=TRUE, center = FALSE)



library(factoextra)
pdf("plots/PCAsixmer.scree.pdf")
fviz_eig(PCA.sixmer)
dev.off()

# Eigenvalues
eig.val <- get_eigenvalue(PCA.sixmer)
eig.val[1:10,]


# Results for Variables
res.var <- get_pca_var(PCA.sixmer)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

## first 10 PCs for plotting##
PCA.df <- data.frame(PCA.sixmer$x[,1:10], metadat.df[rownames(PCA.sixmer$x),])

library(viridis)
ggplot(PCA.df, aes(x=PC1,y=PC2,color=date))+
  geom_point(size=3, alpha=0.6)+
  scale_color_viridis(discrete = TRUE, option = "A")+
  scale_fill_viridis(discrete = TRUE) +
  theme_classic()
ggsave("plots/PCAsixmer.12.date.pdf")  



ggplot(PCA.df, aes(x=PC1,y=PC2,color=country))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  theme_classic()
ggsave("plots/PCAsixmer.12.country.pdf")  

ggplot(PCA.df, aes(x=PC2,y=PC3,color=country))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  theme_classic()
ggsave("plots/PCAsixmer.23.country.pdf")  

###### 9-mer #####
all.ninemer <- sapply(corona_all,count, wordsize=9)
var.all.ninemer <- apply(all.ninemer,1,var)
mat.ninemer <- all.ninemer[var.all.ninemer!=0,]
mat.ninemer <- as.matrix(t(mat.ninemer))

PCA.ninemer <- prcomp(mat.ninemer, scale=TRUE, center = FALSE)


pdf("plots/PCAninemer.scree.pdf")
fviz_eig(PCA.ninemer)
dev.off()

# Eigenvalues
eig.val <- get_eigenvalue(PCA.ninemer)
eig.val[1:10,]


# Results for Variables
res.var <- get_pca_var(PCA.ninemer)
#res.var$coord          # Coordinates
#res.var$contrib        # Contributions to the PCs
#res.var$cos2           # Quality of representation 
round(sort(abs(res.var$contrib[,1]),decreasing=TRUE),4)[1:100]


## first 10 PCs for plotting##
PCA.df <- data.frame(PCA.ninemer$x[,1:10], metadat.df[rownames(PCA.ninemer$x),])

library(viridis)
ggplot(PCA.df, aes(x=PC1,y=PC2,color=date))+
  geom_point(size=3, alpha=0.6)+
  scale_color_viridis(discrete = TRUE, option = "A")+
  scale_fill_viridis(discrete = TRUE) +
  theme_classic()
ggsave("plots/PCAninemer.12.date.pdf")  



ggplot(PCA.df, aes(x=PC1,y=PC2,color=country))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  theme_classic()
ggsave("plots/PCAninemer.12.country.pdf")  

ggplot(PCA.df, aes(x=PC2,y=PC3,color=country))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  ylim(-15,15) + 
  theme_classic()
ggsave("plots/PCAninemer.23.zoom.country.pdf")  


#########################

library(Biostrings)
library(ape)

malign <- readDNAMultipleAlignment(filepath = "data/corona_germany_trimmed.fasta") 

mat.align <- consensusMatrix(malign, baseOnly=FALSE)

var.align <- apply(mat.align,2,var)

mat.align[,which.min(var.align)]

entropy <- function(col) {
  n = sum(col)
  H = -((col["A"]+1)/n * log2((col["A"]+1)/n) + (1+col["C"])/n * log2((1+col["C"])/n) + (1+col["G"])/n * log2((1+col["G"])/n) + (1+col["T"])/n * log2((1+col["T"])/n))
  return(H)
}

## calculate the entropy of the mutations ##
entropy.mat <- apply(mat.align,2,entropy)

### top 100 positions with the highest variability ###
mat.align[c("A","C","G","T","N","-"),order(entropy.mat,decreasing=TRUE)[1:100]]

snps <- order(entropy.mat,decreasing=TRUE)[1:15]

snps.mat <- as.matrix(malign)[,snps]

snps.df <- data.frame(snps.mat)
colnames(snps.df) <- paste("pos",snps,sep="_")

#####################
PCAsixmer.df <- read.delim("data/PCAsixmer.df.dat",header=TRUE,row.names = 1)
PCAsixmer.df$country <- factor(PCAsixmer.df$country, levels=names(sort(table(PCAsixmer.df$country),
                                                                   decreasing=T)))

#length(intersect(rownames(PCAsixmer.df),rownames(snps.df)))

PCAsixmer.df2 <- data.frame(PCAsixmer.df, snps.df[rownames(PCAsixmer.df),])
dim(PCAsixmer.df2)

ggplot(PCAsixmer.df2, aes(x=PC1,y=PC2,color=pos_22316))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  theme_classic()

ggplot(PCAsixmer.df2, aes(x=PC2,y=PC3,color=pos_22316))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  theme_classic()


ggplot(PCAsixmer.df2, aes(x=PC2,y=PC3,color=country, shape=pos_22316))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  theme_classic()
ggsave("plots/PCAsixmer.23.pos22316.pdf")  


PCAninemer.df <- read.delim("data/PCAninemer.df.dat",header=TRUE,row.names = 1)
PCAninemer.df$country <- factor(PCAninemer.df$country, levels=names(sort(table(PCAninemer.df$country),
                                                                       decreasing=T)))

#length(intersect(rownames(PCAsixmer.df),rownames(snps.df)))

PCAninemer.df2 <- data.frame(PCAninemer.df, snps.df[rownames(PCAninemer.df),])
dim(PCAninemer.df2)

ggplot(PCAninemer.df2, aes(x=PC1,y=PC2,color=pos_27794))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  theme_classic()
ggplot(PCAninemer.df2, aes(x=PC1,y=PC2,color=pos_27795))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  theme_classic()
ggplot(PCAninemer.df2, aes(x=PC1,y=PC2,color=pos_1386))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  theme_classic()
ggplot(PCAninemer.df2, aes(x=PC1,y=PC2,color=pos_2834))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  theme_classic()

ggplot(PCAninemer.df2, aes(x=PC1,y=PC2,color=pos_22316))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  theme_classic()

ggplot(PCAninemer.df2, aes(x=PC2,y=PC3,color=pos_22316))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  theme_classic()


ggplot(PCAninemer.df2, aes(x=PC2,y=PC3,color=country, shape=pos_22316))+
  geom_point(size=3, alpha=0.6)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
#  xlim(-10,10)+
#  ylim(-3,3)+
  theme_classic()
ggsave("plots/PCAninemer.23.country.pos22316.pdf")  

ggsave("plots/PCAninemer.23.country.zoom.pos22316.pdf")  

metadat2.df <- read.csv2("data/metadata.all.manual2.csv",header=TRUE,row.names = 1, as.is=TRUE)
rownames(metadat2.df) <- metadat2.df$header
## define location as city or country
metadat2.df$location <- metadat2.df$division
metadat2.df[is.na(metadat2.df$location),"location"] <- metadat2.df[is.na(metadat2.df$location),"country"]
metadat2.df[metadat2.df$location=="unknown","location"] <- "Germany"
metadat2.df[is.na(metadat2.df$division),"division"] <- "unknown"
metadat2.df[is.na(metadat2.df$city),"city"] <- metadat2.df[is.na(metadat2.df$city),"location"]
metadat2.df$city <- as.factor(metadat2.df$city)
metadat2.df$label <- metadat2.df$header
metadat2.df$division <- as.factor(metadat2.df$division)
metadat2.df$location <- as.factor(metadat2.df$location)
metadat2.df$country <- factor(metadat2.df$country, levels=names(sort(table(metadat2.df$country),
                                                                   decreasing=T)))


setdiff(rownames(metadat2.df),rownames(snps.df))
metadat2.df <- data.frame(metadat2.df,snps.df[rownames(metadat2.df),])



coronaTree_all <- read.tree("data/corona_germany_trimmed.fasta.raxml.bestTree.txt")
id1 <- grep("MN908947",coronaTree_all$tip.label)
coronaTree_all$tip.label[id1]
coronaTree_all <- root(coronaTree_all,outgroup = coronaTree_all$tip.label[id1],resolve.root = T)
coronaTree_all <- ladderize(coronaTree_all,right=F)
 
coronaTree_all <- as.treedata(coronaTree_all)

coronaTree <- full_join(coronaTree_all, metadat2.df, by = 'label')
#coronaTree <- root(coronaTree,outgroup = "MN908947_|Severe_acute_respiratory_syndrome_coronavirus_2_isolate_Wuhan-Hu-1|_complete_genome|China",resolve.root = T)

id1 <- grep("MN908947",coronaTree@phylo$tip.label)
coronaTree@phylo <- root(coronaTree@phylo,outgroup = coronaTree@phylo$tip.label[id1],resolve.root = T)

pos.snp <- grep("pos_",colnames(metadat2.df),value=TRUE)

for (pos.s in pos.snp[1:10]){
ggtree(coronaTree) + 
  layout_dendrogram() +
  geom_tippoint(aes_string(colour=pos.s),alpha=0.7, size=3)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
 # ggtitle("Phylogenetic tree by country") +
  scale_size_manual(values = 1) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(colour="grey30", size = 6))
ggsave(paste0("plots/tree_all_by",pos.s,".pdf"),width=10,height=7)
}


to_drop <- metadat2.df[!(metadat2.df$country%in%c("Germany","Wuhan")),"label"]
coronaTree_reduced <- drop.tip(coronaTree, to_drop)

for (pos.s in pos.snp){
  ggtree(coronaTree_reduced) + 
    layout_dendrogram() +
    geom_tippoint(aes_string(colour=pos.s),alpha=0.7, size=4)+
    scale_fill_d3("category20") +
    scale_color_d3("category20") +
    # ggtitle("Phylogenetic tree by country") +
    scale_size_manual(values = 1) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(colour="grey30", size = 8))
  ggsave(paste0("plots/tree_Germany_by",pos.s,".pdf"),width=10,height=7)
}
###################
### Plotting #####

ggplot(metadat2.df,aes(x=country,fill=country))+
  geom_bar(stat="count", color="grey50", alpha=0.5)+
  scale_fill_d3("category20") +
  scale_color_d3("category20") +
  #rotate_x_text(angle = -90, align = 0.5, valign = 0) +
  theme(legend.position="none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill="white")) 
#  theme_classic()
ggsave("plots/countries_barplot.pdf",width=8,height=5)


library(ggrepel)
entropy.df <- data.frame(pos=1:length(entropy.mat),entropy=entropy.mat,large=as.factor(ifelse(entropy.mat>0.5,1,0)))
ggplot(entropy.df,aes(x=pos,y=entropy,fill=large,color=large))+
  geom_bar(stat="identity")+
  geom_text_repel(data=entropy.df[entropy.df$large==1,],aes(label=pos),color="#8B1A1A") + 
  scale_fill_manual(values=c("#999999","#8B1A1A"))+
  scale_color_manual(values=c("#999999","#8B1A1A"))+
  scale_y_discrete(breaks=c(0,0.5,1))+
  theme_void()
ggsave("plots/entropy_barplot.pdf",width=10,height=5)


################################
ref <- readDNAStringSet("data/reference_ncbi.fasta")

eu <- readDNAStringSet("data/ncbi_europe_sequences.fasta")
eu <- eu[which(width(eu)>29000)]
eu <- replaceAmbiguities(eu, new="N")

submat1 <- nucleotideSubstitutionMatrix(match = 0, mismatch = -1, baseOnly = TRUE) 


pa1 <-
  pairwiseAlignment(eu[1], ref, substitutionMatrix = submat1,
                    gapOpening = 0, gapExtension = 1, scoreOnly = F)

pa2 <-
  pairwiseAlignment(eu[3], ref, substitutionMatrix = submat1,
                    gapOpening = 0, gapExtension = 1, scoreOnly = F)



MA <- readDNAMultipleAlignment("data/ncbi_europe_sequences.musclealignment.clw",format="clustal")
alphabetFrequency(MA)

nmatch(MA)