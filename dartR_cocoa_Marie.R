
library(pacman)
pacman::p_load(adegenet, car, diveRsity, genepop, ggplot2, inbreedR, lattice, lme4, lmerTest, 
               magrittr, nlme, pegas, PopGenReport, poppr, RColorBrewer, tidyr, vegan,reshape,ggpubr)
install.packages("devtools")
devtools::install_github("kkeenan02/diveRsity")
# loading data
dat <- read.table("C:/Users/kalousovam/OneDrive - CZU v Praze/Dokumenty/cacao/filtered_loci_metadata.txt", header=TRUE, "")

# create genind object
cocoa.genind <- df2genind(X=dat[,c(7:3246)], ploidy=2, type="codom", sep = "/", NA.char="NA",
                 pop = dat$pop, loc.names = NULL, ind.names = dat$id, ncode = 1)
                 strata = as.data.frame(dat$country)
cocoa.other <- dat[c(1:6)]

#adding coordinates and strata
cocoa.genind$other$xy <- cocoa.other[c(3,4)]

cocoa.genind$other$pop_country <- cocoa.other[c(2,6)]

strata(cocoa.genind)<-other(cocoa.genind)$pop_country

######################################################
############# BASIC STATS ############################
######################################################

# Basic stats with poppr 
# use ? poppr () to check the output details  
# Shannon, Simpson, Hex
poppr_stats<- poppr(cocoa.genind, sample = 1000, strata = ~pop/country)
poppr_stats
write_xlsx(poppr_stats,"poppr2.xlsx") #export table to excel
misssing.data<-info_table(cocoa.genind, plot = TRUE)
misssing.data

# Basic stats with hierfstat 
cocoa.hierfstat <- genind2hierfstat(cocoa.genind, pop = NULL) #from genind to hierfstat

basic.stats.hierfstat<- basic.stats(cocoa.hierfstat,diploid = T, digits = 2)

hierfstat::Hs(cocoa.hierfstat) #mean population gene diversities

# use ? basic.stats () to check the out details
basic.stats.hierfstat$overall #view the overall values
summary(basic.stats.hierfstat$Ho)
summary(basic.stats.hierfstat$Hs)
summary(basic.stats.hierfstat$Fis)
boxplot(basic.stats.hierfstat$Hs)

# global FST and Fis
WC<-wc(cocoa.genind)
WC

# Fst per population
betas<- betas(cocoa.hierfstat, nboot = 200, diploid = T, betaijT = F)
print(betas)

# Fst pairwise 
WC84<- genet.dist(cocoa.genind, method = "WC84") # you can try with different method
cocoa.fst <- pairwise.neifst(cocoa.hierfstat)
library(corrplot)
corrplot(cocoa.fst,is.corr = F, method = "square", type = "lower")
heatmap(cocoa.fst,Colv = "Rowv")

# Basic PCA 

div <- summary(cocoa.genind)
x <- indpca(cocoa.hierfstat) 
plot(x, cex = 0.7)

######################################################
###############   AMOVA   ############################
######################################################

# Loading data from genalex format 
cocoa.genalex <- read.genalex("C:/Users/kalousovam/OneDrive - CZU v Praze/Dokumenty/cacao/genalex.csv")
strata(cocoa.genalex)<-other(cocoa.genind)$pop_country

# in hier = always define from globa to specific e.g ~continent/country/pop etc...
# AMOVA using ade4 method
amova.result <- poppr.amova(cocoa.genalex, hier = ~ country/pop, nperm = 1000)
amova.result
amova.test <- randtest(amova.result) # Test for significance
plot(amova.test)

# AMOVA using pegas method
amova.pegas <- poppr.amova(cocoa.genalex, ~country/pop, method = "pegas")
amova.pegas
amova.pegas$varcomp/sum(amova.pegas$varcomp)

# proportion os shared alleles between pops 
prop.shared.alles<-PopGenReport::pairwise.propShared(cocoa.genind)

#######################################################
###############   MANTEL TEST   #######################
######################################################

# Genetic distance 
cocoa.gen.dist <- dist.gene(x=cocoa.genind@tab, method="pairwise")
# Geographical distance
cocoa.geo.dist <- dist(x=cocoa.genind$other$xy, method="euclidean", diag=TRUE,
                       upper=TRUE)
# Mantel test
cocoa.mantel <- mantel.randtest(m1=cocoa.gen.dist, m2=cocoa.geo.dist, nrepet=1000)
cocoa.mantel # See text output
plot(avocado.mantel, nclass=30)

# Different implementation of Mantel test testing distance classes
cocoa.mantel.cor <- mantel.correlog(D.eco=cocoa.gen.dist,D.geo=cocoa.geo.dist,
                                    XY=NULL, n.class=0, break.pts=NULL, cutoff=FALSE, r.type="pearson",
                                    nperm=1000, mult="holm", progressive=TRUE)
cocoa.mantel.cor # See results for respective classes
summary(cocoa.mantel.cor)

mantel.cocoa.vegan<- mantel(cocoa.geo.dist,cocoa.gen.dist, method = "spearman", permutations = 100, na.rm = T)
mantel.cocoa.vegan

# in conclusion.... p < 0.001, mantel statistic r = 0.2758 

############################################
######GENOTYPE ACCUMULATION CURVE ##########
############################################

gac_cocoa_dart <- poppr::genotype_curve(cocoa.genind, sample = 10, quiet = TRUE)
p<-last_plot()
p+geom_smooth()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))+theme(axis.text.y = element_text(size = 12))+theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))+xlab("Numer of loci sampled")+ylab("Number of multilocus genotypes")

############################################
######### POPULATION STRUCTURE  ############
############################################

# UPGMA /NJ
cocoa.genpop<- genind2genpop(cocoa.genind, process.other=TRUE)

upgma.cocoa.pop <- aboot(x=cocoa.genpop, sample = 100, tree = "upgma", root = T, distance = "nei.dist")
upgma.cocoa.ind <- aboot(x=cocoa.genpop, sample = 100, tree = "upgma", root = T, distance = "nei.dist")

# DAPC
# number of PCA to be retained to further analysis 
best_a_score <- adegenet::optim.a.score(dapc(cocoa.genind, n.pca = 50, n.da = nPop(cocoa.genind))) # shows that the best no. of PCs is 24

# DAPC
cocoa.dapc = dapc(cocoa.genind, var.contrib = TRUE, scale = FALSE, n.pca = 11, n.da = nPop(cocoa.genind)) # set the n.pca to 24, as predicted by the optim.a.score() function above
levels(pop(cocoa.genind))

jBrewColors = brewer.pal(n = 8, name = "Set1")  # create pallette color or use manual codes
scatter(cocoa.dapc, clabel = F, cstar = 1, cell = 1, 
        col = c("#0C5BB0FF","#EE0011FF", "#15983DFF", "#EC579AFF", "#FA6B09FF", "#149BEDFF", "#A1C720FF", "#FEC10BFF", "#16A08CFF","#9A703EFF","black"), 
        cex = 2, solid = 0.8, legend = T, scree.pca = T,posi.pca = "bottomleft") 

# Plot like STRUCTURE
# Data handle 
dapc.results <- as.data.frame(cocoa.dapc$posterior)
dapc.results$pop <- pop(cocoa.genind)
dapc.results$indNames <- rownames(dapc.results)
head(dapc.results, 5)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
head(dapc.results, n = 6)
colnames(dapc.results) <- c("Original_Pop","Tree","Population","Ancestry")
#dapc.results$Original_Pop <- factor(dapc.results$Original_Pop, levels = c("pop1","pop2","pop3")) use this line if you wanna set different order in population list
#just be careful to write the pop names exactly as in data.results object!! otherwise you will have an error 

# create the barplot 
p.1 <- ggplot(dapc.results, aes(x=Tree, y=Ancestry, fill=Population))
p.1 <- p.1 + geom_bar(stat='identity') 
p.1 <- p.1 + scale_fill_manual(values = c("#FEC10BFF","#A1C720FF","#15983DFF","#149BEDFF","#0C5BB0FF","#EE0011FF","#EC579AFF","#FA6B09FF"))
# "#0C5BB0FF", "#EC579AFF", "#15983DFF", "#FA6B09FF", "#EE0011FF", "#149BEDFF", "#A1C720FF", "#FEC10BFF", "#16A08CFF","#9A703EFF","black")) 
p.1 <- p.1 + facet_grid(~dapc.results$Original_Pop, scales = "free")
p.1 <- p.1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))
p.1 <- p.1+ theme(legend.position = "top")
p.1


# PCA
x = tab(cocoa.genind, NA.method = "mean")
percent = pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,12),
        names.arg = round(percent, 1))

pca1 = dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)
ind_coords = as.data.frame(pca1$li)
colnames(ind_coords) = c("Axis1","Axis2","Axis3")
ind_coords$Ind = indNames(cocoa.genind)
ind_coords$Site = cocoa.genind$pop
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))
cols = brewer.pal(nPop(cocoa.genind), "Set1")
cols = c("#0C5BB0FF", "#EC579AFF", "#15983DFF", "#FA6B09FF", "#EE0011FF", "#149BEDFF", "#A1C720FF", "#FEC10BFF", "#16A08CFF","#9A703EFF","black")
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15))
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = T)+
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  labs(x = xlab, y = ylab)+
  ggtitle("")+
  ggtheme

### KMEANS 
## Determine number of cluster
maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(cocoa.genind, n.pca = 50, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)

# define the cluster based on lowest BIC value or elbown method
p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1


my_k <- 3:4 #set the number of cluster. In this example I set for 3 and 4 cluster

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(cocoa.genind, n.pca = 11, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(cocoa.genind, pop = grp_l[[i]]$grp, n.pca = 11 , n.da = my_k[i])
}

my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp

my_pal <- RColorBrewer::brewer.pal(n=9, name = "Dark2")

# create scatter plot 
p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=clist$retro_7)#you must create clist object from pophelper script =)
p2 <- p2 + scale_fill_manual(values=c(paste(clist$retro_7, "90", sep = "")))
p2 <- p2 + theme(legend.position = "top")
p2

# data handle to create bar plot based on posterior probabilities
tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- cocoa.other$pop # yo can change to <- cocoa.other$country if you want to see the groups based on Country
my_df <- tmp

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- cocoa.other$pop # if you change to cocoa.other$country before, you must do it here also, if not you will have an error
  
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

# set populations order 
my_df$Region <- factor(my_df$Region, levels = c("PanzosAV","CobanAV","Criollo","ChibaylLanquin","Ecuador"))

# create barplot 
p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity", alpha=1)
p3 <- p3 + facet_grid(K ~ Region, scales = "free_x", space = "fixed", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")+xlab("Individuals")
p3 <- p3 + theme(legend.position='none')
#p3 <- p3 + scale_color_brewer(palette="Dark2")
p3 <- p3 + scale_fill_manual(values=c(clist$retro_7))#you must create clist object from pophelper script =)
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3 <- p3 + scale_x_discrete(guide = guide_axis(check.overlap = TRUE))
p3

#create multiplot con BIC values, scatterplot and barplot 
ggarrange(ggarrange(p1, p2, ncol = 2, labels = c("A", "B")),p3, nrow = 2,labels = c("", "C"), heights = c(1, 1))



