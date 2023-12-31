### Working on further analyses after filtering

##First,load/call the filtered data

library(dartR)
gl <- gl.load("afrocarpus1_1.Rdata")  #specific sites used as 'pop' while preparing the individual meta file, 70 %call rate

gl_progeny <- gl[gl$other$ind.metrics$Stage!= "parent", ] # extracts and creates a new genlight object containing individuals of the 'progeny' growth stage only
gl_parent <- gl[gl$other$ind.metrics$Stage!= "progeny", ] # extracts and creates a new genlight object containing individuals of the 'parent' growth stage only


## Genetic distance
#a distance matrix can be or  constructed between individuals based on the genetic profiles or between populations based on their allele frequency profiles. 
#Distance matrices can be generated by a number of R packages, the most popular of which are dist() from package stats and vegdist from package vegan. The function gl.dist is a wrapper for those two functions, applying them to allele frequencies calculated for each locus at each population defined in the genlight object.

library (stats)
library (vegan)

gl.dist.ind(gl[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci  # for the whole data
gl.dist.ind(gl_parent[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci  # for the parents data
gl.dist.ind(gl_progeny[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci  # for the progeny data

dist_pop <- gl.dist.pop(gl) # euclidean distance between populations for the whole data
dist_pop2 <- gl.dist.pop(gl_parent) # euclidean distance between populations for the parents data
dist_pop3 <- gl.dist.pop(gl_progeny) # euclidean distance between populations for the progeny data

## F Statistics
#For Fst and Neis Gst, the functions of the StAMPP package allow for the use of parallel computing, which is much faster.
#It also works on genlight objects directly and therefore no conversion is necessary
#stamppFst calculates pairwise Fst values between populations
#
library(StAMPP)
pwfst <-stamppFst(gl[1:20,], nboots=1, percent=95, nclusters=1)
?stamppFst
#For Neis Gst use
pwGst <-stamppNeisD(gl[1:20,]) 
# threw error message
# 
## Principal Coordinates Analysis (PCoA) 

 pc <- gl.pcoa(gl, nfactors=5) # for the whole data
 names(pc)
 pc2 <- gl.pcoa(gl_parent, nfactors=5) # for the parent data
 pc3 <- gl.pcoa(gl_progeny, nfactors=5) # for the progeny data
 
  barplot(pc$eig/sum(pc$eig)*100, ) # percentage of variation that is represented by the axes for the whole data 
  barplot(pc2$eig/sum(pc2$eig)*100, ) # for the parent data
  barplot(pc3$eig/sum(pc3$eig)*100, ) # for the progeny data
  
  
d <- gl.dist.ind(gl, method="euclidean", scale=TRUE) #distance for the whole data
d2 <- gl.dist.ind(gl_parent, method="euclidean", scale=TRUE) #distance for the parent data
d3 <- gl.dist.ind(gl_progeny, method="euclidean", scale=TRUE) #distance for the progeny data
  pca <-  gl.pcoa(d) # pca based on distance matrix for whole data
  pca2 <- gl.pcoa(d2) # pca based on distance matrix for parent data
  pca3 <- gl.pcoa(d3) # pca based on distance matrix for progeny data
  pca
 
 #Plotting the results of PCoA
 #
 gl.pcoa.plot(pc, gl, pop.labels="pop", xaxis=1, yaxis=2) # whole data
 #axis 1 explained only 2.5% of the variation and axis explained only 1.7% of the variation
 
 gl.pcoa.plot(pc2, gl_parent, pop.labels="pop", xaxis=1, yaxis=2) # parent data
 #axis 1 explained only 3.4% of the variation and axis explained only 3.2% of the variation
 
 gl.pcoa.plot(pc3, gl_progeny, pop.labels="pop", xaxis=1, yaxis=2) # whole data
 #axis 1 explained only 3.3% of the variation and axis explained only 1.7% of the variation
 #
 gl.pcoa.plot(pca, gl, pop.labels="pop", xaxis=1, yaxis=2) # whole data
 #axis 1 explained only 2.8% of the variation and axis explained only 2.1% of the variation
 #
 gl.pcoa.plot(pca2, gl_parent, pop.labels="pop", xaxis=1, yaxis=2) #parent
 #axis 1 explained only 4.1% of the variation and axis explained only 3.9% of the variation
 gl.pcoa.plot(pca3, gl_progeny, pop.labels="pop", xaxis=1, yaxis=2) #progeny
 #axis 1 explained only 3.8% of the variation and axis explained only 2.9% of the variation
 #
 
## Neighbour-joining trees
## 
 gl.tree.nj(gl,type = "fan") #whole data
 gl.tree.nj(gl_parent,type = "fan") #parent
 gl.tree.nj(gl_progeny,type = "fan") #progeny
 
 gl.tree.nj(gl,type = "phylogram") #whole data
 gl.tree.nj(gl_parent,type = "phylogram") #parent
 gl.tree.nj(gl_progeny,type = "phylogram") #progeny

 
 ##First,load/call the filterd data
 
 library(dartR)
 gl <- gl.load("afrocarpus5.Rdata")  #using 80% call rate 
 
 
 ## Genetic distance
 #a distance matrix can be or  constructed between individuals based on the genetic profiles or between populations based on their allele frequency profiles. 
 #Distance matrices can be generated by a number of R packages, the most popular of which are dist() from package stats and vegdist from package vegan. The function gl.dist is a wrapper for those two functions, applying them to allele frequencies calculated for each locus at each population defined in the genlight object.
 library (stats)
 library (vegan)
 
 gl.dist.ind(gl[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci
 
 dist_pop <- gl.dist.pop(gl) # euclidean distance between populations
 
 ## F Statistics
 #For Fst and Neis Gst, the functions of the StAMPP package allow for the use of parallel computing, which is much faster.
 #It also works on genlight objects directly and therefore no conversion is necessary
 #stamppFst calculates pairwise Fst values between populations
 #
 library(StAMPP)
 pwfst <-stamppFst(gl[1:20,], nboots=1, percent=95, nclusters=1)
 ?stamppFst
 #For Neis Gst use
 pwGst <-stamppNeisD(gl[1:20,]) 
 # threw error message
 # 
 ## Principal Coordinates Analysis (PCoA) 
 
 pc <- gl.pcoa(gl, nfactors=5)
 names(pc)
 
 barplot(pc$eig/sum(pc$eig)*100, ) # percentage of variation that is represented by the axes 
 
 d <- gl.dist.ind(gl, method="euclidean", scale=TRUE) 
 pca <-  gl.pcoa(d) # pca based on distance matrix
 pca
 
 #Plotting the results of PCoA
 #
 gl.pcoa.plot(pc, gl, pop.labels="pop", xaxis=1, yaxis=2)
 #axis 1 explained only 3.1% of the variation and axis explained only 3.1% of the variation
 
 gl.pcoa.plot(pca, gl, pop.labels="pop", xaxis=1, yaxis=2) 
 #axis 1 explained only 3.4% of the variation and axis explained only 3.2% of the variation
 #
 ## Neighbour-joining trees
 ## 
 gl.tree.nj(gl,type = "fan")
 gl.tree.nj(gl,type = "phylogram")
 
  
 
 ##First,load/call the filterd data
 
 library(dartR)
 gl <- gl.load("afrocarpus6.Rdata")  #using 90% call rate 
 
 
 ## Genetic distance
 #a distance matrix can be or  constructed between individuals based on the genetic profiles or between populations based on their allele frequency profiles. 
 #Distance matrices can be generated by a number of R packages, the most popular of which are dist() from package stats and vegdist from package vegan. The function gl.dist is a wrapper for those two functions, applying them to allele frequencies calculated for each locus at each population defined in the genlight object.
 library (stats)
 library (vegan)
 
 gl.dist.ind(gl[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci
 
 dist_pop <- gl.dist.pop(gl) # euclidean distance between populations
 
 ## F Statistics
 #For Fst and Neis Gst, the functions of the StAMPP package allow for the use of parallel computing, which is much faster.
 #It also works on genlight objects directly and therefore no conversion is necessary
 #stamppFst calculates pairwise Fst values between populations
 #
 library(StAMPP)
 pwfst <-stamppFst(gl[1:20,], nboots=1, percent=95, nclusters=1)
 ?stamppFst
 #For Neis Gst use
 pwGst <-stamppNeisD(gl[1:20,]) 
 # threw error message
 # 
 ## Principal Coordinates Analysis (PCoA) 
 
 pc <- gl.pcoa(gl, nfactors=5)
 names(pc)
 
 barplot(pc$eig/sum(pc$eig)*100, ) # percentage of variation that is represented by the axes 
 
 d <- gl.dist.ind(gl, method="euclidean", scale=TRUE) 
 pca <-  gl.pcoa(d) # pca based on distance matrix
 pca
 
 #Plotting the results of PCoA
 #
 gl.pcoa.plot(pc, gl, pop.labels="pop", xaxis=1, yaxis=2)
 #axis 1 explained only 10.5% of the variation and axis explained only 10.3% of the variation
 
 gl.pcoa.plot(pca, gl, pop.labels="pop", xaxis=1, yaxis=2) 
 #axis 1 explained only 9.3% of the variation and axis explained only 8.9% of the variation
 #
 ## Neighbour-joining trees
 ## 
 gl.tree.nj(gl,type = "fan")
 gl.tree.nj(gl,type = "phylogram")
 
 
 
 ##First,load/call the filterd data
 
 library(dartR)
 gl <- gl.load("afrocarpus2.Rdata")  #using 65% call rate 
 
 
 ## Genetic distance
 #a distance matrix can be or  constructed between individuals based on the genetic profiles or between populations based on their allele frequency profiles. 
 #Distance matrices can be generated by a number of R packages, the most popular of which are dist() from package stats and vegdist from package vegan. The function gl.dist is a wrapper for those two functions, applying them to allele frequencies calculated for each locus at each population defined in the genlight object.
 library (stats)
 library (vegan)
 
 gl.dist.ind(gl[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci
 
 dist_pop <- gl.dist.pop(gl) # euclidean distance between populations
 
 ## F Statistics
 #For Fst and Neis Gst, the functions of the StAMPP package allow for the use of parallel computing, which is much faster.
 #It also works on genlight objects directly and therefore no conversion is necessary
 #stamppFst calculates pairwise Fst values between populations
 #
 library(StAMPP)
 pwfst <-stamppFst(gl[1:20,], nboots=1, percent=95, nclusters=1)
 ?stamppFst
 #For Neis Gst use
 pwGst <-stamppNeisD(gl[1:20,]) 
 # threw error message
 # 
 ## Principal Coordinates Analysis (PCoA) 
 
 pc <- gl.pcoa(gl, nfactors=5)
 names(pc)
 
 barplot(pc$eig/sum(pc$eig)*100, ) # percentage of variation that is represented by the axes 
 
 d <- gl.dist.ind(gl, method="euclidean", scale=TRUE) 
 pca <-  gl.pcoa(d) # pca based on distance matrix
 pca
 
 #Plotting the results of PCoA
 #
 gl.pcoa.plot(pc, gl, pop.labels="pop", xaxis=1, yaxis=2)
 #axis 1 explained only 2.3% of the variation and axis explained only 1.5% of the variation
 
 gl.pcoa.plot(pca, gl, pop.labels="pop", xaxis=1, yaxis=2) 
 #axis 1 explained only 2.9% of the variation and axis explained only 1.9% of the variation
 #
 ## Neighbour-joining trees
 ## 
 gl.tree.nj(gl,type = "fan")
 gl.tree.nj(gl,type = "phylogram")
 
 
 
 
 ##First,load/call the filterd data
 
 library(dartR)
 gl <- gl.load("afrocarpus3.Rdata")  #using 70% call rate, with imputation
 
 
 ## Genetic distance
 #a distance matrix can be or  constructed between individuals based on the genetic profiles or between populations based on their allele frequency profiles. 
 #Distance matrices can be generated by a number of R packages, the most popular of which are dist() from package stats and vegdist from package vegan. The function gl.dist is a wrapper for those two functions, applying them to allele frequencies calculated for each locus at each population defined in the genlight object.
 library (stats)
 library (vegan)
 
 gl.dist.ind(gl[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci
 
 dist_pop <- gl.dist.pop(gl) # euclidean distance between populations
 
 ## F Statistics
 #For Fst and Neis Gst, the functions of the StAMPP package allow for the use of parallel computing, which is much faster.
 #It also works on genlight objects directly and therefore no conversion is necessary
 #stamppFst calculates pairwise Fst values between populations
 #
 library(StAMPP)
 pwfst <-stamppFst(gl[1:20,], nboots=1, percent=95, nclusters=1)
 ?stamppFst
 #For Neis Gst use
 pwGst <-stamppNeisD(gl[1:20,]) 
 # threw error message
 # 
 ## Principal Coordinates Analysis (PCoA) 
 
 pc <- gl.pcoa(gl, nfactors=5)
 names(pc)
 
 barplot(pc$eig/sum(pc$eig)*100, ) # percentage of variation that is represented by the axes 
 
 d <- gl.dist.ind(gl, method="euclidean", scale=TRUE) 
 pca <-  gl.pcoa(d) # pca based on distance matrix
 pca
 
 #Plotting the results of PCoA
 #
 gl.pcoa.plot(pc, gl, pop.labels="pop", xaxis=1, yaxis=2)
 #axis 1 explained only 2.9% of the variation and axis explained only 2.2% of the variation
 
 gl.pcoa.plot(pca, gl, pop.labels="pop", xaxis=1, yaxis=2) 
 #axis 1 explained only 2.9% of the variation and axis explained only 2.2% of the variation
 #
 ## Neighbour-joining trees
 ## 
 gl.tree.nj(gl,type = "fan")
 gl.tree.nj(gl,type = "phylogram")

 
 ##First,load/call the filterd data
 
 library(dartR)
 gl <- gl.load("afrocarpus4.Rdata")  #using 65% call rate, with imputation
 
 
 ## Genetic distance
 #a distance matrix can be or  constructed between individuals based on the genetic profiles or between populations based on their allele frequency profiles. 
 #Distance matrices can be generated by a number of R packages, the most popular of which are dist() from package stats and vegdist from package vegan. The function gl.dist is a wrapper for those two functions, applying them to allele frequencies calculated for each locus at each population defined in the genlight object.
 library (stats)
 library (vegan)
 
 gl.dist.ind(gl[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci
 
 dist_pop <- gl.dist.pop(gl) # euclidean distance between populations
 
 ## F Statistics
 #For Fst and Neis Gst, the functions of the StAMPP package allow for the use of parallel computing, which is much faster.
 #It also works on genlight objects directly and therefore no conversion is necessary
 #stamppFst calculates pairwise Fst values between populations
 #
 library(StAMPP)
 pwfst <-stamppFst(gl[1:20,], nboots=1, percent=95, nclusters=1)
 ?stamppFst
 #For Neis Gst use
 pwGst <-stamppNeisD(gl[1:20,]) 
 # threw error message
 # 
 ## Principal Coordinates Analysis (PCoA) 
 
 pc <- gl.pcoa(gl, nfactors=5)
 names(pc)
 
 barplot(pc$eig/sum(pc$eig)*100, ) # percentage of variation that is represented by the axes 
 
 d <- gl.dist.ind(gl, method="euclidean", scale=TRUE) 
 pca <-  gl.pcoa(d) # pca based on distance matrix
 pca
 
 #Plotting the results of PCoA
 #
 gl.pcoa.plot(pc, gl, pop.labels="pop", xaxis=1, yaxis=2)
 #axis 1 explained only 3.4% of the variation and axis explained only 2% of the variation
 
 gl.pcoa.plot(pca, gl, pop.labels="pop", xaxis=1, yaxis=2) 
 #axis 1 explained only 3.4% of the variation and axis explained only 2% of the variation
 #
 ## Neighbour-joining trees
 ## 
 gl.tree.nj(gl,type = "fan")
 gl.tree.nj(gl,type = "phylogram")
  
 
 ##First,load/call the filtered data
 
 library(dartR)
 gl <- gl.load("afrocarpus7.Rdata")  #specific sites used as 'pop' while preparing the individual meta file, filtering individuals with call rate < 70 %call rate and loci with call rate < 70%
 
 
 ## Genetic distance
 #a distance matrix can be or  constructed between individuals based on the genetic profiles or between populations based on their allele frequency profiles. 
 #Distance matrices can be generated by a number of R packages, the most popular of which are dist() from package stats and vegdist from package vegan. The function gl.dist is a wrapper for those two functions, applying them to allele frequencies calculated for each locus at each population defined in the genlight object.
 
 library (stats)
 library (vegan)
 
 gl.dist.ind(gl[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci
 
 dist_pop <- gl.dist.pop(gl) # euclidean distance between populations
 
 ## F Statistics
 #For Fst and Neis Gst, the functions of the StAMPP package allow for the use of parallel computing, which is much faster.
 #It also works on genlight objects directly and therefore no conversion is necessary
 #stamppFst calculates pairwise Fst values between populations
 #
 library(StAMPP)
 pwfst <-stamppFst(gl[1:20,], nboots=1, percent=95, nclusters=1)
 ?stamppFst
 #For Neis Gst use
 pwGst <-stamppNeisD(gl[1:20,]) 
 # threw error message
 # 
 ## Principal Coordinates Analysis (PCoA) 
 
 pc <- gl.pcoa(gl, nfactors=5)
 names(pc)
 
 barplot(pc$eig/sum(pc$eig)*100, ) # percentage of variation that is represented by the axes 
 
 d <- gl.dist.ind(gl, method="euclidean", scale=TRUE) 
 pca <-  gl.pcoa(d) # pca based on distance matrix
 pca
 
 #Plotting the results of PCoA
 #
 gl.pcoa.plot(pc, gl, pop.labels="pop", xaxis=1, yaxis=2)
 #axis 1 explained only 2.4% of the variation and axis explained only 1.7% of the variation
 
 gl.pcoa.plot(pca, gl, pop.labels="pop", xaxis=1, yaxis=2) 
 
 ## Neighbour-joining trees
 ## 
 gl.tree.nj(gl,type = "fan")
 
 gl.tree.nj(gl,type = "phylogram")
 
 
 
 
 ##First,load/call the filtered data
 
 library(dartR)
 gl <- gl.load("afrocarpus8.Rdata")  #specific sites used as 'pop' while preparing the individual meta file, 90 %call rate
 
 gl_progeny <- gl[gl$other$ind.metrics$Stage!= "parent", ] # extracts and creates a new genlight object containing individuals of the 'progeny' growth stage only
 gl_parent <- gl[gl$other$ind.metrics$Stage!= "progeny", ] # extracts and creates a new genlight object containing individuals of the 'parent' growth stage only
 
 
 ## Genetic distance
 #a distance matrix can be or  constructed between individuals based on the genetic profiles or between populations based on their allele frequency profiles. 
 #Distance matrices can be generated by a number of R packages, the most popular of which are dist() from package stats and vegdist from package vegan. The function gl.dist is a wrapper for those two functions, applying them to allele frequencies calculated for each locus at each population defined in the genlight object.
 
 library (stats)
 library (vegan)
 
 gl.dist.ind(gl[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci  # for the whole data
 gl.dist.ind(gl_parent[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci  # for the parents data
 gl.dist.ind(gl_progeny[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci  # for the progeny data
 
 dist_pop <- gl.dist.pop(gl) # euclidean distance between populations for the whole data
 dist_pop2 <- gl.dist.pop(gl_parent) # euclidean distance between populations for the parents data
 dist_pop3 <- gl.dist.pop(gl_progeny) # euclidean distance between populations for the progeny data
 
 ## F Statistics
 #For Fst and Neis Gst, the functions of the StAMPP package allow for the use of parallel computing, which is much faster.
 #It also works on genlight objects directly and therefore no conversion is necessary
 #stamppFst calculates pairwise Fst values between populations
 #
 library(StAMPP)
 pwfst <-stamppFst(gl[1:20,], nboots=1, percent=95, nclusters=1)
 ?stamppFst
 #For Neis Gst use
 pwGst <-stamppNeisD(gl[1:20,]) 
 # threw error message
 # 
 ## Principal Coordinates Analysis (PCoA) 
 
 pc <- gl.pcoa(gl, nfactors=5) # for the whole data
 names(pc)
 pc2 <- gl.pcoa(gl_parent, nfactors=5) # for the parent data
 pc3 <- gl.pcoa(gl_progeny, nfactors=5) # for the progeny data
 
 barplot(pc$eig/sum(pc$eig)*100, ) # percentage of variation that is represented by the axes for the whole data 
 barplot(pc2$eig/sum(pc2$eig)*100, ) # for the parent data
 barplot(pc3$eig/sum(pc3$eig)*100, ) # for the progeny data
 
 
 d <- gl.dist.ind(gl, method="euclidean", scale=TRUE) #distance for the whole data
 d2 <- gl.dist.ind(gl_parent, method="euclidean", scale=TRUE) #distance for the parent data
 d3 <- gl.dist.ind(gl_progeny, method="euclidean", scale=TRUE) #distance for the progeny data
 pca <-  gl.pcoa(d) # pca based on distance matrix for whole data
 pca2 <- gl.pcoa(d2) # pca based on distance matrix for parent data
 pca3 <- gl.pcoa(d3) # pca based on distance matrix for progeny data
 pca
 
 #Plotting the results of PCoA
 #
 gl.pcoa.plot(pc, gl, pop.labels="pop", xaxis=1, yaxis=2) # whole data
 #axis 1 explained only 6% of the variation and axis explained only 4.9 % of the variation
 
 gl.pcoa.plot(pc2, gl_parent, pop.labels="pop", xaxis=1, yaxis=2) # parent data
 #axis 1 explained only 9.5% of the variation and axis explained only 7.8% of the variation
 
 gl.pcoa.plot(pc3, gl_progeny, pop.labels="pop", xaxis=1, yaxis=2) # whole data
 #axis 1 explained only 6.5% of the variation and axis explained only 5.2% of the variation
 #
 gl.pcoa.plot(pca, gl, pop.labels="pop", xaxis=1, yaxis=2) # whole data
 #axis 1 explained only 5.9% of the variation and axis explained only 4.8% of the variation
 #
 gl.pcoa.plot(pca2, gl_parent, pop.labels="pop", xaxis=1, yaxis=2) #parent
 #axis 1 explained only 10.1% of the variation and axis explained only 8.2% of the variation
 gl.pcoa.plot(pca3, gl_progeny, pop.labels="pop", xaxis=1, yaxis=2) #progeny
 #axis 1 explained only 6.9% of the variation and axis explained only 5.7% of the variation
 #
 
 ## Neighbour-joining trees
 ## 
 gl.tree.nj(gl,type = "fan") #whole data
 gl.tree.nj(gl_parent,type = "fan") #parent
 gl.tree.nj(gl_progeny,type = "fan") #progeny
 
 gl.tree.nj(gl,type = "phylogram") #whole data
 gl.tree.nj(gl_parent,type = "phylogram") #parent
 gl.tree.nj(gl_progeny,type = "phylogram") #progeny
 

 ##First,load/call the filtered data
 
 library(dartR)
 gl <- gl.load("afrocarpus7.Rdata")  #specific sites used as 'pop' while preparing the individual meta file, filtering individuals with call rate < 70 %call rate and loci with call rate < 70%
 
 
 ## Genetic distance
 #a distance matrix can be or  constructed between individuals based on the genetic profiles or between populations based on their allele frequency profiles. 
 #Distance matrices can be generated by a number of R packages, the most popular of which are dist() from package stats and vegdist from package vegan. The function gl.dist is a wrapper for those two functions, applying them to allele frequencies calculated for each locus at each population defined in the genlight object.
 
 library (stats)
 library (vegan)
 
 gl.dist.ind(gl[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci
 
 dist_pop <- gl.dist.pop(gl) # euclidean distance between populations
 
 ## F Statistics
 #For Fst and Neis Gst, the functions of the StAMPP package allow for the use of parallel computing, which is much faster.
 #It also works on genlight objects directly and therefore no conversion is necessary
 #stamppFst calculates pairwise Fst values between populations
 #
 library(StAMPP)
 pwfst <-stamppFst(gl[1:20,], nboots=1, percent=95, nclusters=1)
 ?stamppFst
 #For Neis Gst use
 pwGst <-stamppNeisD(gl[1:20,]) 
 # threw error message
 # 
 ## Principal Coordinates Analysis (PCoA) 
 
 pc <- gl.pcoa(gl, nfactors=5)
 names(pc)
 
 barplot(pc$eig/sum(pc$eig)*100, ) # percentage of variation that is represented by the axes 
 
 d <- gl.dist.ind(gl, method="euclidean", scale=TRUE) 
 pca <-  gl.pcoa(d) # pca based on distance matrix
 pca
 
 #Plotting the results of PCoA
 #
 gl.pcoa.plot(pc, gl, pop.labels="pop", xaxis=1, yaxis=2)
 #axis 1 explained only 2.4% of the variation and axis explained only 1.7% of the variation
 
 gl.pcoa.plot(pca, gl, pop.labels="pop", xaxis=1, yaxis=2) 
 
 ## Neighbour-joining trees
 ## 
 gl.tree.nj(gl,type = "fan")
 
 gl.tree.nj(gl,type = "phylogram")
 
 
 
 
 ##First,load/call the filtered data
 
 library(dartR)
 gl <- gl.load("afrocarpus8.Rdata")  #specific sites used as 'pop' while preparing the individual meta file, 90 %call rate, 157 SNPs
 
 gl_progeny <- gl[gl$other$ind.metrics$Stage!= "parent", ] # extracts and creates a new genlight object containing individuals of the 'progeny' growth stage only
 gl_parent <- gl[gl$other$ind.metrics$Stage!= "progeny", ] # extracts and creates a new genlight object containing individuals of the 'parent' growth stage only
 
 
 ## Genetic distance
 #a distance matrix can be or  constructed between individuals based on the genetic profiles or between populations based on their allele frequency profiles. 
 #Distance matrices can be generated by a number of R packages, the most popular of which are dist() from package stats and vegdist from package vegan. The function gl.dist is a wrapper for those two functions, applying them to allele frequencies calculated for each locus at each population defined in the genlight object.
 
 library (stats)
 library (vegan)
 
 gl.dist.ind(gl[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci  # for the whole data
 gl.dist.ind(gl_parent[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci  # for the parents data
 gl.dist.ind(gl_progeny[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci  # for the progeny data
 
 dist_pop <- gl.dist.pop(gl) # euclidean distance between populations for the whole data
 dist_pop2 <- gl.dist.pop(gl_parent) # euclidean distance between populations for the parents data
 dist_pop3 <- gl.dist.pop(gl_progeny) # euclidean distance between populations for the progeny data
 
 ## F Statistics
 #For Fst and Neis Gst, the functions of the StAMPP package allow for the use of parallel computing, which is much faster.
 #It also works on genlight objects directly and therefore no conversion is necessary
 #stamppFst calculates pairwise Fst values between populations
 #
 library(StAMPP)
 pwfst <-stamppFst(gl[1:20,], nboots=1, percent=95, nclusters=1)
 ?stamppFst
 #For Neis Gst use
 pwGst <-stamppNeisD(gl[1:20,]) 
 # threw error message
 # 
 ## Principal Coordinates Analysis (PCoA) 
 
 pc <- gl.pcoa(gl, nfactors=5) # for the whole data
 names(pc)
 pc2 <- gl.pcoa(gl_parent, nfactors=5) # for the parent data
 pc3 <- gl.pcoa(gl_progeny, nfactors=5) # for the progeny data
 
 barplot(pc$eig/sum(pc$eig)*100, ) # percentage of variation that is represented by the axes for the whole data 
 barplot(pc2$eig/sum(pc2$eig)*100, ) # for the parent data
 barplot(pc3$eig/sum(pc3$eig)*100, ) # for the progeny data
 
 
 d <- gl.dist.ind(gl, method="euclidean", scale=TRUE) #distance for the whole data
 d2 <- gl.dist.ind(gl_parent, method="euclidean", scale=TRUE) #distance for the parent data
 d3 <- gl.dist.ind(gl_progeny, method="euclidean", scale=TRUE) #distance for the progeny data
 pca <-  gl.pcoa(d) # pca based on distance matrix for whole data
 pca2 <- gl.pcoa(d2) # pca based on distance matrix for parent data
 pca3 <- gl.pcoa(d3) # pca based on distance matrix for progeny data
 pca
 
 #Plotting the results of PCoA
 #
 gl.pcoa.plot(pc, gl, pop.labels="pop", xaxis=1, yaxis=2) # whole data
 #axis 1 explained only 6% of the variation and axis explained only 4.9 % of the variation
 
 gl.pcoa.plot(pc2, gl_parent, pop.labels="pop", xaxis=1, yaxis=2) # parent data
 #axis 1 explained only 9.5% of the variation and axis explained only 7.8% of the variation
 
 gl.pcoa.plot(pc3, gl_progeny, pop.labels="pop", xaxis=1, yaxis=2) # whole data
 #axis 1 explained only 6.5% of the variation and axis explained only 5.2% of the variation
 #
 gl.pcoa.plot(pca, gl, pop.labels="pop", xaxis=1, yaxis=2) # whole data
 #axis 1 explained only 5.9% of the variation and axis explained only 4.8% of the variation
 #
 gl.pcoa.plot(pca2, gl_parent, pop.labels="pop", xaxis=1, yaxis=2) #parent
 #axis 1 explained only 10.1% of the variation and axis explained only 8.2% of the variation
 gl.pcoa.plot(pca3, gl_progeny, pop.labels="pop", xaxis=1, yaxis=2) #progeny
 #axis 1 explained only 6.9% of the variation and axis explained only 5.7% of the variation
 #
 
 ## Neighbour-joining trees
 ## 
 gl.tree.nj(gl,type = "fan") #whole data
 gl.tree.nj(gl_parent,type = "fan") #parent
 gl.tree.nj(gl_progeny,type = "fan") #progeny
 
 gl.tree.nj(gl,type = "phylogram") #whole data
 gl.tree.nj(gl_parent,type = "phylogram") #parent
 gl.tree.nj(gl_progeny,type = "phylogram") #progeny
 
 
 
 ##First,load/call the filtered data
 
 library(dartR)
 gl <- gl.load("afrocarpus9.Rdata")  #specific sites used as 'pop' while preparing the individual meta file, 80% loci call rate, 1017 SNPs
 
 gl_progeny <- gl[gl$other$ind.metrics$Stage!= "parent", ] # extracts and creates a new genlight object containing individuals of the 'progeny' growth stage only
 gl_parent <- gl[gl$other$ind.metrics$Stage!= "progeny", ] # extracts and creates a new genlight object containing individuals of the 'parent' growth stage only
 
 
 ## Genetic distance
 #a distance matrix can be or  constructed between individuals based on the genetic profiles or between populations based on their allele frequency profiles. 
 #Distance matrices can be generated by a number of R packages, the most popular of which are dist() from package stats and vegdist from package vegan. The function gl.dist is a wrapper for those two functions, applying them to allele frequencies calculated for each locus at each population defined in the genlight object.
 
 library (stats)
 library (vegan)
 
 gl.dist.ind(gl[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci  # for the whole data
 gl.dist.ind(gl_parent[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci  # for the parents data
 gl.dist.ind(gl_progeny[1:7, 1:100], method="euclidean") # euclidean distance for the first 7 individuals and first 100 loci  # for the progeny data
 
 dist_pop <- gl.dist.pop(gl) # euclidean distance between populations for the whole data
 dist_pop2 <- gl.dist.pop(gl_parent) # euclidean distance between populations for the parents data
 dist_pop3 <- gl.dist.pop(gl_progeny) # euclidean distance between populations for the progeny data
 
 ## F Statistics
 #For Fst and Neis Gst, the functions of the StAMPP package allow for the use of parallel computing, which is much faster.
 #It also works on genlight objects directly and therefore no conversion is necessary
 #stamppFst calculates pairwise Fst values between populations
 #
 library(StAMPP)
 pwfst <-stamppFst(gl[1:20,], nboots=1, percent=95, nclusters=1)
 ?stamppFst
 #For Neis Gst use
 pwGst <-stamppNeisD(gl[1:20,]) 
 # threw error message
 # 
 ## Principal Coordinates Analysis (PCoA) 
 
 pc <- gl.pcoa(gl, nfactors=5) # for the whole data
 names(pc)
 pc2 <- gl.pcoa(gl_parent, nfactors=5) # for the parent data
 pc3 <- gl.pcoa(gl_progeny, nfactors=5) # for the progeny data
 
 barplot(pc$eig/sum(pc$eig)*100, ) # percentage of variation that is represented by the axes for the whole data 
 barplot(pc2$eig/sum(pc2$eig)*100, ) # for the parent data
 barplot(pc3$eig/sum(pc3$eig)*100, ) # for the progeny data
 
 
 d <- gl.dist.ind(gl, method="euclidean", scale=TRUE) #distance for the whole data
 d2 <- gl.dist.ind(gl_parent, method="euclidean", scale=TRUE) #distance for the parent data
 d3 <- gl.dist.ind(gl_progeny, method="euclidean", scale=TRUE) #distance for the progeny data
 pca <-  gl.pcoa(d) # pca based on distance matrix for whole data
 pca2 <- gl.pcoa(d2) # pca based on distance matrix for parent data
 pca3 <- gl.pcoa(d3) # pca based on distance matrix for progeny data
 pca
 
 #Plotting the results of PCoA
 #
 gl.pcoa.plot(pc, gl, pop.labels="pop", xaxis=1, yaxis=2) # whole data
 #axis 1 explained only 2.9% of the variation and axis explained only 2 % of the variation
 
 gl.pcoa.plot(pc2, gl_parent, pop.labels="pop", xaxis=1, yaxis=2) # parent data
 #axis 1 explained only 3.9% of the variation and axis explained only 3.7% of the variation
 
 gl.pcoa.plot(pc3, gl_progeny, pop.labels="pop", xaxis=1, yaxis=2) # whole data
 #axis 1 explained only 4.1% of the variation and axis explained only 2.8% of the variation
 #
 gl.pcoa.plot(pca, gl, pop.labels="pop", xaxis=1, yaxis=2) # whole data
 #axis 1 explained only 3.2% of the variation and axis explained only 2.4% of the variation
 #
 gl.pcoa.plot(pca2, gl_parent, pop.labels="pop", xaxis=1, yaxis=2) #parent
 #axis 1 explained only 5.1% of the variation and axis explained only 4.6% of the variation
 gl.pcoa.plot(pca3, gl_progeny, pop.labels="pop", xaxis=1, yaxis=2) #progeny
 #axis 1 explained only 4.5% of the variation and axis explained only 3.3% of the variation
 #
 
 ## Neighbour-joining trees
 ## 
 gl.tree.nj(gl,type = "fan") #whole data
 gl.tree.nj(gl_parent,type = "fan") #parent
 gl.tree.nj(gl_progeny,type = "fan") #progeny
 
 gl.tree.nj(gl,type = "phylogram") #whole data
 gl.tree.nj(gl_parent,type = "phylogram") #parent
 gl.tree.nj(gl_progeny,type = "phylogram") #progeny
 