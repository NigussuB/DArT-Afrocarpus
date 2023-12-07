####DArTseq data analysis####
install.packages("dartR")
install.packages("devtools")
library(devtools)
install.packages("BiocManager")
BiocManager::install(c("SNPRelate", "qvalue"))
install_github("green-striped-gecko/dartR")
library(dartR)

install.packages("devtools")
library(devtools)
install.packages("BiocManager")
BiocManager::install(c("SNPRelate", "qvalue"))
install_github("green-striped-gecko/dartR")
library(dartR)

#alternative installation (it works!)

install.packages("devtools")
library(devtools)
install_github("green-striped-gecko/dartR")
library(dartR)
gl.install.vanilla.dartR() 

#set working directory

setwd("C:/Users/kalousovam/OneDrive - CZU v Praze/Dokumenty/cacao")

#Reading DArT Files into a genlight Object using read.dart()

gl <- gl.read.dart(filename = "Report_DCa21-5986_SNP_2.csv", ind.metafile = "metadata_ind.csv")
gl

#Filtering (100% reproducibility, >90% call rate, and monomorphic loci)

gl2 <- gl.filter.callrate(gl, method = "loc", threshold = 0.90) #filter out loci with callrate lower than 90%
gl2

gl3 <- gl.filter.callrate(gl2, method="ind", threshold = 0.50) #Filter individuals on call rate (threshold =80% )deleted individuals E17[], E29[], E18[], E11[], E23[], E24[], E15[], E26[]
gl3

gl4 <- gl.filter.reproducibility(gl3, t=1) #filter loci with 100% reproducibility
gl4

gl5 <- gl.filter.monomorphs(gl4, v=0) #filter out monomorphic loci
gl5

m <- as.matrix(gl5) #save genlight object as matrix
write.csv(m,file="filtered_loci_metadata.csv") # saves as .csv

nLoc(gl5) #returns the number of loci.

nPop(gl5) #returns the number of populations to which the individuals are assigned.

indNames(gl5) #returns or sets labels for individuals.

levels(pop(gl5)) #List the population labels

gl.report.hwe(gl5) #Reports departure of Hardy-Weinberg-Equilibrium for every loci per
                   #population or overall

table(pop(gl5)) #table on individuals per population

glG <- gl.drop.pop(gl5, pop.list=c("Ecuador")) #delete population Ecuador (creates gl object with Guatemalan accessions only)

glE <- gl.keep.pop(gl5, pop.list=c("Ecuador")) #create gl object with Ecuadorian accessions only

glE <- gl.recalc.metrics(glE) #When you delete individuals or populations, many of the locus metadata (such as Call Rate) are no longer correct. You can recalculate the locus metadata using the script
glEf <- gl.filter.monomorphs(glE, v=0)                              

########now I have the genlight object ready for further analysis########


####SambaR####
source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.09.txt") #to get the package
getpackages() #to install dependent packages

