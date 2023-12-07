
### Working on the DArT data

## Reading DArT Files into a genlight Object using read.dart()

### first, install the dartR package
install.packages("devtools")
library(devtools)
install_github("green-striped-gecko/dartR")
library(dartR)
gl.install.vanilla.dartR() 

#Reading DArT Files into a genlight Object using read.dart()

gl <- gl.read.dart(filename = "Report_DAfg23-2937_SNP_2.csv", ind.metafile = "metadata_ind.csv")
gl
nLoc(gl) # 10219 initial SNPS


##1. Filtering (95% reproducibility, >70% call rate,  monomorphic loci, and removing individuals with < 50% call rate )
gl2 <- gl.filter.secondaries(gl) # filter out SNPs that share a sequence tag, except one retained at random

nLoc(gl2) #produced 9819 SNPs
gl3 <- gl.filter.rdepth(gl2) # filter out loci with exceptionally low or high read depth (coverage)
nLoc(gl3) #produced 8551 SNPs

gl4 <- gl.filter.reproducibility(gl3, threshold = 0.95) # filter out loci for which the reproducibility (strictly repeatability) is less than threshold = 0.95

gl4 <- gl.recalc.metrics(gl4) # recalculating the locus metrics after filtering because the initial call rate parameter will no longer be accurate after the filtering out individuals or populations
nLoc(gl4) #produced 7519 SNPs

gl5 <- gl.filter.callrate(gl4, method = "ind", threshold = 0.50) #Filter individuals on call rate (threshold =50% )deleted individuals. Individuals deleted (CallRate <=  0.5 ):AG025[Tesso/Borite], AG041[Wonsho Abbo], AG057[Wonsho Gudumale], AG002[Tesso/Sodicho], AG090[Wondogenet college], AG026[Tesso/Borite], AG042[Wonsho Abbo], AG050[Wonsho Abbo], AG058[Wonsho Gudumale], AG083[Wondogenet college], AG027[Tesso/Borite], AG043[Wonsho Abbo], AG060[Wonsho Gudumale], AG029[Tesso/Borite], AG037[Wonsho Abbo], AG053[Wonsho Abbo], AG022[Tesso/Borite], AG030[Tesso/Borite], AG038[Wonsho Abbo], AG079[Wondogenet college], AG048[Wonsho Abbo], AG064[Wonsho Gudumale], AG103[Auger], AG127[ME Gedam], AG135[ME Gedam], AG143[Deban/Amba den], AG151[Deban/Amba den], AG128[ME Gedam], AG136[ME Gedam], AG169[Bera-Tedicho], AG137[ME Gedam], AG106[Auger], AG138[ME Gedam], AG154[Deban/Amba den], AG131[ME Gedam], AG155[Deban/Amba den], AG172[Bera-Tedicho], AG140[ME Gedam], AG148[Deban/Amba den], AG173[Bera-Tedicho], AG181[Bera-Tedicho], AG109[Auger], AG141[Deban/Amba den], AG149[Deban/Amba den], AG174[Bera-Tedicho], AG142[Deban/Amba den], AG150[Deban/Amba den],

gl5 <- gl.recalc.metrics(gl5)
nLoc(gl5) # produced 7519 SNPs

gl6 <- gl.filter.callrate(gl5, method="loc", threshold = 0.70)  #filter out loci with callrate lower than 70%

gl6 <- gl.recalc.metrics(gl6)

nLoc(gl6) # produced 3134 SNPs

gl7 <- gl.filter.monomorphs(gl6, v=0) #filter out monomorphic loci
gl7 <- gl.recalc.metrics(gl7)
gl7
nLoc(gl7) # produced 2745 SNPs

m <- as.matrix(gl7) #save genlight object as matrix
write.csv(m,file="filtered_loci_metadata1_1.csv") # saves as .csv
gl.save(gl7,file="afrocarpus1_1.Rdata") # saves the the data as an R object filtered with 95% reproducibility, 70 % call rate, deleted individuals with less than 50% call rate, and monomorphic loci filtered out. it resulted in 2745 SNPS 


##2. Filtering (95% reproducibility, >65% call rate,  monomorphic loci, and removing individuals with < 50% call rate )
gl2 <- gl.filter.secondaries(gl) # filter out SNPs that share a sequence tag, except one retained at random

nLoc(gl2) #produced 9819 SNPs
gl3 <- gl.filter.rdepth(gl2) # filter out loci with exceptionally low or high read depth (coverage)
nLoc(gl3) #produced 8551 SNPs

gl4 <- gl.filter.reproducibility(gl3, threshold = 0.95) # filter out loci for which the reproducibility (strictly repeatability) is less than threshold = 0.95

nLoc(gl4) #produced 7519 SNPs
gl4 <- gl.recalc.metrics(gl4) # recalculating the locus metrics after filtering because the initial call rate parameter will no longer be accurate after the filtering out individuals or populations

gl5 <- gl.filter.callrate(gl4, method = "loc", threshold = 0.65) #filter out loci with callrate lower than 70%
gl5 <- gl.recalc.metrics(gl5)
nLoc(gl5) # produced 2745 SNPs

gl6 <- gl.filter.callrate(gl5, method="ind", threshold = 0.50) #Filter individuals on call rate (threshold =50% )deleted individuals Individuals deleted (CallRate <=  0.5 ): AG041[Wonsho Abbo], AG083[Wondogenet college], AG022[Tesso/Borite], AG128[ME Gedam], AG137[ME Gedam], AG172[Bera-Tedicho], AG150[Deban/Amba den]

gl6 <- gl.recalc.metrics(gl6)

nLoc(gl6) # produced 2745 SNPs

gl7 <- gl.filter.monomorphs(gl6, v=0) #filter out monomorphic loci
gl7 <- gl.recalc.metrics(gl7)
gl7
nLoc(gl7) # produced 2675 SNPs

m <- as.matrix(gl7) #save genlight object as matrix
write.csv(m,file="filtered_loci_metadata2.csv") # saves as .csv
gl.save(gl7,file="afrocarpus2.Rdata") # saves the the data as an R object filtered with 95% reproducibility, 65 % call rate, deleted individuals with less than 50% call rate, and monomorphic loci filtered out. it resulted in 2675 SNPS 


##3. Filtering (95% reproducibility, >70% call rate,  monomorphic loci, and without removing individuals with  low call rate, rather use imputation )
gl2 <- gl.filter.secondaries(gl) # filter out SNPs that share a sequence tag, except one retained at random

nLoc(gl2) #produced 9819 SNPs
gl3 <- gl.filter.rdepth(gl2) # filter out loci with exceptionally low or high read depth (coverage)
nLoc(gl3) #produced 8551 SNPs

gl4 <- gl.filter.reproducibility(gl3, threshold = 0.95) # filter out loci for which the reproducibility (strictly repeatability) is less than threshold = 0.95

nLoc(gl4) #produced 7519 SNPs
gl4 <- gl.recalc.metrics(gl4) # recalculating the locus metrics after filtering because the initial call rate parameter will no longer be accurate after the filtering out individuals or populations

gl5 <- gl.filter.callrate(gl4, method = "loc", threshold = 0.7) #filter out loci with callrate lower than 70%
gl5 <- gl.recalc.metrics(gl5)
nLoc(gl5) # produced 1734 SNPs

gl6 <- gl.filter.monomorphs(gl5, v=0) #filter out monomorphic loci
gl6 <- gl.recalc.metrics(gl6)
gl6
nLoc(gl6) # produced 1734 SNPs

gl7 <- gl.filter.allna(gl6) #filter out loci that are all missing values (NA)
gl7 <- gl.impute(gl7, method = "neighbour")
gl7 <- gl.recalc.metrics(gl7)

nLoc(gl7) # produced 1734 SNPs

m <- as.matrix(gl7) #save genlight object as matrix
write.csv(m,file="filtered_loci_metadata3.csv") # saves as .csv
gl.save(gl7,file="afrocarpus3.Rdata") # saves the the data as an R object filtered with 95% reproducibility, 70 % call rate, used imputation rather than deleting individuals with lower call rate, and monomorphic loci filtered out. it resulted in 1734 SNPS 


##4. Filtering (95% reproducibility, >65% call rate,  monomorphic loci, and without removing individuals with  low call rate, rather use imputation )
gl2 <- gl.filter.secondaries(gl) # filter out SNPs that share a sequence tag, except one retained at random

nLoc(gl2) #produced 9819 SNPs
gl3 <- gl.filter.rdepth(gl2) # filter out loci with exceptionally low or high read depth (coverage)
nLoc(gl3) #produced 8551 SNPs

gl4 <- gl.filter.reproducibility(gl3, threshold = 0.95) # filter out loci for which the reproducibility (strictly repeatability) is less than threshold = 0.95

nLoc(gl4) #produced 7519 SNPs
gl4 <- gl.recalc.metrics(gl4) # recalculating the locus metrics after filtering because the initial call rate parameter will no longer be accurate after the filtering out individuals or populations

gl5 <- gl.filter.callrate(gl4, method = "loc", threshold = 0.65) #filter out loci with callrate lower than 65%
gl5 <- gl.recalc.metrics(gl5)
nLoc(gl5) # produced 2745 SNPs

gl6 <- gl.filter.monomorphs(gl5, v=0) #filter out monomorphic loci
gl6 <- gl.recalc.metrics(gl6)
gl6
nLoc(gl6) # produced 2745 SNPs

gl7 <- gl.filter.allna(gl6) #filter out loci that are all missing values (NA)
gl7 <- gl.impute(gl7, method = "neighbour")
gl7 <- gl.recalc.metrics(gl7)

nLoc(gl7) # produced 2745 SNPs

gl7 <- gl.filter.monomorphs(gl6, v=0) #filter out monomorphic loci
gl7 <- gl.recalc.metrics(gl7)
gl7
nLoc(gl7) # produced 2745 SNPs

m <- as.matrix(gl7) #save genlight object as matrix
write.csv(m,file="filtered_loci_metadata4.csv") # saves as .csv
gl.save(gl7,file="afrocarpus4.Rdata") # saves the the data as an R object filtered with 95% reproducibility, 650 % call rate, used imputation rather than deleting individuals with lower call rate, and monomorphic loci filtered out. it resulted in 2745 SNPS 


##5. Filtering (95% reproducibility, >80% call rate,  monomorphic loci, and removing individuals with < 50% call rate )
gl2 <- gl.filter.secondaries(gl) # filter out SNPs that share a sequence tag, except one retained at random

nLoc(gl2) #produced 9819 SNPs
gl3 <- gl.filter.rdepth(gl2) # filter out loci with exceptionally low or high read depth (coverage)
nLoc(gl3) #produced 8548 SNPs

gl4 <- gl.filter.reproducibility(gl3, threshold = 0.95) # filter out loci for which the reproducibility (strictly repeatability) is less than threshold = 0.95

nLoc(gl4) #produced 7521 SNPs
gl4 <- gl.recalc.metrics(gl4) # recalculating the locus metrics after filtering because the initial call rate parameter will no longer be accurate after the filtering out individuals or populations

gl5 <- gl.filter.callrate(gl4, method = "loc", threshold = 0.8) #filter out loci with callrate lower than 80%
gl5 <- gl.recalc.metrics(gl5)
nLoc(gl5) # produced 444 SNPs

gl6 <- gl.filter.callrate(gl5, method="ind", threshold = 0.50) #Filter individuals on call rate (threshold =50% )deleted individuals Individuals deleted (CallRate <=  0.5 ): AG041[Wonsho Abbo], AG083[Wondogenet college], AG022[Tesso/Borite], AG128[ME Gedam], AG137[ME Gedam], AG172[Bera-Tedicho], AG150[Deban/Amba den]

gl6 <- gl.recalc.metrics(gl6)

nLoc(gl6) # produced 1734 SNPs

gl7 <- gl.filter.monomorphs(gl6, v=0) #filter out monomorphic loci
gl7 <- gl.recalc.metrics(gl7)
gl7
nLoc(gl7) # produced 433 SNPs

m <- as.matrix(gl7) #save genlight object as matrix
write.csv(m,file="filtered_loci_metadata5.csv") # saves as .csv
gl.save(gl7,file="afrocarpus5.Rdata") # saves the the data as an R object filtered with 95% reproducibility, 80 % call rate, deleted individuals with less than 50% call rate, and monomorphic loci filtered out. it resulted in 433 SNPS 



##6. Filtering (95% reproducibility, >90% call rate,  monomorphic loci, and removing individuals with < 50% call rate )
gl2 <- gl.filter.secondaries(gl) # filter out SNPs that share a sequence tag, except one retained at random

nLoc(gl2) #produced 9819 SNPs
gl3 <- gl.filter.rdepth(gl2) # filter out loci with exceptionally low or high read depth (coverage)
nLoc(gl3) #produced 8549 SNPs

gl4 <- gl.filter.reproducibility(gl3, threshold = 0.95) # filter out loci for which the reproducibility (strictly repeatability) is less than threshold = 0.95

nLoc(gl4) #produced 7517 SNPs
gl4 <- gl.recalc.metrics(gl4) # recalculating the locus metrics after filtering because the initial call rate parameter will no longer be accurate after the filtering out individuals or populations

gl5 <- gl.filter.callrate(gl4, method = "loc", threshold = 0.9) #filter out loci with callrate lower than 90%
gl5 <- gl.recalc.metrics(gl5)
nLoc(gl5) # produced 50 SNPs

gl6 <- gl.filter.callrate(gl5, method="ind", threshold = 0.50) #Filter individuals on call rate (threshold =50% )deleted individuals Individuals deleted (CallRate <=  0.5 ): AG041[Wonsho Abbo], AG083[Wondogenet college], AG022[Tesso/Borite], AG128[ME Gedam], AG137[ME Gedam], AG172[Bera-Tedicho], AG150[Deban/Amba den]

gl6 <- gl.recalc.metrics(gl6)

nLoc(gl6) # produced 50 SNPs

gl7 <- gl.filter.monomorphs(gl6, v=0) #filter out monomorphic loci
gl7 <- gl.recalc.metrics(gl7)
gl7
nLoc(gl7) # produced 49 SNPs

m <- as.matrix(gl7) #save genlight object as matrix
write.csv(m,file="filtered_loci_metadata6.csv") # saves as .csv
gl.save(gl7,file="afrocarpus6.Rdata") # saves the the data as an R object filtered with 95% reproducibility, 80 % call rate, deleted individuals with less than 50% call rate, and monomorphic loci filtered out. it resulted in 433 SNPS 

### Use another population designation, zone/region

#Reading DArT Files into a genlight Object using read.dart()

gl <- gl.read.dart(filename = "Report_DAfg23-2937_SNP_2.csv", ind.metafile = "metadata_ind_1.csv")
gl
nLoc(gl) # 10219 initial SNPS


##1. Filtering (95% reproducibility, >70% call rate,  monomorphic loci, and removing individuals with < 50% call rate )
gl2 <- gl.filter.secondaries(gl) # filter out SNPs that share a sequence tag, except one retained at random

nLoc(gl2) #produced 9819 SNPs
gl3 <- gl.filter.rdepth(gl2) # filter out loci with exceptionally low or high read depth (coverage)
nLoc(gl3) #produced 8551 SNPs

gl4 <- gl.filter.reproducibility(gl3, threshold = 0.95) # filter out loci for which the reproducibility (strictly repeatability) is less than threshold = 0.95

nLoc(gl4) #produced 7519 SNPs
gl4 <- gl.recalc.metrics(gl4) # recalculating the locus metrics after filtering because the initial call rate parameter will no longer be accurate after the filtering out individuals or populations

gl5 <- gl.filter.callrate(gl4, method = "loc", threshold = 0.7) #filter out loci with callrate lower than 70%
gl5 <- gl.recalc.metrics(gl5)
nLoc(gl5) # produced 1734 SNPs

gl6 <- gl.filter.callrate(gl5, method="ind", threshold = 0.50) #Filter individuals on call rate (threshold =50% )deleted individuals Individuals deleted (CallRate <=  0.5 ): AG041[Wonsho Abbo], AG083[Wondogenet college], AG022[Tesso/Borite], AG128[ME Gedam], AG137[ME Gedam], AG172[Bera-Tedicho], AG150[Deban/Amba den]

gl6 <- gl.recalc.metrics(gl6)

nLoc(gl6) # produced 1734 SNPs

gl7 <- gl.filter.monomorphs(gl6, v=0) #filter out monomorphic loci
gl7 <- gl.recalc.metrics(gl7)
gl7
nLoc(gl7) # produced 1692 SNPs

m <- as.matrix(gl7) #save genlight object as matrix
write.csv(m,file="filtered_loci_metadata5.csv") # saves as .csv
gl.save(gl7,file="afrocarpus5.Rdata") # saves the the data as an R object filtered with 95% reproducibility, 70 % call rate, deleted individuals with less than 50% call rate, and monomorphic loci filtered out. it resulted in 1692 SNPS. I used zone/region as my 'pop' while preparing the individual meta file



##7. Filtering (95% reproducibility, >70% call rate,  monomorphic loci, and removing individuals with < 70% call rate )
gl2 <- gl.filter.secondaries(gl) # filter out SNPs that share a sequence tag, except one retained at random

nLoc(gl2) #produced 9819 SNPs
gl3 <- gl.filter.rdepth(gl2) # filter out loci with exceptionally low or high read depth (coverage)
nLoc(gl3) #produced 8551 SNPs

gl4 <- gl.filter.reproducibility(gl3, threshold = 0.95) # filter out loci for which the reproducibility (strictly repeatability) is less than threshold = 0.95

nLoc(gl4) #produced 7519 SNPs
gl4 <- gl.recalc.metrics(gl4) # recalculating the locus metrics after filtering because the initial call rate parameter will no longer be accurate after the filtering out individuals or populations

gl5 <- gl.filter.callrate(gl4, method = "ind", threshold = 0.7) #filter out individuals with callrate lower than 70%. 
#Filter individuals on call rate (threshold =70% )deleted individuals. Individuals deleted (CallRate <=  0.7 ):AG025[Tesso/Borite], AG041[Wonsho Abbo], AG057[Wonsho Gudumale], AG002[Tesso/Sodicho], AG090[Wondogenet college], AG026[Tesso/Borite], AG042[Wonsho Abbo], AG050[Wonsho Abbo], AG058[Wonsho Gudumale], AG003[Tesso/Sodicho], AG083[Wondogenet college], AG027[Tesso/Borite], AG043[Wonsho Abbo], AG059[Wonsho Gudumale], AG060[Wonsho Gudumale], AG029[Tesso/Borite], AG037[Wonsho Abbo], AG053[Wonsho Abbo], AG022[Tesso/Borite], AG030[Tesso/Borite], AG038[Wonsho Abbo], AG079[Wondogenet college], AG023[Tesso/Borite], AG063[Wonsho Gudumale], AG048[Wonsho Abbo], AG064[Wonsho Gudumale], AG103[Auger], AG127[ME Gedam], AG135[ME Gedam], AG143[Deban/Amba den], AG151[Deban/Amba den], AG128[ME Gedam], AG136[ME Gedam], AG169[Bera-Tedicho], AG137[ME Gedam], AG106[Auger], AG130[ME Gedam], AG138[ME Gedam], AG154[Deban/Amba den], AG171[Bera-Tedicho], AG131[ME Gedam], AG155[Deban/Amba den], AG172[Bera-Tedicho], AG180[Bera-Tedicho], AG140[ME Gedam], AG148[Deban/Amba den], AG173[Bera-Tedicho], AG181[Bera-Tedicho], AG109[Auger], AG141[Deban/Amba den], AG149[Deban/Amba den], AG174[Bera-Tedicho], AG142[Deban/Amba den], AG150[Deban/Amba den],
#
#only 52 individuals left
gl5 <- gl.recalc.metrics(gl5)
nLoc(gl5) # produced 7519 SNPs

gl6 <- gl.filter.callrate(gl5, method="loci", threshold = 0.70) # Filtering out loci with call rate < 70%

gl6 <- gl.recalc.metrics(gl6)

nLoc(gl6) # produced 5550 SNPs

gl7 <- gl.filter.monomorphs(gl6, v=0) #filter out monomorphic loci
gl7 <- gl.recalc.metrics(gl7)
gl7
nLoc(gl7) # produced 3557 SNPs

m <- as.matrix(gl7) #save genlight object as matrix
write.csv(m,file="filtered_loci_metadata7.csv") # saves as .csv
gl.save(gl7,file="afrocarpus7.Rdata") # saves the the data as an R object filtered with 95% reproducibility, deleted individuals with less than 70% call rate, filtered loci with 70 % call rate, , and monomorphic loci filtered out. it resulted in 3557 SNPS 

nLoc(gl7) #returns the number of loci.

nPop(gl7) #returns the number of populations to which the individuals are assigned.

indNames(gl7) #returns or sets labels for individuals.
table(pop(gl7)) # table on individuals per population
barplot(table(pop(gl5)), las=2) # barplot on the number of individuals per population

levels(pop(gl7)) #List the population labels

install.packages("HardyWeinberg")
library(HardyWeinberg)
install.packages("ggtern")
library(ggtern)
gl.report.hwe(gl7) #Reports departure of Hardy-Weinberg-Equilibrium for every loci per
#population or overall


table(pop(gl7)) #table on individuals per population

