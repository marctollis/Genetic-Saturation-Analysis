## Analyzing DNA Saturation 
library(ape)
library(ggplot2)
library(dplyr)
library(tidyverse)
# We are going to use the alignments w complete species sampling
setwd <- ("~/Dropbox/turtleBUSCO_ms/saturation/trimal_subset_allSpecies685/")

# this is the list of the filenames for the alignments
in.files <- list.files("~/Dropbox/turtleBUSCO_ms/saturation/trimal_subset_allSpecies685", full.names = TRUE)
filenames <- list.files("~/Dropbox/turtleBUSCO_ms/saturation/trimal_subset_allSpecies685", full.names = FALSE)

# this code attempts to read the dna for each alignment as a separate vector but I will probably abondon this
#aligns <- lapply(in.files, function(x) {
#  dat <- read.dna(x, format = "fasta", as.character = TRUE, skip = 0)
#})


#aligns <- list.files("~/Dropbox/turtleBUSCO_ms/saturation/trimal_subset_allSpecies685/trimal_subset_allSpecies685", full.names = FALSE)


coeffs<- lapply(in.files, function(x) { # this function takes the list of alignments, 
  dat <- read.dna(x, format = "fasta", as.character = TRUE, skip = 0) #loads them into ape
  dat <- as.DNAbin(dat) # and reads them as dnabins
  dist <- dist.dna(dat, model = "raw") #then estimate the uncorrected distances
  dist.corrected <- dist.dna(dat, model = "TN93") # corrected distances
  lm_coef<-coef(lm(dist~dist.corrected)) # calculates the coefficients between
  as.numeric(lm_coef[2]) # saves coefficients as numeric data
  })

# here i attempt to merge the columns of the two lists and make a dataframe
fn <- as.data.frame(filenames)
#fn1 <- t(fn)
sat <- as.data.frame(coeffs)
cfs <- t(sat)
colnames(cfs) <- c('saturation_cf')
df <- as.data.frame(cfs)
bound <- cbind(fn,df)
colnames(bound) <- c('alignment','saturation_cf')
plot <- ggplot(df, aes(x=saturation_cf)) + geom_histogram() + xlab("Slope (Raw distance ~ TN93 distance)") + ylab("Number of Alignments")
plot



mean(bound$saturation_cf)
unsat <- subset(bound, saturation_cf > mean(bound$saturation_cf))
mean(unsat$saturation_cf)
saturated <- subset(bound, saturation_cf < mean(bound$saturation_cf))
mean(saturated$saturation_cf)

range(bound$saturation_cf)
quantile(bound$saturation_cf, probs = c(0.2,0.8))

# print lists of saturated and unsaturated genes
lapply(unsat$alignment, write, "unsat.txt", append=TRUE, ncolumns=1000)
lapply(saturated$alignment, write, "sat.txt", append=TRUE, ncolumns=1000)

#let's compare saturation plots between AApos1_2 and AApos3
pos1_2 <- read.dna("~/Dropbox/turtleBUSCO_ms/saturation/macse/sat_unsat_AApos1_2.fas", format = "fasta", as.character = TRUE, skip = 0) #loads them into ape
pos1_2 <- as.DNAbin(pos1_2) # and reads them as dnabins
dist1_2 <- dist.dna(pos1_2, model = "raw") #then estimate the uncorrected distances
dist.corrected1_2 <- dist.dna(pos1_2, model = "TN93") # corrected distances
lm_coef1<-coef(lm(dist1_2~dist.corrected1_2)) # calculates the coefficients between
lm_coef1
###Make plot###
plot(dist1_2~dist.corrected1_2, pch=20, col="red", xlab="TrN model distance", ylab="Uncorrected genetic distance", main="First and Second Codon Positions")
abline(0,1, lty=2)
abline(lm(dist1_2~dist.corrected1_2), lwd=3)
#lm_coef<-coef(lm(dist1_2~dist.corrected1_2))
text(0.1,0.05,bquote(y == .(lm_coef1[2])*x))

pos3 <- read.dna("~/Dropbox/turtleBUSCO_ms/saturation/macse/sat_unsat_AApos3.fas", format = "fasta", as.character = TRUE, skip = 0) #loads them into ape
pos3 <- as.DNAbin(pos3) # and reads them as dnabins
dist3 <- dist.dna(pos3, model = "raw") #then estimate the uncorrected distances
dist.corrected3 <- dist.dna(pos3, model = "TN93") # corrected distances
lm_coef3<-coef(lm(dist3~dist.corrected3)) # calculates the coefficients between
lm_coef3
###Make plot###
plot(dist3~dist.corrected3, pch=20, col="red", xlab="TrN model distance", ylab="Uncorrected genetic distance", main="Third Codon Positions")
abline(0,1, lty=2)
abline(lm(dist1_2~dist.corrected1_2), lwd=3)
#lm_coef<-coef(lm(dist1_2~dist.corrected1_2))
text(0.1,0.05,bquote(y == .(lm_coef3[2])*x))