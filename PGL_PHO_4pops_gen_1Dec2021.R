#### https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html

setwd("C:/Users/hwu8/Dropbox/My PC (DESKTOP-OOS7MNJ)/Desktop/PGL_PHO_52/52samples_2spp_4pops_no_outgroup")

library(vcfR)
library(poppr)
library(ape)
#library(wesanderson)
library(adegenet)

#homoeophylla.VCF <- read.vcfR("PGL_PHO.vcf")

## Scanning file to determine attributes.
## File attributes:
##   meta lines: 10
## header_line: 11
## variant count: 50984
## column count: 61
## Meta line 10 read in.
## All meta lines processed.
## gt matrix initialized.
## Character matrix gt created.
## Character matrix gt rows: 50984
## Character matrix gt cols: 61
## skip: 0
## nrows: 50984
## row_num: 0
## Processed variant: 50984
## All variants processed

#homoeophylla.VCF

## ***** Object of Class vcfR *****
## 52 samples
## 1147 CHROMs
## 50,984 variants
## Object size: 30.7 Mb
## 0 percent missing data
## *****        *****         *****

# Modify VCF file - Did this on the server (phoebe), outside of R (see below)
#system("vcftools --vcf homoeophylla_glabra.vcf --max-missing 0.5 --maf 0.05 --out homoeophylla_glabra_out --recode 2>&1", intern = T)

## (base) twidhelm@phoebe:~/P.glabra/homoeophylla_project/PGL_PHO_outfiles$ vcftools --vcf PGL_PHO.vcf --max-missing 0.5 --
## maf 0.05 --out PGL_PHO_52.vcf --recode 2>&1

## VCFtools - 0.1.15
## (C) Adam Auton and Anthony Marcketta 2009

## Parameters as interpreted:
##   --vcf PGL_PHO.vcf
## --maf 0.05
## --max-missing 0.5
## --out PGL_PHO_52.vcf
## --recode

## After filtering, kept 52 out of 52 Individuals
## Outputting VCF file...
## After filtering, kept 11689 out of a possible 50984 Sites
## Run Time = 0.00 seconds

# Read in modified VCF file
homoeophylla_glabra.VCF <- read.vcfR("PGL_PHO_52.recode.vcf")

## Scanning file to determine attributes.
## File attributes:
##   meta lines: 10
## header_line: 11
## variant count: 11689
## column count: 61
## Meta line 10 read in.
## All meta lines processed.
## gt matrix initialized.
## Character matrix gt created.
## Character matrix gt rows: 11689
## Character matrix gt cols: 61
## skip: 0
## nrows: 11689
## row_num: 0
## Processed variant: 11689
## All variants processed

homoeophylla_glabra.VCF

## ***** Object of Class vcfR *****
## 52 samples
## 1090 CHROMs
## 11,689 variants
## Object size: 7.6 Mb
## 0 percent missing data
## *****        *****         *****

# Read in pop file
pop.data <- read.table("popfile_52_samples_4_pops", sep = "\t", header = FALSE)

# We can now check that all the samples in the VCF and the population data frame are included:
all(colnames(homoeophylla_glabra.VCF@gt)[-1] == pop.data$AccessID)

## [1] TRUE

# Convert to a genlight object
gl.homoeophylla <- vcfR2genlight(homoeophylla_glabra.VCF)

## Warning message:
## In vcfR2genlight(homoeophylla_glabra.VCF) :
##   Found 127 loci with more than two alleles.
## Objects of class genlight only support loci with two alleles.
## 127 loci will be omitted from the genlight object.

# Specify ploidy
ploidy(gl.homoeophylla) <- 1

pop(gl.homoeophylla) <- pop.data$V2

gl.homoeophylla

## /// GENLIGHT OBJECT /////////

## // 52 genotypes,  11,562 binary SNPs, size: 1.7 Mb
## 138717 (23.07 %) missing data

## // Basic content
## @gen: list of 52 SNPbin
## @ploidy: ploidy of each individual  (range: 1-1)

## // Optional content
## @ind.names:  52 individual labels
## @loc.names:  11562 locus labels
## @chromosome: factor storing chromosomes of the SNPs
## @position: integer storing positions of the SNPs
## @pop: population of each individual (group size range: 12-14)
## @other: a list containing: elements without names  



### Population genetic analyses for GBS data
# Distance Matrix

homoeophylla.dist <- poppr::bitwise.dist(gl.homoeophylla)

tree <- aboot(gl.homoeophylla, 
              tree = "upgma", 
              distance = bitwise.dist, 
              sample = 100, 
              showtree = F, 
              cutoff = 100, 
              quiet = T)

#tree <- root(tree, outgroup = "15848_PFR")
#tree <- ladderize(tree, right = FALSE)

#cols <- wes_palette("Darjeeling1", n = nPop(gl.homoeophylla))
cols <- c("skyblue4", "plum4", "seagreen4", "goldenrod2")
plot.phylo(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.homoeophylla)], type = "phylogram")
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
legend('topleft', legend = c("AUS", "CHI", "NZ", "PHO"), fill = cols, border = FALSE, bty = "n", cex = 1.3)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")

### Minimum spanning networks

library(igraph)

homoeophylla.dist <- bitwise.dist(gl.homoeophylla)
homoeophylla.msn <- poppr.msn(gl.homoeophylla, homoeophylla.dist, showplot = FALSE, include.ties = T)

node.size <- rep(2, times = nInd(gl.homoeophylla))
names(node.size) <- indNames(gl.homoeophylla)
vertex.attributes(homoeophylla.msn$graph)$size <- node.size

set.seed(9)
plot_poppr_msn(gl.homoeophylla, homoeophylla.msn , palette = c("skyblue4", "plum4", "seagreen4", "goldenrod2"), 
               gadj = 70, inds = "cool", nodelab = 1000)

### Principal components analysis

#homoeophylla.pca <- glPca(gl.homoeophylla, nf = 4)
homoeophylla.pca <- glPca(gl.homoeophylla, nf = 3)
barplot(100*homoeophylla.pca$eig/sum(homoeophylla.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

homoeophylla.pca.scores <- as.data.frame(homoeophylla.pca$scores)
homoeophylla.pca.scores$pop <- pop(gl.homoeophylla)

library(ggplot2)
set.seed(9)
p <- ggplot(homoeophylla.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=3)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()

p

# Obtain variation contribution from eigenvalues for PC1 & PC2
var_frac <- homoeophylla.pca$eig/sum(homoeophylla.pca$eig)
signif(sum(var_frac[1]) * 100, 3) ## 20.2%
signif(sum(var_frac[2]) * 100, 3) ## 14.2%

### DAPC

#homoeophylla.dapc <- dapc(gl.homoeophylla, n.pca = 16, n.da = 2)
homoeophylla.dapc <- dapc(gl.homoeophylla, n.pca = 10, n.da = 2)

scatter(homoeophylla.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomright", scree.pca = TRUE,
        posi.pca = "bottomleft", posi.da = "topleft", cleg = 0.75)

compoplot(homoeophylla.dapc, col = cols, posi = 'top')

dapc.results <- as.data.frame(homoeophylla.dapc$posterior)
dapc.results$pop <- pop(gl.homoeophylla)
dapc.results$indNames <- rownames(dapc.results)

library(reshape2)
dapc.results <- melt(dapc.results)

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p




