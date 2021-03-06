# PGL-PHO_sp_pair
## Codes used to study species pairs of *Psuedocyphellaria glabra* and *P. homoeophylla* with RADseq data

### [RAxML code](./PGL_PHO_PFR_RAxML.sh)
Command for RAxML analysis:
Input data: [https://www.dropbox.com/s/fphyg64uqk5xbbv/PGL_PHO_PFR.phy?dl=0]

### [Population genomic analyses](./PGL_PHO_4pops_gen_1Dec2021.R)
R codes for population genomic analyses: Data conversion, distance matrix, minimum spanning network (MSN), principal components analysis (PCA), and discriminant analysis of principal components analysis (DAPC):
Input data: [VCF file](./PGL_PHO_52.recode.vcf) &
[Population file](./popfile_52_samples_4_pops)

### [STRUCTURE analysis](./STRUCTURE_PGL_PHO.py)
Python code for STRUCTURE analysis:
Input data [.hdf5 file](./PGL_PHO.snps.hdf5)

