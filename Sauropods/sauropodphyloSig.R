#Code below is for plotting continuous trait values on phylogeny scaled against geologic time, modified from https://lukejharmon.github.io/ilhabela/instruction/2015/07/05/plotting-methods/

#setwd
require(ape)
require(phytools)
require(scatterplot3d)
require(rgl)
library(geiger)
library(phangorn)
library(RColorBrewer)
require(viridisLite)

#apply branch lengths using code provided by Graeme Lloyd
source("http://www.graemetlloyd.com/pubdata/functions_2.r") 
clade.tree<-read.nexus("tree_SauropodMass12_27_22.nex")
ages<-read.table("sauropod_ages.txt",row.names=1)
branchLengthsTree<-date.phylo(clade.tree, ages, rlen=1, method="equal")

# Save the tree as a landscape PNG
#png("phylogenetic_tree.png", width = 2400, height = 1400, res = 300)  # 8x4.6 inches at 300 dpi

# Plot the time-scaled tree
plotTree(branchLengthsTree, fsize = 0.3, ftype = "i", lwd = 1)

write.tree(branchLengthsTree, file="branchLengthsTree.phy")

sauropod.tree<-read.tree("branchLengthsTree.phy")

sauropod.data<-read.csv("sauropod_masses.csv",row.names=1)
sauropod.data  

bodyMass<-setNames(sauropod.data$body_mass,rownames(sauropod.data))

#phlyogenetic signal with Blomberg's K
phylosig(sauropod.tree,bodyMass)

#randomization test
bmassK<-phylosig(sauropod.tree,bodyMass,test=TRUE,nsim=10000)
bmassK
plot(bmassK)

#phlyogenetic signal with Pagel's lambda
phylosig(sauropod.tree,bodyMass,method="lambda")
  
bmasslambda<-phylosig(sauropod.tree,bodyMass,method="lambda",test=TRUE,nsim=10000)
bmasslambda
plot(bmasslambda)
