#### Amy Willis, November 2016
#### R code accompanying:
#### "Uncertainty in phylogenetic tree estimates", Willis & Bell, 2016

#### Please feel free to contact me with any issues or questions! 

# Folder containing everything contained in github.com/adw96/TreeUncertainty
# Please update this according to your clone!
mywd <- "/Users/adw96/Documents/Phylogenetics/TreeUncertainty/TreeUncertaintyAnalysis/git/TreeUncertainty/"
setwd(mywd)

####################################
######### Section 4.2
####################################

## I have collated the OrthoMam  trees and rates for (your) convenience,
## and calculated the log maps. The original trees and rates are available in
## Data/OrthoMaM_trees.txt and Data/OrthoMaM_rates.txt, the logmaps are in
## Data/orthomam_logmaps.txt, the log map of the base tree (weighted Frechet mean) 
## is in Data/orthomam_wm_lm.txt and 
## the base tree itself is in Data/orthomam_wm_tree.txt. 

## Load in the data
source("Data/orthomam_logmaps.txt")
source("Data/orthomam_wm_lm.txt")
mammals_mean_tree <- tree0
tree_info <- c(read.table("Data/OrthoMaM_rates.txt")$V1) # rates

observations_mammals <- matrix(NA, ncol = length(tree1), nrow = 574)
for (i in 1:dim(observations_mammals)[1])  observations_mammals[i, ] <- get(paste("tree", i, sep=""))
dim(observations_mammals)

## Calculating the BHV  distances took ~16 hours, so I saved the
## distances into Data/bhv_dist_matrices.csv
bhv_dist <- as.matrix(read.csv("Data/bhv_dist_matrices.csv", sep = " ", header = F))
require(MASS)
bhv_fit <- isoMDS(bhv_dist, k=2) # k is the number of dim

## ML estimates
mean_estimate <- apply(observations_mammals/tree_info, 2, sum)/sum(1/tree_info) 
sigma_squared_hat <- sum(apply((observations_mammals - mean_estimate)^2/tree_info, 2, sum))/(dim(observations_mammals)[1])

#### Construct data frame for ggplot
new_x <- svd(t(observations_mammals))$v[, 1:2] ## 2-dim approx of observations_mammals: 1000 x 2
colnames(new_x) = c("x", "y")
mammals_df <- data.frame(new_x, "Gene" = paste("Gene", 1:dim(observations_mammals)[1], sep=""), "Type" = "obs")

## make ellipses, rotate them
angles <- (0:100) * 2 * pi/100; unit.circle2 <- (cbind(cos(angles), sin(angles)))
n <- dim(observations_mammals)[1] ## number of trees
m <- dim(observations_mammals)[2] ## dimension of logmaps

A <- t(observations_mammals) 
R <- svd(A)$u%*%solve(diag(svd(A)$d))[, 1:2]
radius <- sqrt(m * stats::qf(0.95, m, n-1))
for (i in 1:n) {
  sigma_hat1 <- diag(rep(sigma_squared_hat*tree_info[i], m)) ## cov of group
  my_check <- sqrt(sigma_hat1) # since diagonal
  L1_inv <- my_check * radius
  c1 <- observations_mammals[i, ]
  svd1 <- svd(t(R)%*%L1_inv)
  smart_e1 <- t(c(c1%*%R) + svd1$u %*% diag(svd1$d) %*%t(unit.circle2))
  colnames(smart_e1) <- colnames(new_x)
  new_df <- data.frame(smart_e1, "Gene" = paste("Gene", i, sep=""), "Type" = "ellipse")
  mammals_df <- rbind(mammals_df, new_df)
}

require(ggplot2)
gg_base <- ggplot(mammals_df, aes(x, y, Gene, Type)) + guides(fill=FALSE) +
  xlab("First Principal Component") + ylab("Second Principal Component") +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))+ theme(legend.position="none")

r <- gg_base +
  geom_point(data = subset(mammals_df, Type == "obs"), aes(col=Gene), alpha=1) +
  guides(fill=FALSE)



gg_base <- ggplot(mammals_df, aes(x, y, Gene, Type)) + guides(fill=FALSE) +
  xlab("First Principal Component") + ylab("Second Principal Component") +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))+ theme(legend.position="none")

df_mds <- data.frame("x"= bhv_fit$points[,1], "y"= bhv_fit$points[,2], "Gene" = paste("Gene", 1:574, sep="")) 
tt <- ggplot(df_mds, aes(x,y,Gene)) + guides(fill=FALSE) +
  xlab("First MDS Coordinate") + ylab("Second MDS Coordinate") +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))+ theme(legend.position="none") + 
  geom_point(aes(col=Gene), alpha=1)

#### Figure 1!
require(gridExtra)
require(grid)
require(lattice)
grid.arrange(tt + coord_fixed(xlim = c(-10,10), ylim = c(-10,10), ratio = 1), 
             r + coord_fixed(ratio = 1), ncol=2)


q <- gg_base +
  geom_polygon(aes(fill=Gene), alpha=0.1) +
  geom_point(data = subset(mammals_df, Type == "obs"), aes(col=Gene), alpha=1) +
  geom_point(data = subset(mammals_df, Type == "obs"), shape =1 )+  scale_shape(solid = FALSE)+ guides(fill=FALSE)

tree0_b <- mean_estimate
deviations2 <- apply(observations_mammals, 1, function(x) sqrt(sum((x-tree0_b)^2)))
data1 <- data.frame(tree_info, "Gene" =  paste("Gene", 1:574, sep=""), deviations2)
s <-  ggplot(data1, aes(tree_info, deviations2)) + 
  geom_point(aes(col=Gene), alpha=1) + 
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) + 
  xlab("Relative evolutionary rate") + ylab("Residuals (BHV)")+ theme(legend.position="none")


#### Figure 2!
qsize = 1.3
grid.arrange(q + coord_fixed(xlim = c(-qsize,qsize), ylim = c(-qsize,qsize), ratio = 1), 
             s + coord_fixed(xlim = c(0,5.5), ylim = c(0,3.5), ratio = 5.5/3.5), ncol=2)

####################################
######### Section 5.2
####################################

## I similarly log mapped the trees for your convenience.
## If you would like the original scripts or the original trees (there are 10,000)
##  please shoot me an e-mail! :)

source("Data/mapped_terrapene_trees_no_base_all.R")
gene_tree_files <- as.character(read.table("Data/TerrapeneGeneNames.txt")[,1])

for (i in 1:1e7) {
  if (class(try(get(paste("tree", i, sep="")), silent=TRUE)) == "try-error") {
    break;
  }
}
number_trees <- i-1
observations_turtles <- matrix(NA, ncol = length(tree1), nrow = number_trees)
for (i in 1:number_trees) observations_turtles[i, ] <- get(paste("tree", i, sep=""))
apply(observations_turtles, 2, mean)
apply(observations_turtles, 2, sd)

turtles_pcr <- prcomp(observations_turtles, center=T, scale=T) # data are in same units 
projections <- turtles_pcr$x

gene_names <- rep(unlist(lapply(strsplit(gene_tree_files, "_"), function(x) x[2])), each=100)
first_projections <- data.frame("PC1"=projections[,1], "PC2"=projections[,2], "PC3"=projections[,3], "Genes"=gene_names)

#### make gorgeous figure
## make data frame with ellipses
## PCA of all data
new_x <- svd(t(observations_turtles))$v[, 1:2] ## 2-dim approx of observations_turtles_turtles: 1000 x 2
colnames(new_x) = c("x", "y")
turtles_df <- data.frame(new_x, "Gene" = gene_names, "Type" = "obs")

## make ellipses
# find rotation
angles <- (0:100) * 2 * pi/100; unit.circle2 <- (cbind(cos(angles), sin(angles)))

## Projecting the ellipsoids  in accordance with the procedure outlined in
## SuppMat.pdf
A <- t(observations_turtles) ## 7 x 1000 matrix
R <- svd(A)$u%*%solve(diag(svd(A)$d))[, 1:2]
radius <- sqrt(7 * stats::qf(0.95, 7, 99))
for (i in 1:10) {
  y <- observations_turtles[1:100 + 100*(i-1), ]
  sigma_hat1 <- cov(y) ## cov of group
  my_check <- try(chol(sigma_hat1), silent = T)
  if (class(my_check) == "try-error") {
    Q <- chol(sigma_hat1, pivot = T)
    pivot <- attr(Q, "pivot")
    my_check <- Q[, order(pivot)] # recover m
  }
  L1_inv <- my_check * radius
  c1 <- colMeans(y)
  svd1 <- svd(t(R)%*%L1_inv)
  smart_e1 <- t(c(c1%*%R) + svd1$u %*% diag(svd1$d) %*%t(unit.circle2))
  colnames(smart_e1) <- colnames(new_x)
  new_df <- data.frame(smart_e1, "Gene" = unique(gene_names)[i], "Type" = "ellipse")
  turtles_df <- rbind(turtles_df, new_df)
}

levels(turtles_df$Gene)[levels(turtles_df$Gene)=="cytb"] <- "Cyt-b"
levels(turtles_df$Gene)[levels(turtles_df$Gene)=="Gapd"] <- "GAPD"
levels(turtles_df$Gene)[levels(turtles_df$Gene)=="vim"] <- "Vim"
levels(turtles_df$Gene)
relevel(turtles_df$Gene, ref="Cyt-b")

gene1 <- turtles_df$Gene
unique(gene1)
tmp <- factor(turtles_df$Gene, levels = unique(gene1)[c(3, 1:2, 4:10)])
turtles_df$Gene <- tmp


levels(turtles_df$Gene)[levels(turtles_df$Gene)!="Cyt-b"] <- paste("N:", levels(turtles_df$Gene)[levels(turtles_df$Gene)!="Cyt-b"] )
levels(turtles_df$Gene)[levels(turtles_df$Gene)=="Cyt-b"] <- "M: Cyt-b"

## Figure 3
xlim1 <- 0.25
ggplot(turtles_df, aes(x, y, Gene, Type)) +
  geom_point(data = subset(turtles_df, Type == "obs"), aes(x,y, col=Gene)) + 
  geom_polygon(data = subset(turtles_df, Type == "ellipse"), aes(col=Gene, fill = Gene), alpha=0.1) +
  xlab("First Principal Component") + ylab("Second Principal Component")


####################################
######### Section 6: "Distortions..."
####################################
require(phangorn)
require(ape)

number_of_leaves <- 50
set.seed(1)
base_tree <- rtree(number_of_leaves, rooted=F, br = 1)
trees<- nni(base_tree)


### Credit: Katie Everson, http://www.kmeverson.org/blog/visualizing-tree-space-in-r
tree.dist.matrix <- function(trees, treenames=names(trees)){
  N <- length(trees)
  
  if(N != length(treenames)){
    stop("Names and tree list must be the same length")
  }
  
  #Create an empty matrix for results
  RF <- matrix(0, N, N)
  
  for(i in 1:(N-1)){
    #print(paste("Tree", i, "of", N))        
    for(j in (i+1):N){
      RFd <- RF.dist(trees[[i]],trees[[j]])
      if(RFd==0) RFd = 0.000000001
      RF[i,j]<-RF[j,i]<-RFd
    }
  }
  
  #Row and column names
  rownames(RF) <- treenames
  colnames(RF) <- treenames
  
  RF
}

subset_distances <- tree.dist.matrix(trees, treenames = paste("tree", 1:length(trees), sep=""))
fit <- isoMDS(subset_distances, k=3) # k is the number of dim
# write.tree(base_tree, "equal_distance_base.txt")
# write.tree(trees, "equal_distance.txt")

## code to calculate logmaps is available in my other github directory, or you can contact me :)

# system("java -jar /Users/adw96/Documents/Phylogenetics/TreeUncertainty/TreeUncertaintyAnalysis/logmap_base.jar equal_distance_base.txt equal_distance_base.txt > equal_distance_maps_base.R")
source("Data/equal_distance_maps_base.R"); tree_base_coords <- tree1
# system("java -jar /Users/adw96/Documents/Phylogenetics/TreeUncertainty/TreeUncertaintyAnalysis/logmap_base.jar equal_distance_base.txt equal_distance.txt > equal_distance_maps2.R")
source("Data/equal_distance_maps2.R")
logmap_observations <- matrix(NA, ncol = length(tree1), nrow = length(trees))
for (i in 1:length(trees))  logmap_observations[i, ] <- get(paste("tree", i, sep=""))
equally_distant_pcr <- prcomp(logmap_observations, center = T, scale = TRUE) 

## Figure 4
require(scatterplot3d)
par(mfrow=c(1,2))
scatterplot3d(rbind(c(0,0,0), fit$points[,1:3]), pch = 19, color = c("red", rep("green4", dim(logmap_observations[,1:3])[1])), main="", xlab="First MDS coordinate", ylab="Second MDS coordinate", zlab="Third MDS coordinate", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), zlim=c(-1.5,1.5), cex.lab=0.8, cex.axis=0.6)
scatterplot3d(rbind(tree_base_coords[1:3], logmap_observations[,1:3]), pch = 19, color = c("red", rep("green4", dim(logmap_observations[,1:3])[1])), main="", xlab="First log map coordinate", ylab="Second log map coordinate", zlab="Third log map coordinate", xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1), cex.lab=0.8, cex.axis=0.6)


####################################
######### Section 6: Topology
####################################

## look at orthomam_wm_tree to tell what each branch is doing
## (need to map up with orthomam_wm_lm)
hist(observations_mammals[,40]) ## platypus from other marsupials
hist(observations_mammals[,23]) ## human from gorilla and chimp
tree0 <- mammals_mean_tree

files_mammals <- as.character(read.table("Data/orthomam_gene_names.txt")[[1]])
mammals_gene_names <- unlist(lapply(strsplit(files_mammals, "_"), FUN = function(x) x[2]))
mammals_gene_names <- unlist(lapply(strsplit(mammals_gene_names, ".xml"), FUN = function(x) x[1]))

contentious_data_frame <- data.frame("Branch1" = observations_mammals[,40] + tree0[40], 
                                     "Branch2" = observations_mammals[,23] + tree0[23],
                                     "Genes"= mammals_gene_names)

## Figure 5
ggplot(contentious_data_frame, aes(x = Branch1, y = Branch2, label = Genes)) + 
  geom_hline(aes(yintercept=0)) +  geom_vline(aes(xintercept=0)) +
  geom_point(aes(), size = 0.7) + 
  scale_x_continuous("Platypus clade", limits = c(-2,2)) +
  scale_y_continuous("Human clade", limits = c(-0.3,0.3)) +
  geom_text(data = subset(contentious_data_frame, Branch2 < -0.2+ tree0[23]),  nudge_x = 0.6) +
  geom_text(data = subset(contentious_data_frame, Genes == "DHCR24"),  nudge_x = 0.2, nudge_y = -0.04) +
  geom_text(data = subset(contentious_data_frame, Branch1 < -1.7+ tree0[40]),  nudge_x = 0.4, nudge_y = 0.05) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))


####################################
######### Summary statistics
####################################

## proportion of trees supporting the branches in figure 5
mean(contentious_data_frame$Branch1 > 0)
mean(contentious_data_frame$Branch2 > 0)

## proportion of trees in cone path
percentage_cone_path <- function(observations_matrix) {
  number_negative_coordinates <- apply(X = observations_matrix, 1, function(x) sum(x < 0))
  mean(number_negative_coordinates / dim(observations_matrix)[2] == 1)
}
percentage_cone_path(observations_mammals) # 0, unsurprising, since there are so many branches
percentage_cone_path(observations_turtles) # 23%

## proportion of trees with same topology
percentage_concordant <- function(observations_matrix) {
  number_nonnegative_coordinates <- apply(X = observations_matrix, 1, function(x) sum(x >= 0))
  mean(number_nonnegative_coordinates / dim(observations_matrix)[2] == 1)
}
percentage_concordant(observations_mammals) # 0 
percentage_concordant(observations_turtles) # 6% 
# Fascinating! Let's look at the difference in the distribution
hist(apply(X = observations_mammals, 1, function(x) sum(x >= 0))/ dim(observations_mammals)[2])
hist(apply(X = observations_turtles, 1, function(x) sum(x >= 0))/ dim(observations_turtles)[2])
# This is highly multivariate information as well! 
# I wonder how we could use these histograms to get better pictures of topological differences...
