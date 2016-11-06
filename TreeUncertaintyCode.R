#### Amy Willis, 2016
#### R code accompanying
#### "Uncertainty in phylogenetic tree estimates", Willis & Bell, 2016
#### Please feel free to contact me with any issues or questions! 

# folder containing everything on github
mywd <- "/Users/adw96/Documents/Phylogenetics/TreeUncertainty/TreeUncertaintyAnalysis/git/TreeUncertainty/"
setwd(mywd)

## Section 4.2
## I have collated the OrthoMam  trees and rates for (your) convenience,
## and calculated the log maps. The original trees and rates are available in
## Data/OrthoMaM_trees.txt and Data/OrthoMaM_rates.txt, the logmaps are in
## Data/orthomam_logmaps.txt, the log map of the base tree (weighted Frechet mean) 
## is in Data/orthomam_wm_lm.txt and 
## the base tree itself is in Data/orthomam_wm_tree.txt. 

## Load in the data
source("Data/orthomam_logmaps.txt")
source("Data/orthomam_wm_lm.txt")
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

## Section 5.2

## Section 6