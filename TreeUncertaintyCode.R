#### Amy Willis, 2016
#### R code accompanying
#### "Uncertainty in phylogenetic tree estimates", Willis & Bell, 2016
#### Please feel free to contact me with any issues or questions! 

## Section 4.2
n <- 2
m <- 1
Sigma <- matrix(c(2,-1,-1,2), nrow=2)
mu = c(5,5)
nn <- 200
x <- mvrnorm(nn, mu, Sigma)
plot(x, xlim=c(0,10), ylim=c(0,10), pch=16)


## Section 5.2

## Section 6