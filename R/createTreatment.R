library(Matrix)

# set a destination for output
outdir <- "../../Data/"

# probability of treatment
p <- 0.1
# number of Monte Carlo simulation trials
m <- 10000

set.seed(1)

# School 1
z <- matrix(rbinom(153*m, 1, p),nrow=m)
z1 <- Matrix(z, sparse = TRUE) 

# School 2
z <- matrix(rbinom(147*m, 1, p),nrow=m)
z2 <- Matrix(z, sparse = TRUE) 

# School 3
z <- matrix(rbinom(159*m, 1, p),nrow=m)
z3 <- Matrix(z, sparse = TRUE) 

# School 4
z <- matrix(rbinom(76*m, 1, p),nrow=m)
z4 <- Matrix(z, sparse = TRUE) 
 
save(z1,z2,z3,z4, file = "z.RData")

