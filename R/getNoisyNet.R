# This function creates the observed (noisy) network.

# data.true: adjacency matrix of the true network
# Nv: number of nodes
# alpha: probability of having an egde in the observed network when non-edge  in the true network
# beta: probability of not having an egde in the observed network when having an edge in the true network

getNoisyNet <- function(data.true,Nv,alpha,beta)
{
    data.obs <- rsparsematrix(Nv, Nv, nnz = rbinom(1,Nv^2,alpha), rand.x = function(n) return(1))
    ind1 <- which(data.true==1)
    data.obs[ind1] <- rbinom(length(ind1),1,1-beta)
    diag(data.obs)<- 0
    data.obs <-  Matrix::forceSymmetric(data.obs,uplo="L")
    return(data.obs)
}
