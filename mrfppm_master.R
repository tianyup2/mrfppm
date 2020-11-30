rm(list = ls(all = TRUE))
library(Rcpp)
library(spatstat)

# sourceCpp("test_ker.cpp")
# sourceCpp("test_Rsamp.cpp")
library(doParallel)
library(doSNOW)
library(utils)
library(fossil)
library(orthogonalsplinebasis)
library(fdapace)
library(rlang)
library(mrfmfmGP)
memory.limit(size = 3000000)

log_vnt <- function(n,m,lambda,gamma){
  # assume gamma=1, t ranges from 1 to n+1
  vnt_list=rep(0,n+m)
  nTrunc=1e2
  for(t in 1:(n+m)){
    r=-Inf
    for(k in t:(t+nTrunc-1)){
      b=lgamma(k+1)-lgamma(k-t+1)-(lgamma(n+k*gamma)-lgamma(k*gamma))+dpois(k-1,lambda,log = TRUE)
      m=max(b,r)
      r=log(exp(r-m)+exp(b-m))+m
    }
    vnt_list[t]=r
  }
  return(vnt_list)
}

updateClusterAssign_GP <- function(data_y, data_x, data_t, adj_mat,
                                   n, log_Vnt, Up_c, m, result,
                                   a0, b0,
                                   a1, b1,
                                   a2, b2,
                                   Mu0, Lambda0,
                                   lambda_mrf, gamma_mfm){
  clusterAssign_update_GP(data_y, data_x, data_t, adj_mat, 
                          n, result$nCluster, Up_c, m,
                          result$clusterAssign, result$clusterSize, log_Vnt,
                          result$beta, result$sigma2c, result$sigma2l, result$multip,
                          a0, b0, a1, b1, a2, b2,
                          Mu0, Lambda0, lambda_mrf, gamma_mfm)
  return(result)
}

updateClusterParameters_GP <- function(data_y, data_x, data_t, n, n_rep, mh_st, result,
                                       a0, b0, a1, b1, a2, b2,
                                       Mu0, Lambda0){
  clusterParameter_update_GP(data_y, data_x,data_t,
                             n, n_rep, mh_st,
                             result$nCluster, result$clusterAssign,
                             result$beta, result$sigma2c, result$sigma2l, result$multip,
                             a0, b0, a1, b1, a2, b2,
                             Mu0, Lambda0)
  return(result)
}

marginalLikMRFMFM_GP<- function(n_ml,n_burn,
                                data_y,data_x, data_t, adj_mat,
                                clusterAssign, Up_C,
                                a0, b0,
                                a1, b1,
                                a2, b2,
                                Mu0, Lambda0,
                                lambda_mrf=0, lambda_mfm=1,gamma_mfm=1,ifProg=TRUE){
  n=length(data_y)
  
  nCluster <- length(unique(clusterAssign))
  clusterSize <- as.integer(c(table(clusterAssign),rep(0,Up_C-nCluster)))
  log_Vnt <- log_vnt(n,1,lambda_mfm,gamma_mfm)
  # initialize:
  result <- list()
  result$ml_max <- -Inf
  result$log_ml_list <- c(0,-Inf)
  result$clusterAssign <- as.integer(clusterAssign)
  result$clusterSize <- as.integer(clusterSize)
  result$nCluster <- as.integer(nCluster)
  if(ifProg){
    pb <- txtProgressBar(min = 0, max = n_ml, style = 3)
  }
  for(i in 1:n_ml){
    ml_mrf_mfm_GP(n_burn,i,result$ml_max, result$log_ml_list,
                  data_y, data_x, data_t, adj_mat,
                  n, as.integer(Up_C), result$nCluster,
                  result$clusterAssign, result$clusterSize, log_Vnt,
                  a0, b0, a1, b1, a2, b2,
                  Mu0, Lambda0, lambda_mrf, gamma_mfm)
    
    if(ifProg){
      setTxtProgressBar(pb, i)
    }
  }
  if(ifProg){
    close(pb)
  }
  # return(result$log_ml_list[2]-log(n_ml))
  return(result$log_ml_list[2]-log(n_ml))
}

mrf_mfm_GP <- function(data_y,data_x, data_t, n_rep, mh_st, adj_mat, init_C, Up_C, n_iter,
                       a0, b0, a1, b1, a2, b2, Mu0, Lambda0, lambda_mrf=0, gamma_mfm=1, lambda_mfm=1,ifProg=TRUE){
  # data_y: list of response, each element should be in column vector
  # data_x: list of covariate, each element should be in n_i*p matrix
  # adj_mat: n*n matrix, where n is the number of observations
  # init_C: the initial number of clusters, which should be smaller than Up_C
  # Up_C: Upper bound of the number of clusters, should be smaller equal to n
  
  n=length(data_y)
  d=length(Mu0)
  # Initialization
  log_Vnt <- log_vnt(n,1,lambda_mfm,gamma_mfm)
  
  clusterAssign <- as.integer(c(1:init_C,
                                sample(1:init_C, size = n - init_C, replace = TRUE)))
  
  nCluster <- as.integer(init_C)
  clusterSize <- as.integer(c(table(clusterAssign),rep(0,Up_C-nCluster)))
  
  sigma2c_samp <- c(rep(1,nCluster),rep(0,Up_C-nCluster))
  beta_samp <- matrix(0,nrow = d,ncol = Up_C)
  sigma2l_samp <- c(rep(1,nCluster),rep(0,Up_C+1-nCluster))
  multip_samp <- c(rep(0.1,Up_C),0)
  
  result<-list()
  result$beta <- beta_samp
  result$sigma2c <- sigma2c_samp
  result$sigma2l <- sigma2l_samp
  result$multip <- multip_samp
  result$nCluster <- nCluster
  result$clusterAssign <- clusterAssign
  result$clusterSize <- clusterSize
  
  History <- vector("list", n_iter)
  if(ifProg){
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  }
  for(i_iter in 1:n_iter)
  {
    # Update clusterAssign
    result <- updateClusterAssign_GP(data_y, data_x, data_t, adj_mat,
                                     n, log_Vnt, as.integer(Up_C), m, result,
                                     a0, b0, a1, b1, a2, b2,
                                     Mu0, Lambda0,
                                     lambda_mrf, gamma_mfm)
    # update clusterParameters:
    result <- updateClusterParameters_GP(data_y, data_x, data_t, n, n_rep, mh_st, result,
                                         a0, b0, a1, b1, a2, b2,
                                         Mu0, Lambda0)
    
    clusterAssign_temp <- duplicate(result$clusterAssign,shallow = TRUE)
    nCluster_temp <- duplicate(result$nCluster[1],shallow = TRUE)
    beta_temp <- duplicate(result$beta[,1:nCluster_temp],shallow = TRUE)
    sigma2c_temp <- duplicate(result$sigma2c[1:nCluster_temp],shallow = TRUE)
    sigma2l_temp <- duplicate(result$sigma2l[1:nCluster_temp],shallow = TRUE)
    multip_temp <- duplicate(result$multip[1:nCluster_temp],shallow = TRUE)
    clusterSize_temp <- duplicate(clusterSize[1:nCluster_temp],shallow = TRUE)
    
    History[[i_iter]] <- list(zout=clusterAssign_temp,zsize=clusterSize_temp,zn=nCluster_temp,
                              beta_out=beta_temp,sigma2c_out=sigma2c_temp,sigma2l_out=sigma2l_temp,multip_out=multip_temp)
    
    if(ifProg){
      setTxtProgressBar(pb, i_iter)
    }
  }
  if(ifProg){
    close(pb)
  }
  return(History)
}

updateClusterAssign_GP_DP <- function(data_y, data_x, data_t, adj_mat,
                                      n, Up_c, m, result,
                                      a0, b0,
                                      a1, b1,
                                      a2, b2,
                                      Mu0, Lambda0,
                                      lambda_mrf, alpha){
  clusterAssign_update_GP_DP(data_y, data_x, data_t, adj_mat, 
                             n, result$nCluster, Up_c, m,
                             result$clusterAssign, result$clusterSize, 
                             result$beta, result$sigma2c, result$sigma2l, result$multip,
                             a0, b0, a1, b1, a2, b2,
                             Mu0, Lambda0, lambda_mrf, alpha)
  return(result)
}

mrf_mfm_GP_DP <- function(data_y,data_x, data_t, n_rep, mh_st, adj_mat, init_C, Up_C, n_iter,
                          a0, b0, a1, b1, a2, b2, Mu0, Lambda0, lambda_mrf=0, alpha=1,ifProg=TRUE){
  # data_y: list of response, each element should be in column vector
  # data_x: list of covariate, each element should be in n_i*p matrix
  # adj_mat: n*n matrix, where n is the number of observations
  # init_C: the initial number of clusters, which should be smaller than Up_C
  # Up_C: Upper bound of the number of clusters, should be smaller equal to n
  
  n=length(data_y)
  d=length(Mu0)
  # Initialization
  
  clusterAssign <- as.integer(c(1:init_C,
                                sample(1:init_C, size = n - init_C, replace = TRUE)))
  
  nCluster <- as.integer(init_C)
  clusterSize <- as.integer(c(table(clusterAssign),rep(0,Up_C-nCluster)))
  
  sigma2c_samp <- c(rep(1,nCluster),rep(0,Up_C-nCluster))
  beta_samp <- matrix(0,nrow = d,ncol = Up_C)
  sigma2l_samp <- c(rep(1,nCluster),rep(0,Up_C+1-nCluster))
  multip_samp <- c(rep(0.1,Up_C),0)
  
  result<-list()
  result$beta <- beta_samp
  result$sigma2c <- sigma2c_samp
  result$sigma2l <- sigma2l_samp
  result$multip <- multip_samp
  result$nCluster <- nCluster
  result$clusterAssign <- clusterAssign
  result$clusterSize <- clusterSize
  
  History <- vector("list", n_iter)
  if(ifProg){
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  }
  for(i_iter in 1:n_iter)
  {
    # Update clusterAssign
    result <- updateClusterAssign_GP_DP(data_y, data_x, data_t, adj_mat,
                                        n, as.integer(Up_C), m, result,
                                        a0, b0, a1, b1, a2, b2,
                                        Mu0, Lambda0,
                                        lambda_mrf, alpha)
    # update clusterParameters:
    result <- updateClusterParameters_GP(data_y, data_x, data_t, n, n_rep, mh_st, result,
                                         a0, b0, a1, b1, a2, b2,
                                         Mu0, Lambda0)
    
    clusterAssign_temp <- duplicate(result$clusterAssign,shallow = TRUE)
    nCluster_temp <- duplicate(result$nCluster[1],shallow = TRUE)
    beta_temp <- duplicate(result$beta[,1:nCluster_temp],shallow = TRUE)
    sigma2c_temp <- duplicate(result$sigma2c[1:nCluster_temp],shallow = TRUE)
    sigma2l_temp <- duplicate(result$sigma2l[1:nCluster_temp],shallow = TRUE)
    multip_temp <- duplicate(result$multip[1:nCluster_temp],shallow = TRUE)
    clusterSize_temp <- duplicate(clusterSize[1:nCluster_temp],shallow = TRUE)
    
    History[[i_iter]] <- list(zout=clusterAssign_temp,zsize=clusterSize_temp,zn=nCluster_temp,
                              beta_out=beta_temp,sigma2c_out=sigma2c_temp,sigma2l_out=sigma2l_temp,multip_out=multip_temp)
    
    if(ifProg){
      setTxtProgressBar(pb, i_iter)
    }
  }
  if(ifProg){
    close(pb)
  }
  return(History)
}

getDahl <- function(MFMfit, burn)
{
  ################################################################
  
  ## Input: MFMfit = the result from GWCRP ##
  ##        burn = the number of burn-in iterations ##
  
  ## Output:
  ##         zout = estimated clustering configuration, a n by 1 vector##
  ##         Qout = estimated probability matrix, a k by k matrix ##
  
  #################################################################
  iters <- MFMfit[-(1:burn)]
  n <- length(iters[[1]][[1]])
  niters <- length(iters)
  membershipMatrices <- lapply(iters, function(x){
    clusterAssign <- x$zout
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  DahlAns
}

# simulation:
label <- unlist(read.table("label.txt",header = FALSE))
# labelc2 <- unlist(read.table("labelc2.txt",header = FALSE))
adj_mat <- read.csv("us_adj.csv",row.names = 1)
idx <- c(1:7,9:49);
adj_mat <- adj_mat[idx,idx];
adj_mat <- matrix(as.double(as.matrix(adj_mat)),nrow = 48,ncol = 48)

# x matrix, dim=2:
nt <- 20
set.seed(5)
ti <- matrix(seq(-1,1,length.out = nt),ncol = 1)

beta0 <- matrix(c(0,1,
                  28,1,
                  -28,1),ncol = 3)

# beta0 <- matrix(c(0,5,
#                   -20,4,
#                   20,6),ncol = 3)

sigma2c <- c(1,1,1)*36
l2 <- c(1,1,1)*10
multip <- c(1,1,1)*0.1
col_l <- c("red","green","blue")
refct <- diag(c(1,1)*1e-1)

data_y <- list()
data_x <- list()
data_t <- list()
matrix_temp <- matrix(0,nrow = 48,ncol = 2)
for(i in 1:48){
  # s <- seq(-5,5,length.out = nt)
  # s <- rt(nt,df = 1)
  s <- runif(n = nt,min = -5,max = 5)
  
  phi0 <- cbind(matrix(1,nrow = nt,ncol = 1),s)
  if(label[i]==1 & runif(1)<1){
    beta_temp <- t(mvtnorm::rmvnorm(1,beta0[,label[i]],refct))
  }else{
    beta_temp <- beta0[,label[i],drop=FALSE]
  }
  # beta_temp <- beta0[,label[i],drop=FALSE]
  
  mean_temp <- t(mvtnorm::rmvnorm(1,phi0%*%beta_temp,SEker(ti,sigma2c = sigma2c[label[i]],
                                                           sigma2l = l2[label[i]],multip = multip[label[i]])))
  
  data_y[[i]] <- mean_temp
  data_x[[i]] <- phi0
  data_t[[i]] <- ti
  matrix_temp[i,] <- solve(t(phi0)%*%phi0)%*%t(phi0)%*%mean_temp
  if(i==1){
    plot(matrix_temp[i,1],matrix_temp[i,2],col=col_l[label[i]],xlim=c(-40,40),ylim=c(-10,20))
  }else{
    points(matrix_temp[i,1],matrix_temp[i,2],col=col_l[label[i]],xlim=c(-40,40),ylim=c(-10,20))
  }
  # if(i==1){
  #   plot(data_t[[i]],data_y[[i]],"l",col=col_l[label[i]],ylim=c(-15,15))
  # }else{
  #   lines(data_t[[i]],data_y[[i]],col=col_l[label[i]])
  # }
}

n_ml <- 1e6
n_burn <- 1e4
n_iter<-1e3
n_rep <- 30
m <- 10
mh_st <- 1e-2
init_C <- 48
Up_C <- 48
gamma_mfm <- 1
lambda_mfm <- 10
alpha <- 1

a0 <- 1e-1
b0 <- 1
a1 <- 2
b1 <- 1
a2 <- 2
b2 <- 1

Mu0 <- matrix(0,nrow = dim(phi0)[2],ncol = 1)
Lambda0 <- 1e-6*diag(dim(phi0)[2])

lambda_mrf <- 0

set.seed(1)

mrf_mfm_out <- mrf_mfm_GP_DP(data_y,data_x, data_t, n_rep, mh_st, adj_mat, init_C, Up_C, n_iter,
                             a0, b0, a1, b1, a2, b2, Mu0, Lambda0, 
                             lambda_mrf=lambda_mrf, alpha = alpha,ifProg=TRUE)

# mrf_mfm_out <- mrf_mfm_GP(data_y,data_x, data_t, n_rep, mh_st, adj_mat, init_C, Up_C, n_iter,
#                           a0, b0, a1, b1, a2, b2, Mu0, Lambda0, 
#                           lambda_mrf=lambda_mrf, gamma_mfm=gamma_mfm, lambda_mfm=lambda_mfm,ifProg=TRUE)

mrf_mfm_Dahl <- getDahl(mrf_mfm_out,5e2)
rand.index(mrf_mfm_Dahl$zout,label)
library(usmap)
library(ggplot2)
library(ggsci)
stateinfo <- readRDS("state_fips_name.rds")
# clst <- c()
# j <- 1
# count <- 1
# idx_rv <- c(2,9,12)
# for(i in 1:51){
#   if(i==idx_rv[count]){
#     clst[i] <- 0
#     if(count < 3){
#       count <- count + 1
#     }
#   }else{
#     clst[i] <- labelc2[j]
#     j <- j+1
#   }
# }
clst <- c()
j <- 1
count <- 1
idx_rv <- c(2,9,12)
for(i in 1:51){
  if(i==idx_rv[count]){
    clst[i] <- 0
    if(count < 3){
      count <- count + 1
    }
  }else{
    clst[i] <- mrf_mfm_Dahl$zout[j]
    j <- j+1
  }
}
stateinfo$clst <- as.factor(clst)

plot_usmap(data=stateinfo,exclude = stateinfo$abbr[idx_rv],values = "clst", labels = TRUE)+
  scale_fill_lancet(name="clst")+
  labs(title="",y="partition result")

marginalLikMRFMFM_GP(1e6,1e4,data_y,data_x,data_t,adj_mat,
                     mrf_mfm_Dahl$zout, Up_C,
                     a0, b0, a1, b1, a2, b2,
                     Mu0, Lambda0,
                     lambda_mrf=lambda_mrf,lambda_mfm=lambda_mfm, gamma_mfm=gamma_mfm)