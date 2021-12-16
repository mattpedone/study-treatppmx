library(parallel)
library(doParallel)
library(dplyr)
library(tibble)
library(tidyselect)
library(ggplot2)
library(vweights)
library(treatppmx)
set.seed(121)

K <- 10
nobs <- 50

gendata <- function(n = 1000, pro = c(0.2,0.5,0.3), dim = 3){
  if(length(pro) != 3){
    stop("The length of probabilities must be equal to 3!")
  }
  data <- matrix(nrow = n, ncol = dim)
  for(i in 1:n){
    u <- runif(1)
    if(u < pro[1]){
      #cat(i,"ciao1","\n")
      data[i,] <- mvtnorm::rmvnorm(1, mean = rep(-2.1, dim), sigma = diag(.5, dim))
    }
    else{
      if(u < (pro[1] + pro[2])){
        #cat(i,"ciao2","\n")
        data[i,] <- mvtnorm::rmvnorm(1, mean = rep(0, dim), sigma = diag(.5, dim))
      }
      else{
        #cat(i,"ciao3","\n")
        data[i,] <- mvtnorm::rmvnorm(1, mean = rep(2.3, dim), sigma = diag(.5, dim))
      }
    }
  }
  return (data)
}

X <- gendata(n = nobs, dim = 25)

#X <- data.frame(simupats)
#X <- X[,c(1:25)]
nobs <- nrow(X)

#X <- data.frame(scale(matchRTComp[,16:38]))
#nobs <- nrow(X)

vec_par <- c(0.0, 2.0, .5, 1.0, 2.0, 2.0, 0.1)
#double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
iterations <- 1000; burnin <- 0; thinning <- 1

#beta <- 4.4185; sigma <- .25; theta <- 19.233
beta <- 48.4185; sigma <- .25; theta <- 19.233
#beta <- 1.0; sigma <- .7553; theta <- 19.233

nout <- (iterations-burnin)/thinning

cor_all <- parallel::detectCores()-1#cores to be allocated
registerDoParallel(cores = cor_all)

res_121 <- foreach(sub = 1:K, .combine = rbind) %dopar%
  {
    out_ppmx_prior <- prior_ppmx(X = X, PPMx = 1, cohesion = 2, alpha = beta, #*sigma,
                                 sigma = sigma, similarity = 2, consim = 2,
                                 similparam = vec_par, calibration = 0,
                                 coardegree = 1, iter = iterations, burn = burnin,
                                 thin = thinning, nclu_init = 30)
  }

res_221 <- foreach(sub = 1:K, .combine = rbind) %dopar%
  {
    out_ppmx_prior <- prior_ppmx(X = X, PPMx = 1, cohesion = 1, alpha = beta, #*sigma,
                                 sigma = sigma, similarity = 2, consim = 2,
                                 similparam = vec_par, calibration = 2,
                                 coardegree = 2, iter = iterations, burn = burnin,
                                 thin = thinning, nclu_init = 30)
  }

res_222 <- foreach(sub = 1:K, .combine = rbind) %dopar%
  {
    out_ppmx_prior <- prior_ppmx(X = X, PPMx = 1, cohesion = 2, alpha = beta, #*sigma,
                                 sigma = sigma, similarity = 2, consim = 2,
                                 similparam = vec_par, calibration = 2,
                                 coardegree = 2, iter = iterations, burn = burnin,
                                 thin = thinning, nclu_init = 30)
  }

res_0 <- foreach(sub = 1:K, .combine = rbind) %dopar%
  {
    out_ppmx_prior <- prior_ppmx(X = X, PPMx = 0, cohesion = 1, alpha = theta,
                                 sigma = sigma, similarity = 1, consim = 1,
                                 similparam = vec_par, calibration = 0,
                                 coardegree = 1, iter = iterations, burn = burnin,
                                 thin = thinning, nclu_init = 30)
  }

res <- t(rbind(apply(res_121/nout, 2, mean), #apply(res_122/nout, 2, mean), 
               apply(res_221/nout, 2, mean), apply(res_222/nout, 2, mean), 
               apply(res_0/nout, 2, mean), apply(res_121/nout, 2, sd), 
               #apply(res_122/nout, 2, sd), 
               apply(res_221/nout, 2, sd), 
               apply(res_222/nout, 2, sd), apply(res_0/nout, 2, sd), 
               c(round(vweights::computepnclu(nobs, sigma, beta*sigma), 5))))

res <- (res[which(rowSums(res) != 0),])
res <- rownames_to_column(as.data.frame(res), var = "cluster")
newres <- cbind(similarity = c(rep("NGG-dd-nocal", nrow(res)), 
                               #rep("NGG-dd-cal", nrow(res)), 
                               rep("DP-dd-coa2", nrow(res)), 
                               rep("NGG-dd-coa2", nrow(res)), 
                               rep("DP", nrow(res)), rep("NGG", nrow(res))),
                cluster = rep(res[,1], 5),
                freq = c(unlist(res[, c(2:5, 10)])),
                sd = c(unlist(res[, c(6:9)]), rep(0, nrow(res))))

dfres <- newres %>%
  as_tibble() %>% # lo trasformo in un dataframe
  #tibble::rownames_to_column(var = "cluster") %>%
  mutate(cluster = as.numeric(cluster)) %>%
  mutate(freq = as.numeric(freq)) %>%
  mutate(sd = as.numeric(sd))

ggplot(dfres, aes(x=cluster, y=freq, group=similarity, color=similarity)) +
  geom_errorbar(aes(ymin=freq-sd, ymax=freq+sd), width=.1) +
  geom_line() +
  geom_point()+
  scale_color_brewer(palette="Paired")+
  theme_minimal() +
  #xlab(expression(c)) +
  #ylab(expression(P(C[n] == c)))
  labs(x = expression(c), y = expression(P(C[n] == c)), 
       color = expression(' '))

#ggsave("output/prior-ppmx/plot_.pdf")

#save(res, file = "output/prior-ppmx/results_.RData")

