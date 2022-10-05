#' @importFrom stats optim
NULL
#-----------------------------------------------------#
# Functions for reference estimation                  #
#-----------------------------------------------------#

## Check zero counts
if.zero <- function(x) return(sum(x == 0) == 0)

## Estimating cell-type specific mean expression
meanEst <- function(ct, logX, cts){
  Indx <- which(cts == ct)
  res  <- rowMeans(logX[, Indx])
  return(res)
}
## Estimating cell-type specific var expression
varEst <- function(ct, logX, cts){
  Indx <- which(cts == ct)
  res  <- apply(logX[, Indx], MARGIN = 1, var)
  return(res)
}
## Merge results from lapply
mergeList <- function(x) t(do.call(rbind, x))


#' @title Estimates cell-type-specific mean expression and variability.
#'
#' @description This function is designed to estimate the cell-type
#' specific proportions and measure intra-tumor heterogeneity from
#' multiple mixed samples from patients.
#'
#' @param X A cell-type-specific reference matrix. It is a \eqn{G}
#' by \eqn{M} matrix, where \eqn{G} is the number of genes and \eqn{M}
#' is the number of samples for reference.
#' @param cts A vector to label cell type in the reference matrix.
#' It is a \eqn{M} by 1 vector.
#'
#' @return mu_gk A matrix of estimated cell-type-specific mean expression.
#' It is a \eqn{G} by \eqn{K} matrix, where \eqn{K} is the number
#' of cell types of interests.
#' @return var_gk A matrix pf estimated cell-type-specific variability.
#' It is a \eqn{G} by \eqn{K} matrix.
#'
#' @author Peng Yang, Ziyi Li.
#'
#' @export
refEst <- function(X, cts){

  ## Check input values
  if(ncol(X) != length(cts)){
    stop("column of X must be the length of cts.")
  }

  if(any(is.na(X))){
    stop("matrix X must not contain any NA entries.")
  }

  if(any(X < 0)){
    stop("matrix X should be non-negative.")
  }

  ## Log-transformation of X
  if(if.zero(X)){
    logX = log(X)
  }else{
    message("Adding 1 to X to ensure valid log transformation.")
    logX = log(X + 1)
  }

  K = length(unique(cts)) ## Number of cell types

  res_mean <- mergeList(lapply(X = 1:K, FUN = meanEst,
                               logX = logX, cts = cts))

  res_var  <- mergeList(lapply(X = 1:K, FUN = varEst,
                               logX = logX, cts = cts))

  colnames(res_mean) <- colnames(res_var) <- paste('CT', 1:K)

  return(list(mu_gk  = res_mean,
              var_gk = res_var))
}



#-----------------------------------------------------#
# Functions for general use                           #
#-----------------------------------------------------#
## Find the tilde xi
tilde_xi <- function(xi_sk) xi_sk/rowSums(xi_sk)

## Find the xi_s0
xi_sum <- function(xi_sk) rowSums(xi_sk)

## Find sample IDs in patient i
SampleID <- function(I_i, i) which(I_i == i)

## Make a vector 1 by K to G by K
vec2mat <- function(vec, p){
  res <- do.call(rbind, lapply(X = 1:p,
                               FUN = function(x, vec){vec},
                               vec = vec))
  return(res)
}

## Make a S by G matrix to a tensor of S by G by P
mat2ten <- function(mat, p){
  res <- array(0, dim = c(dim(mat), p))
  res[, , 1:p] <- mat
  return(res)
}

## Compute number of samples for each patient
findnSample <- function(I_i){
  res <- unlist(lapply(X = 1:N, FUN = function(x) sum(I_i == x)))
  return(res)
}

## Function for column variance
colVar <- function(x) apply(x, 2, var)

## Transfer an int to a tensor with dimenstion of p
int2ten <- function(int, p){
  res    <- array(int, p)
  return(res)
}

## Sample-wise variability
sampleVar <- function(logY){
  res        <- apply(logY, 1, var)
  names(res) <- paste('Gene', 1:dim(logY)[1])
  return(res)
}

#-----------------------------------------------------#
# Define loss functions:                              #
#   Q1:  ExpQ_logY                                    #
#   Q2:  ExpQ_logW                                    #
#   Q3:  ExpQ_logH                                    #
#   Q4:  ExpQ_logQW                                   #
#   Q5:  ExpQ_logQH                                   #
#-----------------------------------------------------#
#-----------------------------------------------------#
# Functions for loss ExpQ_logY                        #
#-----------------------------------------------------#

## Expectation of w_sgj and w_sgl, where j != l
## if condition for grad of gamma and tau
ExpQ_w_jl <- function(gamma_gk, tau_gk, j, l, if.tau){
  res <- exp(gamma_gk[,j] + gamma_gk[,l] +
               (tau_gk[,j]^2 + tau_gk[,l]^2)/2)
  if(if.tau) res <- res*tau_gk[,l]
  return(res)
}

## Cov of h_js and h_ls, where j != l
Cov_h_jl <- function(xi_k_tilde, xi_0, j, l,
                     if.xi_lk, if.xi_k){

  res <- -xi_k_tilde[l]*xi_k_tilde[j]/(xi_0 + 1)
  if(if.xi_lk){
    xi_k <- xi_k_tilde * xi_0
    res  <- xi_k[l]*xi_k[j]*(3*xi_0 + 2) /
      (xi_0^3*(xi_0 + 1)^2)

  }
  if(if.xi_k){
    xi_k <- xi_k_tilde * xi_0
    res  <- xi_k[l]*xi_0*(xi_0 + 1) /
      (xi_0^3 * (xi_0 + 1)^2)

  }
  return(res)
}

## Product of ExpQ_w and Cov_h
Prod_ExpQ_Cov <- function(l, j, gamma_gk, tau_gk, xi_k_tilde, xi_0,
                          if.tau, if.xi_lk, if.xi_k){
  res <- ExpQ_w_jl(gamma_gk, tau_gk, j, l, if.tau) *
    Cov_h_jl(xi_k_tilde, xi_0, j, l, if.xi_lk, if.xi_k)

  return(res)
}
Cross_Prod_ExpQ_Cov <- function(j, gamma_gk, tau_gk, xi_k_tilde, xi_0,
                                if.tau, if.xi_lk, if.xi_k){
  Indx_l <- c(1:K)[-j]
  res <- do.call(rbind, lapply(X = Indx_l, FUN = Prod_ExpQ_Cov,
                               j = j, gamma_gk = gamma_gk, tau_gk = tau_gk,
                               xi_k_tilde = xi_k_tilde, xi_0 = xi_0,
                               if.tau = if.tau, if.xi_lk = if.xi_lk,
                               if.xi_k = if.xi_k))
  return(colSums(res))
}

## sum across k = j where j != l
ExpQ_Cov <- function(gamma_gk, tau_gk, xi_k_tilde, xi_0,
                     if.tau = FALSE,  if.xi_lk = FALSE, if.xi_k = FALSE){

  res <- do.call(rbind, lapply(X = 1:K, FUN = Cross_Prod_ExpQ_Cov,
                               gamma_gk = gamma_gk, tau_gk = tau_gk,
                               xi_k_tilde = xi_k_tilde, xi_0 = xi_0,
                               if.tau = if.tau, if.xi_lk = if.xi_lk,
                               if.xi_k = if.xi_k))
  return(colSums(res))
}

## Var sum across k
VarQ_wh_k <- function(gamma_gk, tau_gk, xi_k_tilde, xi_0){
  tmp1 <- exp(2*gamma_gk + 2*tau_gk^2)
  tmp2 <- vec2mat((xi_k_tilde*(1-xi_k_tilde)/(xi_0 + 1) + xi_k_tilde^2), nG)
  tmp3 <- exp(2*gamma_gk + tau_gk^2)*vec2mat(xi_k_tilde^2, nG)
  return(rowSums(tmp1*tmp2 - tmp3))
}

## Expectation of sum_k w_sgk * h_ks
ExpQ_wh <- function(gamma_igk, tau_gk, xi_sk, I_i){
  ## Store data
  res <- array(0, dim = c(nS, nG))
  ## Compute normalized xi
  xi_sk_tilde <- tilde_xi(xi_sk = xi_sk)

  for(s in 1:nS){
    ## Patient ID
    i <- I_i[s]
    ## Sum across cell types
    res[s, ] <- rowSums(exp(gamma_igk[i, , ] + tau_gk^2) *
                          vec2mat(xi_sk_tilde[s, ], nG))

  }

  return(res)
}

## Variance of sum_k w_sgk * h_ks
VarQ_wh <- function(gamma_igk, tau_gk, xi_sk, I_i){
  ## Store data
  res <- array(0, dim = c(nS, nG))
  ## Compute normalized xi
  xi_sk_tilde <- tilde_xi(xi_sk = xi_sk)
  ## Compute sum xi_i
  xi_s0 <- xi_sum(xi_sk = xi_sk)

  for(s in 1:nS){
    ## Patient ID
    i <- I_i[s]
    ## Sum across cell types k
    res[s, ] <- VarQ_wh_k(gamma_gk = gamma_igk[i, , ], tau_gk = tau_gk,
                          xi_k_tilde = xi_sk_tilde[s, ], xi_0 = xi_s0[s])  +
      ## Sum across cell type l != k
      ExpQ_Cov(gamma_gk = gamma_igk[i, , ], tau_gk = tau_gk,
               xi_k_tilde = xi_sk_tilde[s, ], xi_0 = xi_s0[s],
               if.tau = FALSE)
  }

  return(res)

}

# Functions for loss ExpQ_logY
ExpQ_logY <- function(logY, gamma_igk, tau_gk, xi_sk, lambda_sg, I_i){
  res <- 0
  ## Compute partA
  ExpQ_wh_res  <- ExpQ_wh(gamma_igk = gamma_igk, tau_gk = tau_gk,
                          xi_sk = xi_sk, I_i = I_i)
  VarQ_wh_res  <- VarQ_wh(gamma_igk = gamma_igk, tau_gk = tau_gk,
                          xi_sk = xi_sk, I_i = I_i)
  partA        <- VarQ_wh_res / ExpQ_wh_res^2
  ## Compute partB
  partB        <- (t(logY) - log(ExpQ_wh_res) + 1/2 * partA)^2
  # print(dim(lambda_sg));
  # print(dim(partA))
  res <- - sum(t(lambda_sg)/2*(partA + partB))  #+ sum(log(lambda_sg/(2*pi)))/2
  return(res)
}

#-----------------------------------------------------#
# Functions for loss ExpQ_logW                        #
#-----------------------------------------------------#

ExpQ_logW <- function(mu_gk, gamma_igk, tau_gk,
                      alpha_gk, beta_gk, rho_gk, I_i){
  ## Compute alpha_n and number of samples per patient
  # alpha_n <-  alpha_gk + length(unique(I_i))/2
  alpha_n <-  alpha_gk + length(c(I_i))/2
  nI_i    <- findnSample(I_i)
  res <- contSum <- contSumI <- 0
  ## Compute the constant sum across patients
  for(i in 1:N){
    contSum <- contSum + (nI_i[i] - 1)*tau_gk^2 +
      rho_gk*nI_i[i]/(nI_i[i] + rho_gk) *
      (1/nI_i[i] * tau_gk^2 + (gamma_igk[i, , ] - mu_gk)^2)
    # contSumI <- contSumI + log(nI_i[i] + rho_gk)
  }
  ## Compute the loss across patient i
  res <- - sum(alpha_n * log(beta_gk + 1/2*contSum))
  ## Compute the constants
  # res <- res + sum(log(gamma(alpha_n)) -
  #                    log(gamma(alpha_gk)) + alpha_gk*log(beta_gk) +
  #                    N/2*log(rho_gk) - N/2*log(2*pi) - contSumI/2)

  return(res)

}

#-----------------------------------------------------#
# Functions for loss ExpQ_logH                        #
#-----------------------------------------------------#
## Define multivariate beta function
# mvgamma <- function(alpha){
#   res <- unlist(lapply(alpha, FUN = function(x) gamma(x)))
#   return(prod(res))
# }
# mvbeta <- function(alpha){
#   res <- mvgamma(alpha)/gamma(sum(alpha))
#   return(res)
# }
lmvbeta <- function(alpha){
  res <- sum(lgamma(alpha))- lgamma(sum(alpha))
  return(res)
}

## Functions for loss ExpQ_logH
ExpQ_logH <- function(xi_sk, C_i, pi_k, I_i){
  res <- 0
  for(s in 1:nS){
    i <- I_i[s]
    res <- res - lmvbeta(C_i[i]*pi_k) +
      sum((C_i[i]*pi_k - 1)*(digamma(xi_sk[s, ]) - digamma(sum(xi_sk[s, ]))))
  }
  return(res)
}

#-----------------------------------------------------#
# Functions for loss ExpQ_QlogW                       #
#-----------------------------------------------------#

ExpQ_QlogW <- function(tau_gk){
  res <- -nS*sum(log(tau_gk^2))
  return(res)
}

#-----------------------------------------------------#
# Functions for loss ExpQ_QlogH                       #
#-----------------------------------------------------#

ExpQ_QlogH <- function(xi_sk){
  res <- 0
  for(s in 1:nS){
    res <- res - lmvbeta(xi_sk[s, ]) +
      sum((xi_sk[s, ] - 1) *(digamma(xi_sk[s, ]) - digamma(sum(xi_sk[s, ]))))
  }
  return(res)
}

#-----------------------------------------------------#
# Define the gradient of loss functions:              #
#   G1:  Grad_gamma: nS by nG by K                    #
#   G2:  Grad_tau                                     #
#   G3:  grad_xi                                      #
#-----------------------------------------------------#
#-----------------------------------------------------#
# Functions for graident of  gamma                    #
#-----------------------------------------------------#
## Part 1
## grident of gamma from ExpQ_logY
#-----------------------------------------------------#
## gradient of ExpQ_wh
grad_ExpQ_wh_gamma <- function(gamma_igk, tau_gk, xi_sk, I_i){
  ## Store data
  res <- array(0, dim = c(nS, nG, K))
  ## Compute normalized xi
  xi_sk_tilde <- tilde_xi(xi_sk = xi_sk)
  ## Compute sum xi_i
  xi_s0 <- xi_sum(xi_sk = xi_sk)

  ## Across samples
  for(s in 1:nS){
    i <- I_i[s]
    res[s, , ] <- exp(gamma_igk[i, ,] + tau_gk^2/2)*vec2mat(xi_sk_tilde[s, ], nG)
  }
  # res <-
  return(res)
}

## gradient of Var_wh for cell type k
grad_VarQ_wh_k_gamma <- function(gamma_gk, tau_gk, xi_k_tilde, xi_0){
  tmp1 <- 2*exp(2*gamma_gk + 2*tau_gk^2)
  tmp2 <- vec2mat(xi_k_tilde*(1-xi_k_tilde)/(xi_0 + 1) + xi_k_tilde^2, nG)
  tmp3 <- 2*exp(2*gamma_gk + tau_gk^2)*vec2mat(xi_k_tilde^2, nG)
  return(tmp1*tmp2 - tmp3)
}

## gradient of Cov_wh k != l
grad_CovQ_wh_gamma <- function(gamma_gk, tau_gk, xi_k_tilde, xi_0){

  res <- do.call(rbind, lapply(X = 1:K, FUN = Cross_Prod_ExpQ_Cov,
                               gamma_gk = gamma_gk, tau_gk = tau_gk,
                               xi_k_tilde = xi_k_tilde, xi_0 = xi_0,
                               if.tau = FALSE, if.xi_lk = FALSE,
                               if.xi_k = FALSE))
  return(t(res))
}

## gradient of Var_wh: nS, nG, K
grad_VarQ_wh_gamma <- function(gamma_igk, tau_gk, xi_sk, I_i){
  ## Store data
  res <- array(0, dim = c(nS, nG, K))
  ## Compute normalized xi
  xi_sk_tilde <- tilde_xi(xi_sk = xi_sk)
  ## Compute sum xi_i
  xi_s0 <- xi_sum(xi_sk = xi_sk)

  ## Across samples
  for(s in 1:nS){
    ## Patient ID
    i <- I_i[s]
    ## Sum across cell types k
    res[s, , ] <- grad_VarQ_wh_k_gamma(gamma_gk = gamma_igk[i, , ], tau_gk = tau_gk,
                                       xi_k_tilde = xi_sk_tilde[s, ], xi_0 = xi_s0[s])  +
      ## Sum across cell type l != k
      grad_CovQ_wh_gamma(gamma_gk = gamma_igk[i, , ], tau_gk = tau_gk,
                         xi_k_tilde = xi_sk_tilde[s, ], xi_0 = xi_s0[s])
  }
  return(res)
}

## grident of gamma from ExpQ_logY
grad_logY_gamma <- function(logY, lambda_sg, gamma_igk, tau_gk, xi_sk, I_i){
  ## Store data
  res <- array(0, dim = c(N, nG, K))

  ## partial results for gradient
  grad_ExpQ_res <- grad_ExpQ_wh_gamma(gamma_igk, tau_gk, xi_sk, I_i)
  grad_VarQ_res <- grad_VarQ_wh_gamma(gamma_igk, tau_gk, xi_sk, I_i)
  ## partial results for original
  ExpQ_res      <- mat2ten(ExpQ_wh(gamma_igk, tau_gk, xi_sk, I_i), K)
  VarQ_res      <- mat2ten(VarQ_wh(gamma_igk, tau_gk, xi_sk, I_i), K)

  ## Compute two parts for the gradient
  partA <- (grad_VarQ_res*ExpQ_res - 2*grad_ExpQ_res*VarQ_res)/(ExpQ_res^3)
  partB <- -(mat2ten(t(logY), K) - log(ExpQ_res) + 1/2*VarQ_res/ExpQ_res^2) *
    (2*grad_ExpQ_res/ExpQ_res - partA)

  ## Across samples
  for(i in 1:N){
    Index_s <- which(I_i == i)
    ## Part A sum across sample
    res[i, , ] <- - colSums(partA[Index_s, , ] *
                              mat2ten(t(lambda_sg[,Index_s]), K), dims = 1) / 2 -
      ## Part B sum across sample
      colSums(partB[Index_s, , ] * mat2ten(t(lambda_sg[,Index_s]), K), dims = 1) / 2
  }

  return(res)

}

#-----------------------------------------------------#
## Part 2
## grident of gamma from ExpQ_logW
#-----------------------------------------------------#
## grident of gamma from ExpQ_logW
grad_logW_gamma <- function(mu_gk, gamma_igk, tau_gk,
                            alpha_gk, beta_gk, rho_gk, I_i){
  ## Compute alpha_n and number of samples per patient
  # alpha_n  <- length(unique(I_i))/2 + alpha_gk
  alpha_n  <- length(c(I_i))/2 + alpha_gk
  nI_i     <- findnSample(I_i)
  ## denominator constant sum over patient i
  denSum <- 0
  for(i in 1:N){
    denSum <- denSum + (nI_i[i] - 1)*tau_gk^2 +
      rho_gk*nI_i[i]/(nI_i[i] + rho_gk) *
      (1/nI_i[i] * tau_gk^2 + (gamma_igk[i, , ] - mu_gk)^2)
  }
  ## Store data
  res <- array(0, c(N, nG, K))
  for(i in 1:N){
    rho_term   <- rho_gk * nI_i[i] / (nI_i[i] + rho_gk)
    res[i, , ] <- - alpha_n * (rho_term*(gamma_igk[i, , ] - mu_gk)) /
      (beta_gk + 1/2*denSum)
  }
  return(res)
}

#-----------------------------------------------------#
# Functions to update gamma                           #
#-----------------------------------------------------#

## Define total loss function
loss_gamma <- function(Curr_gamma, logY, lambda_sg, mu_gk, tau_gk,
                       xi_sk = xi_sk, alpha_gk, beta_gk, rho_gk, I_i){
  # print(head(Curr_gamma, 100))
  ## Reshape the para
  cur_gamma                 <- array(0, c(N, nG, K))
  cur_gamma[1:N, 1:nG, 1:K] <- c(Curr_gamma)
  ## Sum loss from logY and logW
  res <- ExpQ_logY(logY, gamma_igk = cur_gamma,
                   tau_gk, xi_sk, lambda_sg, I_i) +
    ExpQ_logW(mu_gk, gamma_igk = cur_gamma, tau_gk,
              alpha_gk, beta_gk, rho_gk, I_i)
  return(res)
}
## Define total gradient function
grad_gamma <- function(Curr_gamma, logY, lambda_sg, mu_gk, tau_gk,
                       xi_sk = xi_sk, alpha_gk, beta_gk, rho_gk, I_i){
  ## Reshape the para
  cur_gamma                 <- array(0, c(N, nG, K))
  cur_gamma[1:N, 1:nG, 1:K] <- Curr_gamma
  ## Sum gradient from logY and logW
  res <- grad_logY_gamma(logY, lambda_sg, gamma_igk = cur_gamma,
                         tau_gk, xi_sk, I_i) +
    grad_logW_gamma(mu_gk, gamma_igk = cur_gamma,
                    tau_gk, alpha_gk, beta_gk, rho_gk, I_i)
}

Update_gamma <- function(Curr_gamma, logY, lambda_sg, xi_sk, mu_gk, tau_gk, alpha_gk, beta_gk, rho_gk, I_i){
  ## Solve by L-BFGS-B
  opt <- optim(par = Curr_gamma, fn = loss_gamma, gr = grad_gamma,
               logY = logY, lambda_sg = lambda_sg,
               mu_gk = mu_gk, tau_gk = tau_gk, xi_sk = xi_sk,
               alpha_gk = alpha_gk, beta_gk = beta_gk,
               rho_gk = rho_gk, I_i = I_i,
               method = 'L-BFGS-B', lower = -20, upper = 20,
               control = list(fnscale = -1))
  ## Return the restuls
  return(opt)
}

#-----------------------------------------------------#
# Functions for graident of tau                       #
#-----------------------------------------------------#
## Part 1
## grident of tau from ExpQ_logY
#-----------------------------------------------------#
## gradient of ExpQ_wh
grad_ExpQ_wh_tau <- function(gamma_igk, tau_gk, xi_sk, I_i){
  ## Store data
  res <- array(0, dim = c(nS, nG, K))
  ## Compute normalized xi
  xi_sk_tilde <- tilde_xi(xi_sk = xi_sk)
  ## Compute sum xi_i
  xi_s0 <- xi_sum(xi_sk = xi_sk)

  ## Across samples
  for(s in 1:nS){
    i <- I_i[s]
    res[s, , ] <- tau_gk*exp(gamma_igk[i, ,] + tau_gk^2/2) *
      vec2mat(xi_sk_tilde[s, ], nG)
  }
  # res <-
  return(res)
}

## gradient of Var_wh for cell type k
grad_VarQ_wh_k_tau <- function(gamma_gk, tau_gk, xi_k_tilde, xi_0){
  tmp1 <- 4*tau_gk*exp(2*gamma_gk + 2*tau_gk^2)
  tmp2 <- vec2mat(xi_k_tilde*(1-xi_k_tilde)/(xi_0 + 1) + xi_k_tilde^2, nG)
  tmp3 <- 2*tau_gk*exp(2*gamma_gk + tau_gk^2)*vec2mat(xi_k_tilde^2, nG)
  return(tmp1*tmp2 - tmp3)
}

## gradient of Cov_wh k != l
grad_CovQ_wh_tau <- function(gamma_gk, tau_gk, xi_k_tilde, xi_0,
                             if.tau = TRUE){

  res <- do.call(rbind, lapply(X = 1:K, FUN = Cross_Prod_ExpQ_Cov,
                               gamma_gk = gamma_gk, tau_gk = tau_gk,
                               xi_k_tilde = xi_k_tilde, xi_0 = xi_0,
                               if.tau = if.tau, if.xi_lk = FALSE,
                               if.xi_k = FALSE))
  return(t(res))
}

## gradient of Var_wh: nS, nG, K
grad_VarQ_wh_tau <- function(gamma_igk, tau_gk, xi_sk, I_i){
  ## Store data
  res <- array(0, dim = c(nS, nG, K))
  ## Compute normalized xi
  xi_sk_tilde <- tilde_xi(xi_sk = xi_sk)
  ## Compute sum xi_i
  xi_s0 <- xi_sum(xi_sk = xi_sk)

  ## Across samples
  for(s in 1:nS){
    ## Patient ID
    i <- I_i[s]
    ## Sum across cell types k
    res[s, , ] <- grad_VarQ_wh_k_tau(gamma_gk = gamma_igk[i, , ], tau_gk = tau_gk,
                                     xi_k_tilde = xi_sk_tilde[s, ], xi_0 = xi_s0[s])  +
      ## Sum across cell type l != k
      grad_CovQ_wh_tau(gamma_gk = gamma_igk[i, , ], tau_gk = tau_gk,
                       xi_k_tilde = xi_sk_tilde[s, ], xi_0 = xi_s0[s])
  }
  return(res)
}

## grident of tau from ExpQ_logY
grad_logY_tau <- function(logY, lambda_sg, gamma_igk, tau_gk, xi_sk, I_i){
  ## Store data
  res <- array(0, dim = c(nG, K))

  ## partial results for gradient
  grad_ExpQ_res <- grad_ExpQ_wh_tau(gamma_igk, tau_gk, xi_sk, I_i)
  grad_VarQ_res <- grad_VarQ_wh_tau(gamma_igk, tau_gk, xi_sk, I_i)
  ## partial results for original
  ExpQ_res      <- mat2ten(ExpQ_wh(gamma_igk, tau_gk, xi_sk, I_i), K)
  VarQ_res      <- mat2ten(VarQ_wh(gamma_igk, tau_gk, xi_sk, I_i), K)

  ## Compute two parts for the gradient
  partA <- (grad_VarQ_res*ExpQ_res - 2*grad_ExpQ_res*VarQ_res)/(ExpQ_res^3)
  partB <- -(mat2ten(t(logY), K) - log(ExpQ_res) + 1/2*VarQ_res/ExpQ_res^2) *
    (2*grad_ExpQ_res/ExpQ_res - partA)

  ## Across samples
  res <- - colSums(partA * mat2ten(t(lambda_sg), K), dims = 1) / 2 -
    ## Part B sum across sample
    colSums(partB * mat2ten(t(lambda_sg), K), dims = 1) / 2

  return(res)

}

#-----------------------------------------------------#
## Part 2
## grident of tau from ExpQ_logW
#-----------------------------------------------------#
grad_logW_tau <- function(mu_gk, gamma_igk, tau_gk, alpha_gk, beta_gk, rho_gk, I_i){
  ## Compute alpha_n and number of samples per patient
  alpha_n  <- length(c(I_i))/2 + alpha_gk
  nI_i     <- findnSample(I_i)
  ## denominator constant sum over patient i
  denSum <- 0; numSum <- 0
  for(i in 1:N){
    numSum <- numSum + (nI_i[i] - 1)*tau_gk + rho_gk /
      (nI_i[i] + rho_gk) * tau_gk
    denSum <- denSum + (nI_i[i] - 1)*tau_gk^2 +
      rho_gk*nI_i[i]/(nI_i[i] + rho_gk) *
      (1/nI_i[i] * tau_gk^2 + (gamma_igk[i, , ] - mu_gk)^2)
  }
  ## Across genes and cell types
  res <- -alpha_n * numSum / (beta_gk + 1/2 * denSum)
  return(res)
}

#-----------------------------------------------------#
## Part 3
## grident of tau from ExpQ_QlogW
#-----------------------------------------------------#
grad_QlogW_tau <- function(tau_gk, I_i){
  ## Compute number of samples per patient
  nI_i     <- findnSample(I_i)
  ## Compute the gradient
  res <- 0
  for(i in 1:N){
    res <- res - 1/2*nI_i[i]/(pi*tau_gk)
  }
  return(res)
}


#-----------------------------------------------------#
# Functions to update tau                             #
#-----------------------------------------------------#

## Define total loss function
loss_tau <- function(Curr_tau, logY, lambda_sg, mu_gk, gamma_igk,
                     xi_sk, alpha_gk, beta_gk, rho_gk, I_i){
  # print(head(Curr_gamma, 100))
  ## Reshape the para
  curr_tau            <- array(0, c(nG, K))
  curr_tau[1:nG, 1:K] <- c(Curr_tau)
  ## Sum loss from logY and logW
  res <- ExpQ_logY(logY, gamma_igk, tau_gk = curr_tau, xi_sk, lambda_sg, I_i) +
    ExpQ_logW(mu_gk, gamma_igk, tau_gk = curr_tau, alpha_gk, beta_gk, rho_gk, I_i) -
    ExpQ_QlogW(tau_gk = curr_tau)
  return(res)
}

## Define total gradient function
grad_tau <- function(Curr_tau, logY, lambda_sg, mu_gk, gamma_igk,
                     xi_sk, alpha_gk, beta_gk, rho_gk, I_i){
  ## Reshape the para
  curr_tau            <- array(0, c(nG, K))
  curr_tau[1:nG, 1:K] <- Curr_tau
  ## Sum gradient from logY and logW
  res <- grad_logY_tau(logY, lambda_sg, gamma_igk, tau_gk = curr_tau, xi_sk, I_i)  +
    grad_logW_tau(mu_gk, gamma_igk, tau_gk = curr_tau, alpha_gk, beta_gk, rho_gk, I_i) -
    grad_QlogW_tau(tau_gk = curr_tau, I_i)
}

Update_tau <- function(Curr_tau, logY, lambda_sg, mu_gk, gamma_igk,
                       xi_sk, alpha_gk, beta_gk, rho_gk, I_i){
  ## Solve by L-BFGS-B
  opt <- optim(par = Curr_tau, fn = loss_tau, gr = grad_tau,
               logY = logY, lambda_sg = lambda_sg,
               mu_gk = mu_gk, gamma_igk = gamma_igk,
               xi_sk = xi_sk, I_i = I_i,
               alpha_gk = alpha_gk, beta_gk = beta_gk,
               rho_gk = rho_gk,
               method = 'L-BFGS-B', lower = 1e-2, upper = 10,
               control = list(fnscale = -1))
  ## Return the restuls
  return(opt)
}

#-----------------------------------------------------#
# Functions for graident of xi                        #
#-----------------------------------------------------#
## Part 1
## grident of xi from ExpQ_logY
#-----------------------------------------------------#

## gradient of ExpQ_wh
grad_ExpQ_wh_xi <- function(gamma_igk, tau_gk, xi_sk, I_i){
  ## Store data
  res <- array(0, dim = c(nS, nG, K))
  ## Compute normalized xi
  xi_sk_tilde <- tilde_xi(xi_sk = xi_sk)
  ## Compute sum xi_i
  xi_s0 <- xi_sum(xi_sk = xi_sk)

  ## Across samples
  for(s in 1:nS){
    i <- I_i[s]
    res[s, , ] <- exp(gamma_igk[i, ,] + tau_gk^2/2) /
      int2ten(xi_s0[s], c(nG, K)) -
      t(vec2mat(rowSums(exp(gamma_igk[i, ,] + tau_gk^2/2) *
                          vec2mat(xi_sk[s, ], nG) / int2ten(xi_s0[s], c(nG, K))^2), K))
  }
  # res <-
  return(res)
}

## gradient of var_wh for cell type k w.r.t xi
grad_VarQ_wh_xi_k <- function(gamma_gk, tau_gk, xi_k, xi_0){
  tmp1 <- exp(2*gamma_gk + 2*tau_gk^2)
  tmp2 <- xi_0*(xi_0 + 1)*(xi_0 - 2*xi_k)/(xi_0^3*(xi_0 + 1)^2) +
    2*xi_k*xi_0/xi_0^3
  tmp2 <- vec2mat(tmp2, nG)
  tmp3 <- 2*exp(2*gamma_gk + tau_gk^2)*vec2mat(xi_k*xi_0/xi_0^3, nG)
  return(tmp1*tmp2 - tmp3)
}

## gradient of Var_wh across cell type k
grad_VarQ_wh_k_xi <- function(gamma_gk, tau_gk, xi_k, xi_0){
  tmp1 <- exp(2*gamma_gk + 2*tau_gk^2)
  tmp2 <- (xi_0*(xi_0 + 1)*xi_k - (3*xi_0 + 2)*xi_k*(xi_0 - xi_k)) /
    (xi_0^3*(xi_0 + 1)^2) - 2*xi_k^2/xi_0^3
  tmp2 <- vec2mat(tmp2, nG)
  tmp3 <- 2*exp(2*gamma_gk + tau_gk^2)*vec2mat(xi_k^2/xi_0^3, nG)
  return(rowSums(tmp1*tmp2 + tmp3))
}

## gradient of Cov_wh k != l
grad_CovQ_wh_xi <- function(gamma_gk, tau_gk, xi_k, xi_0){
  ## Compute the xi_k tilde
  xi_k_tilde = xi_k/sum(xi_k)
  ## Graident sum k, l with k != l
  res_lk <- do.call(rbind, lapply(X = 1:K, FUN = Cross_Prod_ExpQ_Cov,
                                  gamma_gk = gamma_gk, tau_gk = tau_gk,
                                  xi_k_tilde = xi_k_tilde, xi_0 = xi_0,
                                  if.tau = FALSE, if.xi_lk = TRUE,
                                  if.xi_k = FALSE))
  res_lk <- colSums(res_lk)
  ## Graident sum k with k != l
  res_k <- do.call(rbind, lapply(X = 1:K, FUN = Cross_Prod_ExpQ_Cov,
                                 gamma_gk = gamma_gk, tau_gk = tau_gk,
                                 xi_k_tilde = xi_k_tilde, xi_0 = xi_0,
                                 if.tau = FALSE, if.xi_lk = FALSE,
                                 if.xi_k = TRUE))

  res   <-  vec2mat(res_lk, K) - res_k

  return(t(res))
}

## gradient of Var_wh: nS, nG, K
grad_VarQ_wh_xi <- function(gamma_igk, tau_gk, xi_sk, I_i){
  ## Store data
  res <- array(0, dim = c(nS, nG, K))
  ## Compute sum xi_i
  xi_s0 <- xi_sum(xi_sk = xi_sk)

  ## Across samples
  for(s in 1:nS){
    ## Patient ID
    i <- I_i[s]
    ## Sum across cell types k
    res[s, , ] <- grad_VarQ_wh_xi_k(gamma_gk = gamma_igk[i,,], tau_gk = tau_gk,
                                    xi_k = xi_sk[s, ], xi_0 = xi_s0[s]) +
      grad_VarQ_wh_k_xi(gamma_gk = gamma_igk[i,,], tau_gk = tau_gk,
                        xi_k = xi_sk[s, ], xi_0 = xi_s0[s]) +
      ## Sum across cell type l != k
      grad_CovQ_wh_xi(gamma_gk = gamma_igk[i,,], tau_gk = tau_gk,
                      xi_k = xi_sk[s, ], xi_0 = xi_s0[s])
  }
  return(res)
}

## grident of gamma from ExpQ_logY
grad_logY_xi <- function(logY, lambda_sg, gamma_igk, tau_gk, xi_sk, I_i){
  ## Store data
  res <- array(0, dim = c(nS, K))

  ## partial results for gradient
  grad_ExpQ_res <- grad_ExpQ_wh_xi(gamma_igk, tau_gk, xi_sk, I_i)
  grad_VarQ_res <- grad_VarQ_wh_xi(gamma_igk, tau_gk, xi_sk, I_i)
  ## partial results for original
  ExpQ_res      <- mat2ten(ExpQ_wh(gamma_igk, tau_gk, xi_sk, I_i), K)
  VarQ_res      <- mat2ten(VarQ_wh(gamma_igk, tau_gk, xi_sk, I_i), K)

  ## Compute two parts for the gradient
  partA <- (grad_VarQ_res*ExpQ_res - 2*grad_ExpQ_res*VarQ_res)/(ExpQ_res^3)
  partB <- -(mat2ten(t(logY), K) - log(ExpQ_res) + 1/2*VarQ_res/ExpQ_res^2) *
    (2*grad_ExpQ_res/ExpQ_res - partA)

  ## Across samples
  for(s in 1:nS){
    ## Part A and B sum across genes
    res[s, ] <- - colSums((partA[s, , ] + partB[s, , ]) * t(vec2mat(lambda_sg[,s], K))) / 2
  }

  return(res)

}

#-----------------------------------------------------#
## Part 2
## grident of xi from ExpQ_logH
#-----------------------------------------------------#
## grident of xi from ExpQ_logW
grad_logH_xi <- function(xi_sk, C_i, pi_k, I_i){
  res <- array(0, c(nS, K))
  for(s in 1:nS){
    i        <- I_i[s]
    res[s, ] <- (C_i[i] * pi_k - 1)*psigamma(xi_sk[s,], 1) -
      sum((C_i[i] * pi_k - 1)*psigamma(sum(xi_sk[s, ]), 1))
  }
  return(res)
}

#-----------------------------------------------------#
## Part 3
## grident of xi from ExpQ_logQH
#-----------------------------------------------------#
## grident of xi from ExpQ_logQH
grad_logQH_xi <- function(xi_sk){
  res <- array(0, c(nS, K))
  for(s in 1:nS){
    res[s, ] <- (xi_sk[s, ] - 1)*psigamma(xi_sk[s,], 1) -
      sum((xi_sk[s, ] - 1) * psigamma(sum(xi_sk[s,]), 1))
  }
  return(res)
}


#-----------------------------------------------------#
# Functions to update xi                              #
#-----------------------------------------------------#

## Define total loss function
loss_xi <- function(Curr_xi, logY, lambda_sg, gamma_igk, tau_gk, pi_k, C_i, I_i){
  ## Reshape the para
  curr_xi            <- array(0, c(nS, K))
  curr_xi[1:nS, 1:K] <- c(Curr_xi)
  ## Sum loss from logY and logW
  res <- ExpQ_logY(logY, gamma_igk, tau_gk, xi_sk = curr_xi, lambda_sg, I_i) +
    ExpQ_logH(xi_sk = curr_xi, C_i, pi_k, I_i) -
    ExpQ_QlogH(xi_sk = curr_xi)

  if(!is.finite(res)){
    cat('Likelihood infinite \n')
    cat('logY: ', ExpQ_logY(logY, gamma_igk, tau_gk, xi_sk = curr_xi, lambda_sg, I_i), '\n',
        'logH: ', ExpQ_logH(xi_sk = curr_xi, C_i, pi_k, I_i), '\n',
        'QlogH: ',ExpQ_QlogH(xi_sk = curr_xi), '\n')
    res <- -1e10
  }
  return(res)
}

## Define total gradient function
grad_xi <- function(Curr_xi, logY, lambda_sg, gamma_igk, tau_gk, pi_k, C_i, I_i){
  ## Reshape the para
  curr_xi            <- array(0, c(nS, K))
  curr_xi[1:nS, 1:K] <- c(Curr_xi)
  ## Sum gradient from logY and logW
  res <- grad_logY_xi(logY, lambda_sg, gamma_igk, tau_gk, xi_sk = curr_xi, I_i)  +
    grad_logH_xi(xi_sk = curr_xi, C_i, pi_k, I_i) -
    grad_logQH_xi(xi_sk = curr_xi)
  if(sum(!is.finite(res)) > 0){
    cat('Gradient infinite \n')
    cat('logY: ', grad_logY_xi(logY, lambda_sg, gamma_igk, tau_gk, xi_sk = curr_xi, I_i) , '\n',
        'logH: ', grad_logH_xi(xi_sk = curr_xi, C_i, pi_k, I_i), '\n',
        'QlogH: ',grad_logQH_xi(xi_sk = curr_xi), '\n')
    res[which(!is.finite(res), arr.ind = T)] <- 0
  }
  return(res)
}


Update_xi <- function(Curr_xi, logY, lambda_sg, gamma_igk, tau_gk, pi_k, C_i, I_i){
  ## Solve by L-BFGS-B
  maxXi = ifelse(K > 6, 5, 20)
  opt <- optim(par = Curr_xi, fn = loss_xi, gr = grad_xi,
               logY = logY, lambda_sg = lambda_sg,
               gamma_igk = gamma_igk, tau_gk = tau_gk,
               pi_k = pi_k, C_i = C_i, I_i = I_i,
               method = 'L-BFGS-B', lower = 1e-6, upper = maxXi,
               control = list(fnscale = -1))
  ## Return the restuls
  return(opt)
}

#-----------------------------------------------------#
# Functions to update C_i and pi_k                    #
#-----------------------------------------------------#

Update_pik <- function(xi_sk){
  return(colMeans(xi_sk))
}

Update_Ci <- function(xi_sk, pi_k, I_i){
  res  <- rep(0, N)
  ## compute the ITH-score acorss patient
  for(i in 1:N){
    sampleIndx <- which(I_i == i)
    poolVar    <- apply(X = xi_sk[sampleIndx, ], MARGIN = 2, FUN = var)
    res[i]     <- sum(pi_k*(1-pi_k)/sum(poolVar)) #- 1
    # cov_hsk    <- cov(xi_sk[sampleIndx], )
  }
  return(res)
}
#
# Update_Ci <- function(xi_sk, pi_k, I_i){
#   res <- rep(0, N)
#   ## compute the ITH-score acorss patient
#   for(i in 1:N){
#     sampleIndx <- which(I_i == i)
#     # poolVar    <- apply(X = xi_sk[sampleIndx, ], MARGIN = 2, FUN = var)
#     # res[i]     <- sum(pi_k*(1-pi_k)/sum(poolVar)) #- 1
#     # res[i]     <- mean(rowSums(xi_sk[sampleIndx, ]))
#     cov_hsk    <- cov(xi_sk[sampleIndx, ])
#     tmp <- 0
#     for(j in 1:K){
#       for(l in j:K){
#         tmp <- tmp + (as.numeric(j == l)*pi_k[j] - pi_k[j]*pi_k[l])
#       }
#     }
#     res[i] = tmp / sum(cov_hsk[upper.tri(cov_hsk, diag = T)])
#     # res[i] = tmp / sum(diag(cov_hsk))
#   }
#   return(res)
# }

#-----------------------------------------------------#
# Iteratively update variational parameters           #
# with initialized valurs                             #
#-----------------------------------------------------#

## Check convergence
IFconverge <- function(obj_vec, t, eps){
  if(t == 1){
    return(FALSE)
  }else{
    if(abs(obj_vec[t] - obj_vec[t - 1]) < eps){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
}

#-----------------------------------------------------#
# Initialize cell-type proportions by SVM             #
#-----------------------------------------------------#

DoSVM <- function(Y, W,
                  nu.v = c(0.25, 0.5, 0.75)){

  requireNamespace('e1071')

  est.lm <- list()
  nui <- 1

  for (nu in nu.v) {

    est.m <- matrix(nrow = ncol(Y), ncol = ncol(W))

    for (s in 1:ncol(Y)) {
      svm.o <- e1071::svm(x = W,
                          y = Y[, s],
                          scale = TRUE,
                          type = "nu-regression",
                          kernel = "linear",
                          nu = nu)
      coef.v <- t(svm.o$coefs) %*% svm.o$SV
      coef.v[which(coef.v < 0)] <- 0
      total <- sum(coef.v)
      coef.v <- coef.v / total
      est.m[s, ] <- coef.v
    }
    est.lm[[nui]] <- est.m
    nui <- nui + 1
  }

  H = matrix(nrow = ncol(Y), ncol = ncol(W))

  #### select best nu
  rmse.m <- matrix(NA, nrow = ncol(Y), ncol = length(nu.v))
  for (nui in 1:length(nu.v)) {
    reconst.m <- W %*% t(est.lm[[nui]])
    for (s in 1:ncol(Y)) {
      rmse.m[s, nui] <- sqrt(mean((Y[, s] - reconst.m[, s]) ^ 2))
    }
  }
  nu.idx <- apply(rmse.m, 1, which.min)
  H <- est.m
  for (s in 1:ncol(Y)) {
    H[s, ] <- est.lm[[nu.idx[s]]][s, ]
  }
  return(H)
}


#' @title Update function for all parameters.
#'
#' @description This function is designed to estimate the cell-type
#' specific proportions and measure intra-tumor heterogeneity from
#' multiple mixed samples from patients.
#'
#' @param logY log-trans input matrix.
#' @param lambda_sg observed variability.
#' @param xi_sk0 initialized cell-type proportions.
#' @param tau_gk0 initialized cell-type variability.
#' @param gamma_igk0 initialized patient-specific mean expression.
#' @param mu_gk cell-type-specific mean expression.
#' @param alpha_gk gamma parameter for variability.
#' @param beta_gk gamma parameter for variability.
#' @param rho_gk gamma parameter for variability.
#' @param I_i sample index for patient
#' @param C_i0 initialized ITH-S.
#' @param pi_k0 nitialized mean cell-type proportions.
#' @param eps convergence criteria.
#' @param maxIters the maximum number of iterations used.
#'
#' @return res A list of estimated parameters across iterations.
#'
#' @author Peng Yang, Ziyi Li.
#'
#' @export
update_Allpara <- function(logY, lambda_sg,
                           xi_sk0, tau_gk0, gamma_igk0,
                           mu_gk, alpha_gk, beta_gk, rho_gk,
                           I_i, C_i0, pi_k0,
                           eps = 1e-4, maxIters = 1e1){

  ## Initiated values for varitional parameters
  gamma_curr <- gamma_igk0
  tau_curr   <- tau_gk0
  xi_curr    <- xi_sk0
  ## Initiated values for common paramters
  Ci_curr    <- C_i0
  pik_curr   <- pi_k0

  ## Store the objections
  obj_xi <- obj_gamma <- obj_tau <- obj_tot <- c()

  ## Store paramters
  iter_xi <- iter_gamma <- iter_tau <- list()
  iter_pik <- iter_Ci <- list()

  ## Iteratively update all the paramters
  for(t in 1:maxIters){
    # t = 1
    ## Update xi
    xi_new <- Update_xi(Curr_xi = xi_curr,
                        gamma_igk = gamma_curr, tau_gk = tau_curr,
                        pi_k = pik_curr,  C_i = Ci_curr,
                        logY = logY, lambda_sg = lambda_sg,
                        I_i = I_i)
    ## Update xi and store the objection
    xi_curr <- xi_new$par
    obj_xi  <- c(obj_xi, xi_new$value)

    ## Update gamma
    gamma_new <- Update_gamma(Curr_gamma = gamma_curr,
                              tau_gk = tau_curr, xi_sk = xi_curr,
                              # pi_k = pik_curr, C_i = Ci_curr,
                              logY = logY, lambda_sg = lambda_sg,
                              mu_gk = mu_gk, rho_gk = rho_gk,
                              alpha_gk = alpha_gk, beta_gk = beta_gk,
                              I_i = I_i)
    ## Update gamma and store the objection
    gamma_curr <- gamma_new$par
    obj_gamma  <- c(obj_gamma, gamma_new$value)



    ## Update tau
    tau_new <- Update_tau(Curr_tau = tau_curr,
                          gamma_igk = gamma_curr, xi_sk = xi_curr,
                          logY = logY, lambda_sg = lambda_sg,
                          mu_gk = mu_gk, alpha_gk = alpha_gk,
                          beta_gk = beta_gk, rho_gk = rho_gk,
                          I_i = I_i)
    ## Update gamma and store the objection
    tau_curr <- tau_new$par
    obj_tau  <- c(obj_tau, tau_new$value)

    ## Update pi_k
    # pik_curr <- Update_pik(xi_sk = xi_curr/rowSums(xi_curr))

    ## Update C_i
    # Ci_curr <- Update_Ci(xi_sk = xi_curr,
    #                      pi_k = pik_curr, I_i = I_i)

    ## Store the updating parameters
    iter_xi[[t]]    <- xi_curr
    iter_gamma[[t]] <- gamma_curr
    iter_tau[[t]]   <- tau_curr
    iter_pik[[t]]   <- pik_curr
    iter_Ci[[t]]    <- Ci_curr

    ## Total objective
    obj_tot <- c(obj_tot, obj_xi[t] + obj_gamma[t] + obj_tau[t])

    ## Print objective value for each iteration
    cat('Iteration:', t, ' objective value: ', obj_tot[t], '\n')

    ## Check convergence
    if(IFconverge(obj_tot, t, eps)){
      cat('Converged', '\n')
      break
    }
  }

  return(list(iter_xi    = iter_xi,
              iter_gamma = iter_gamma,
              iter_tau   = iter_tau,
              iter_pik   = iter_pik,
              iter_Ci    = iter_Ci,
              iter_obj   = obj_tot))

}

#' @title Estimates cell-type-specific proportions and measure intra-tumor
#' heterogeneity.
#'
#' @description This function is designed to estimate the cell-type
#' specific proportions and measure intra-tumor heterogeneity from
#' multiple mixed samples from patients.
#'
#' @param Y A mixed gene expression matrix. It is a \eqn{G}
#' by \eqn{S} matrix, where \eqn{G} is the number of genes and \eqn{S}
#' is the number of samples for mixed gene expression matrix.
#' @param reference A list contains estimated cell-type-specific mean
#' expression and variability.
#' @param sampIndex A vector to label patient samples.
#' It is a \eqn{S} by 1 vector.
#' @param rho A scale factor. The default is 1.
#' @param alpha_init A hyper-parameter. The default is 1.
#' @param eps The convergence criterion. The default is 10^(-4).
#' @param maxIters The maximum number of iterations for optimization. The default
#' is 50.
#'
#' @return A list contains the estimated outcomes
#' \itemize{
#'     \item CT - A matrix of estimated cell-type-specific proportions.
#'     It is a \eqn{S} by \eqn{K} matrix.
#'     \item ITH - A vector of estimated intratumor heterogeneity class.
#'     It is a \eqn{N} by 1 vector, where \eqn{N} is the number of patients.
#'     \item meta.data - Metadata over iterations.
#' }
#'
#' @author Peng Yang, Ziyi Li.
#'
#' @export
ICeITH <- function(Y, reference, sampIndex, rho = 1, alpha_init = 1,
                   eps = 1e-4, maxIters = 50){

  ## Reference
  mu_gk     <- reference$mu_gk
  lambda_gk <- 1/reference$var_gk

  ## Check input values

  if(any(duplicated(colnames(mu_gk))) | any(duplicated(colnames(lambda_gk)))){
    stop("reference matrix has duplicated colnames.")
  }

  if(any(colnames(mu_gk) == "")){
    stop("colnames of reference matrix are cell type labels and cannot be empty.")
  }

  if(nrow(Y)!=nrow(mu_gk)){
    stop("matrix Y and reference do not have the same number of genes.")
  }

  if(!identical(rownames(Y), rownames(mu_gk))){
    stop("gene names of Y and reference are different.")
  }

  if(any(is.na(Y))){
    stop("matrix Y must not contain any NA entries.")
  }

  if(any(Y < 0)){
    stop("matrix Y should be non-negative.")
  }

  ## Log-transformation of Y
  if(if.zero(Y)){
    logY = log(Y)
  }else{
    message("Adding 1 to Y to ensure valid log transformation.")
    logY = log(Y + 1)
  }

  # nS   = dim(Y)[2]                  ## # of tumor samples
  # nG   = dim(Y)[1]                  ## # of genes
  # K    = dim(mu_gk)[2]              ## # of cell types
  # N    = length(unique(sampIndex))  ## # of patients
  assign('nS', dim(Y)[2],      envir = .GlobalEnv)
  assign('nG', dim(Y)[1],      envir = .GlobalEnv)
  assign('K',  dim(mu_gk)[2],  envir = .GlobalEnv)
  assign('N',  length(unique(sampIndex)), envir = .GlobalEnv)

  ## Estimating sample-wsie variablility
  lambda_sg <- array(data = sampleVar(logY),
                     dim = c(nG, nS),
                     dimnames = list(paste('Gene', 1:nG),
                                     paste('Sample', 1:nS)))

  ## Choose hyper-parameters
  alpha_gk <- array(data = alpha_init,
                    dim = c(nG, K),
                    dimnames = list(paste('Gene', 1:nG),
                                    paste('CT', 1:K)))
  beta_gk  <- alpha_gk * reference$var_gk

  rho_gk <- array(data = rho,
                  dim = c(nG, K),
                  dimnames = list(paste('Gene', 1:nG),
                                  paste('CT', 1:K)))

  ## Initiate cell type proportions
  # pi_k0 <- rdirichlet(n = 1, alpha = 2*seq(1:K))
  # h_sk0 <- rdirichlet(n = nS, alpha = 10*pi_k0)
  # h_sk0 <- array(0, c(nS, K))
  # for(s in 1:nS){
  #   h_sk0[s,] <- abs(lm(logY[,s] ~ mu_gk - 1)$coeff)
  # }
  h_sk0 = DoSVM(Y = Y, W = exp(mu_gk + 1/lambda_gk/2))
  h_sk0 = h_sk0/rowSums(h_sk0)
  ## Initiate mean cell type proportions
  pi_k0 <- colMeans(h_sk0)

  ## Initiate ITH-S
  C_i0 <- Update_Ci(xi_sk = h_sk0, pi_k = pi_k0, I_i = sampIndex)

  ## Mean expression for patient/gene/cell-type-specific mean expression
  gamma_igk0 <- array(data = 0,
                      dim = c(N, nG, K))
  for(i in 1:N) gamma_igk0[i, , ] <- mu_gk

  ## Variability for patient/gene/cell-type-specific mean expression
  tau_gk0 <- 1/lambda_gk

  ## Estimating varitional parameters
  res_allpara <- update_Allpara(logY, lambda_sg,
                                h_sk0, tau_gk0, gamma_igk0,
                                mu_gk, alpha_gk, beta_gk, rho_gk,
                                sampIndex, C_i0, pi_k0,
                                eps = eps, maxIters = maxIters)

  ## Estimated cell type proportions
  CT        <- res_allpara$iter_xi[[length(res_allpara$iter_xi)]] /
    rowSums(res_allpara$iter_xi[[length(res_allpara$iter_xi)]])
  ## ITH class
  ITH.class <- ifelse(res_allpara$iter_Ci[[length(res_allpara$iter_Ci)]] >
                     median(res_allpara$iter_Ci[[length(res_allpara$iter_Ci)]]),
                               'Low ITH', 'High ITH')
  ## Return the list
  res   <- list(CT        = CT,
                ITH       = as.factor(ITH.class),
                meta.data = res_allpara)

  ## Return the results
  return(res)

}
