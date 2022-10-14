#' @importFrom gtools rdirichlet
#-------------------------------------------------------------------#
# utility functions                                                 #
#-------------------------------------------------------------------#

Identify_Closest <- function(x, pts2compare){
  return(which.min(abs(pts2compare-x)))
}

Id_Close_appVec <- function(x, pts2compare){
  out = apply(X=matrix(x,ncol=1), MARGIN = 1, FUN = Identify_Closest,
              pts2compare=pts2compare)
  return(out)
}

Id_Close_appMat <- function(X, pts2compare){
  out = apply(X = X, MARGIN = 2, FUN = Id_Close_appVec,
              pts2compare=pts2compare)
  return(out)
}

#-------------------------------------------------------------------#
# simulation functions                                              #
#-------------------------------------------------------------------#
#' @title Function to simulate multi-region gene expression data.
#'
#' @description This function is designed to generate the cell-type-specific
#' reference profile and the multi-region gene expression data to perform
#' reference-based deconvolution by ICeITH model.
#'
#' @param K Number of cell types.
#' @param G Number of marker genes.
#' @param lowS The minimal number of samples per subject.
#' @param maxS The maximal number of samples per subject.
#' @param N Number of patients.
#' @param rho_gk Scale parameter to control the information borrowing
#'                from the reference profile.
#' @param nRef Number of samples for the reference profile.
#' @param lowC Lowest intratumor heterogeneity score.
#' @param maxC Highest intratumor heterogeneity score.
#'
#' @return A list contains the simulated paramters
#' \itemize{
#'     \item Y - Mixed gene expression.
#'     \item ct_g - Cell-type-specific genes.
#'     \item ct_s - Cell-type-specific samples from the reference profile.
#'     \item I_i - Sample index to subjects.
#'     \item W_igk - Hidden cell-type-specific gene expression for each sample.
#'     \item mu_gk - True reference mean expression.
#'     \item mu_igk - True patient-specific mean expression.
#'     \item lambda_gk - True gene and cell-type-specific variability.
#'     \item h_sk - Cell-type-specific proportions.
#'     \item ITH - Patient-specific intratumor level.
#'     \item pi_k - Mean cell-type-specific proportions.
#' }
#'
#' @export
#'
sim_func <- function(K = 6,      ## Number of cell types
                     G = 250,    ## Number of genes
                     # S = 30,     ## Number of tumor samples
                     lowS = 4,   ## lowest # of samples per patient
                     maxS = 8,   ## highest # of samples per patient
                     N = 20,     ## Number of patients
                     rho_gk = 1,
                     nRef = 100, ## Number of samples for reference
                     lowC = 0.1,
                     maxC = 5){  ## Highest ITH-S

  requireNamespace('gtools')

  ## Generate the Index of tumor samples
  if(lowS != maxS){
    nS_i <- sample(lowS:maxS, size = N, replace = T)
  }else{
    nS_i <- rep(lowS, N)
  }

  I_i  <- rep(1:N, nS_i)
  S    <- length(I_i)
  ## Generate the Index for cell type in reference samples
  NCT             <- floor(nRef/K)*K
  CT_Index        <- rep(K, nRef)
  CT_Index[1:NCT] <- rep(1:K, each = floor(nRef/K))

  #---------------------------------------------------------#
  # Step 1:                                                 #
  # Simulate Purified Reference                             #
  #---------------------------------------------------------#
  ct.gene <- floor(G/K)
  ct.gene <- rep(1:K, each = ct.gene)
  ct.gene <- c(ct.gene, rep(K, G - length(ct.gene)))

  ## Randomly assign mean expression to each cell type
  NML <- runif(n = G, min = 0, max = 1)
  mu_gk <- array(data = 0,
                 dim = c(G, K),
                 dimnames = list(paste('Gene', seq(1:G)),
                                 paste('CT', seq(1:K))))

  for(g in 1:G){

    if(NML[g] < 1/3){

      mu_gk[g,]   = runif(n = K, min = 2, max = 4)
      ct          = ct.gene[g]
      mu_gk[g,ct] = runif(n = 1, 3.5, 5)

    }else if(NML[g] < 2/3){

      mu_gk[g,]   = runif(n = K, min = 4, max = 6)
      ct          = ct.gene[g]
      mu_gk[g,ct] = runif(n = 1, 5.5, 7)

    }else{

      mu_gk[g,]   = runif(n = K, min = 6, max = 8)
      ct          = ct.gene[g]
      mu_gk[g,ct] = runif(n = 1, 7.5, 9)

    }

  }

  #-------------------------------------------------#
  # Mean-Variance Relationship                      #
  #-------------------------------------------------#
  data(mean_var_relation)
  sqrt_lambda_id = array(0, dim = c(G, K))
  sqrt_lambda_id = Id_Close_appMat(X = mu_gk, pts2compare = mean_var_relation$x)
  sqrt_lambda    = array(data = 0,
                         dim = c(G, K),
                         dimnames = list(paste('Gene', seq(1:G)),
                                         paste('CT', seq(1:K))))

  for(i in 1:K){
    sqrt_lambda[,i] = mean_var_relation$y[sqrt_lambda_id[,i]]
  }

  sqrt_lambda = sqrt_lambda + rnorm(G, 0, 0.085)
  lambda_gk   = 1/(sqrt_lambda)^2


  ## Simulate mean expression for each gene and cell type belong to patient i
  mu_igk <- array(data = 0,
                  dim = c(N, G, K),
                  dimnames = list(paste('Patient', seq(1:N)),
                                  paste('Gene', seq(1:G)),
                                  paste('CT', seq(1:K))))
  # rho_gk <- 1

  for(i in 1:N){
    # i = 1
    mu_igk[i,,] <- mu_gk

    for(k in 1:K){

      mu_igk[i,,k] <- mu_igk[i,,k] + rnorm(G, mean = 0, sd = 1/sqrt(rho_gk*lambda_gk[,k]))

    }
  }

  #-------------------------------------------------#
  # Simulated Cell Type-specific Expression         #
  #-------------------------------------------------#
  ## Reference level
  X_gr <- array(data = 0,
                dim = c(G, nRef),
                dimnames = list(paste('Gene', seq(1:G)),
                                paste('CT', CT_Index)))
  for(g in 1:G){
    for(r in 1:nRef){
      ## Cell type index
      ct        <- CT_Index[r]
      ## Simulate expresssion
      X_gr[g,r] <- exp(rnorm(n = 1,
                             mean = mu_gk[g, ct],
                             sd = 1/sqrt(lambda_gk[g, ct])))
    }
  }

  #-------------------------------------------------#
  # Simulated Hidden Cell Type-specific Expression  #
  # for each patient                                #
  #-------------------------------------------------#
  W_igk <- list()
  for(i in 1:N){
    nS    <- sum(I_i == i)
    W_sgk <- array(data = 0,
                   dim = c(nS, G, K),
                   dimnames = list(paste('Sample', seq(1:nS)),
                                   paste('Gene', seq(1:G)),
                                   paste('CT', 1:K)))
    ## Simulate acorss samples from one patient
    for(g in 1:G){
      for(k in 1:K){
        W_sgk[,g,k] <- exp(rnorm(n = nS,
                                 mean = mu_igk[i,g,k],
                                 sd = 1/sqrt(lambda_gk[g,k])))
      }
    }
    W_igk[[i]] = W_sgk
  }
  names(W_igk) = paste('Patient', seq(1:N))
  #---------------------------------------------------------#
  # Step 2:                                                 #
  # Simulate cell type proportions                          #
  #---------------------------------------------------------#

  ## simulate ITH-S
  C_i  <- runif(n = N, min = lowC, max = maxC)
  # C_i <- rgamma(n = N, shape = a_0, rate = b_0)
  ## simulate pi_k
  if(K==2){
    aval = c(6, 4)
  } else if(K==3){
    aval = c(4, 2.5, 1.5)*10/8
  } else if(K==4){
    aval = c(4, 3, 2, 1)
  } else{
    aval = c(4, 3, 2, rep(1/(K-3), K-3))
  }
  pi_k <- gtools::rdirichlet(n = 1, alpha = aval)
  ## Simulate true cell type fractions
  h_sk <- list()
  for(i in 1:N){
    nS                  <- sum(I_i == i)
    h_sk[[i]]           <- rdirichlet(n = nS, alpha = C_i[i] * pi_k)
    rownames(h_sk[[i]]) <- paste('Sample', seq(1:nS))
    colnames(h_sk[[i]]) <- paste('CT', seq(1:K))
  }
  names(h_sk) = paste('Patient', seq(1:N))


  #---------------------------------------------------------#
  # Step 3:                                                 #
  # Simulate mixtured gene expression                       #
  #---------------------------------------------------------#
  Y <- array(data = 0,
             dim = c(G, S),
             dimnames = list(paste('Gene', seq(1:G)),
                             paste('Sample', seq(1:S))))
  for(i in 1:N){
    I_indx <- which(I_i == i)
    nS  <- length(I_indx)
    for(s in 1:nS){
      Y[,I_indx[s]] <- W_igk[[i]][s,,] %*% h_sk[[i]][s,]
    }
  }

  ITH <- ifelse(C_i > median(C_i), 'Low ITH', 'High ITH')

  return(list(Y = Y,                   ## Mixtured tumor samples
              X_gr = X_gr,             ## Reference samples
              ct_g = ct.gene,          ## CT-specific genes
              ct_s = CT_Index,         ## CT-specific samples
              I_i = I_i,               ## Sample Index
              W_igk = W_igk,           ## Sample-specific hidden expression
              mu_gk = mu_gk,           ## True reference expression
              mu_igk = mu_igk,         ## True patient-specific expression
              lambda_gk = lambda_gk,   ## True gene and CT-specific variability
              h_sk = h_sk,             ## Cell-type-specific proportions
              ITH = ITH,               ## ITH-S
              pi_k = pi_k))            ## Mean CT-specific proportions


}
