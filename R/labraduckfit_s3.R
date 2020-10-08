#' Create labraduckfit object 
labraduckfit <- function(D, N, T, coord_system, iter=NULL,  
                       alr_base=NULL, ilr_base=NULL,
                       Eta=NULL, Sigma=NULL, Sigma_default=NULL, 
                       Y=NULL, upsilon=NULL, Xi=NULL, Xi_default=NULL,
                       gamma=NULL, F=NULL, G=NULL, W=NULL, W_scale=NULL, M0=NULL, C0=NULL,
                       observations=NULL,
                       init=NULL){
  m <- new_labraduckfit(D, N, T, coord_system, iter, alr_base, ilr_base,
                      Eta, Sigma, Sigma_default, 
                      Y, upsilon, Xi, Xi_default,
                      gamma, F, G, W, W_scale, M0, C0,
                      observations,
                      init)
  verify(m)
  return(m)
}


# internal function 
new_labraduckfit <- function(D, N, T, coord_system, iter=NULL,  
                       alr_base=NULL, ilr_base=NULL,
                       Eta=NULL, Sigma=NULL, Sigma_default=NULL, 
                       Y=NULL, upsilon=NULL, Xi=NULL, Xi_default=NULL,
                       gamma=NULL, F=NULL, G=NULL, W=NULL, W_scale=NULL, M0=NULL, C0=NULL,
                       observations=NULL,
                       init=NULL){
  m <- new_pibblefit(D, N, Q=0, coord_system, iter, alr_base, ilr_base,
                      Eta, Lambda=NULL, Sigma, Sigma_default, 
                      Y, X=NULL, upsilon, Theta=NULL, Xi, Xi_default, Gamma=NULL, 
                      init, names_categories=NULL, names_samples=NULL, names_covariates=NULL)
  m$T <- T
  m$gamma <- gamma
  m$F <- F
  m$G <- G
  m$W <- W
  m$W_scale <- W_scale
  m$M0 <- M0
  m$C0 <- C0
  m$observations <- observations
  class(m) <- c("labraduckfit", "pibblefit")
}

# silent for now
verify.labraduckfit <- function(m,...){
  # verify.labraduckfit(m)
  # stopifnot(is.integer(m$P))
  # P <- m$P
  # ifnotnull(m$VCScale, check_dims(m$VCScale, P, "VCScale"))
  # ifnotnull(m$U, check_dims(m$U, c(P*m$Q, m$Q), "U"))
  # ifnotnull(m$ellinit, check_dims(m$P, P, "ellinit"))
}

# reinstate later
req.labraduckfit <- function(m, r){
  # present <- sapply(m[r], is.null)
  # if(any(present)){
  #   stop("maltipoofit object does not contain required components:", r[present])
  # }
}