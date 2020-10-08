
#' Interface to fit DLM models
#' 
#' Labraduck - exchangeable dynamic linear models in stray
#'
#' @param Y D x N matrix of counts (if NULL uses priors only)
#' @param upsilon dof for inverse wishart prior (numeric must be > D) 
#'   (default: D+3)
#' @param Xi (D-1)x(D-1) prior covariance matrix
#'   (default: ALR transform of diag(1)*(upsilon-D)/2 - this is 
#'   essentially iid on "base scale" using Aitchison terminology)
#' @param F Observation-level transformation matrix (see model)
#' @param G System-level transformation matrix (see model)
#' @param W System-level noise covariance (of dim Q) over rows of eta transpose (see model)
#' @param scale_init Scale initialization system (W) and observation level (gamma) noise
#' @param M0 Prior mean over system matrix
#' @param C0 Prior covariance over system-level noise
#' @param observations Vector timepoints indicating offset from first observation (labeled 1)
#'   (ex: 1, 2, 3, 6, 7, 8 where there is a gap in observations)
#' @param init (D-1) x Q initialization for Eta for optimization 
#' @param pars character vector of posterior parameters to return
#' @param m object of class labraduckfit 
#' @param ... 
#' 
#' @details the full model is given by:
#'    \deqn{Y_j \sim Multinomial(Pi_j)} 
#'    \deqn{Pi_j = Phi^{-1}(Eta_j)}
#'    \deqn{Eta_t^T = F_t^T Theta_t + v_t^T, v_t \sim N_{D-1}(0, gamma_t Sigma)}
#'    \deqn{Theta_t^T = G_t Theta_{t-1} + w_t, w_t \sim N_{Q x D-1}(0, W_t, Sigma)}
#'    \deqn{Sigma \sim InvWish(upsilon, Xi)}
#' @return an object of class labraduckfit
#' @md
#' @name labraduck_fit
NULL

#' @rdname labraduck_fit
#' @export
labraduck <- function(Y=NULL, upsilon=NULL, Xi=NULL, F, G, W, M0, C0, observations, gamma_scale=0, W_scale=0, 
                    init=NULL, pars=c("Eta", "Sigma", "Thetas_filtered", "Thetas_smoothed", "M_star", "Eta_DLM"),
                    ...){
  args <- list(...)
  N <- try_set_dims(c(ncol(Y), args[["N"]]))
  D <- try_set_dims(c(nrow(Y), nrow(Xi)+1, ncol(Xi)+1, args[["D"]]))
  T <- max(observations) # presumed integers!
  if (any(c(N, D) <=0)) stop("N and D must all be greater than 0 (D must be greater than 1)")
  if (D <= 1) stop("D must be greater than 1")
  
  ## construct default values ##
  # for priors
  if (is.null(upsilon)) upsilon <- D+3  # default is minimal information 
                                        # but with defined mean
  if (is.null(Xi)) {
    # default is iid on base scale
    # G <- cbind(diag(D-1), -1) ## alr log-constrast matrix
    # Xi <- 0.5*G%*%diag(D)%*%t(G) ## default is iid on base scale
    Xi <- matrix(0.5, D-1, D-1) # same as commented out above 2 lines
    diag(Xi) <- 1               # same as commented out above 2 lines
    Xi <- Xi*(upsilon-D) # make inverse wishart mean Xi as in previous lines 
  }
  # add defaults for the DLM parameters (TODO)
   
  # check dimensions
  check_dims(upsilon, 1, "upsilon")
  check_dims(Xi, c(D-1, D-1), "Xi")
  
  # set number of iterations 
  n_samples <- args_null("n_samples", args, 2000)
  use_names <- args_null("use_names", args, TRUE)
  
  # This is the signal to sample the prior only
  if (is.null(Y)){
    # create pibblefit object and pass to sample_prior then return
    # untested (TODO)
    if(gamma_scale > 0) {
      gamma_init = gamma_scale
    } else {
      gamma_init = 1
    }
    if(W_scale > 0) {
      W_scale_init = W_scale
    } else {
      W_scale_init = 1
    }
    out <- labraduckfit(N=N, D=D, T=T, coord_system="alr", alr_base=D, upsilon=upsilon, Xi=Xi,
      gamma=gamma_init, F=F, G=G, W=W, W_scale=W_scale_init, M0=M0, C0=C0, observations=observations)
    out <- sample_prior(out, n_samples=n_samples, pars=pars, use_names=use_names)
    return(out)
  } else {
    if(is.null(init)) init <- random_pibble_init(Y)   # initialize init 
  }

  # for optimization and laplace approximation
  calcGradHess <- args_null("calcGradHess", args, TRUE)
  b1 <- args_null("b1", args, 0.9)
  b2 <- args_null("b2", args, 0.99)
  step_size <- args_null("step_size", args, 0.003)
  epsilon <- args_null("epsilon", args, 10e-7)
  eps_f <- args_null("eps_f", args, 1e-10)
  eps_g <- args_null("eps_g", args, 1e-4)
  max_iter <- args_null("max_iter", args, 10000)
  verbose <- args_null("verbose", args, FALSE)
  verbose_rate <- args_null("verbose_rate", args, 10)
  decomp_method <- args_null("decomp_method", args, "cholesky")
  eigvalthresh <- args_null("eigvalthresh", args, 0)
  jitter <- args_null("jitter", args, 0)
  multDirichletBoot <- args_null("multDirichletBoot", args, -1.0)
  optim_method <- args_null("optim_method", args, "lbfgs")
  useSylv <- args_null("useSylv", args, TRUE)
  ncores <- args_null("ncores", args, -1)
  apply_smoother <- args_null("apply_smoother", args, TRUE)
  
  ## fit collapsed model ##
  fitc <- optimLabraduckCollapsed(Y, upsilon, Xi, F, G, W, M0, C0, observations, init,
                                gamma_scale, W_scale, n_samples, 
                                calcGradHess, b1, b2, step_size, epsilon, eps_f, 
                                eps_g, max_iter, verbose, verbose_rate, 
                                decomp_method, optim_method, eigvalthresh, 
                                jitter, multDirichletBoot, 
                                useSylv, ncores)

  timerc <- parse_timer_seconds(fitc$Timer)

  gamma_scale <- fitc$gamma_scale
  W_scale <- fitc$W_scale

  out <- list()

  # if n_samples=0 or if hessian fails, then use MAP eta estimate for 
  # uncollapsing and unless otherwise specified against, use only the 
  # posterior mean for Lambda and Sigma 
  if (is.null(fitc$Samples)) {
    fitc$Samples <- add_array_dim(fitc$Pars, 3)
    ret_mean <- args_null("ret_mean", args, TRUE)
    if (ret_mean && n_samples>0){
      warning("Laplace Approximation Failed, using MAP estimate of eta", 
              " to obtain Posterior mean of Lambda and Sigma", 
              " (i.e., not sampling from posterior distribution of Lambda or Sigma)")
    }
    if (!ret_mean && n_samples > 0){
      warning("Laplace Approximation Failed, using MAP estimate of eta", 
              "but ret_mean was manually specified as FALSE so sampling", 
              "from posterior of Lambda and Sigma rather than using posterior mean")
    }
  } else {
    ret_mean <- args_null("ret_mean", args, FALSE)
  }
  
  seed <- args_null("seed", args, sample(1:2^15, 1))
  ## uncollapse collapsed model ##
  if(ret_mean) {
    apply_smoother <- FALSE
  }

  cat("Loglik:",fitc$LogLik,"\n")

  # ret_mean overrides sample returning for now
  fitu <- uncollapseLabraduck(fitc$Samples, F, G, W, gamma_scale, W_scale, upsilon, Xi, M0, C0, observations, seed=seed, ret_mean=ret_mean, apply_smoother=apply_smoother, ncores=ncores)

  timeru <- parse_timer_seconds(fitu$Timer)
  
  timer <- c(timerc, timeru)
  timer <- timer[which(names(timer)!="Overall")]
  timer <- c(timer, 
             "Overall" = unname(timerc["Overall"]) +  unname(timeru["Overall"]), 
             "Uncollapse_Overall" = timeru["Overall"])
  
  # Marginal Likelihood Computation -- check this (TODO)
  # d <- D^2 + N*D + D*Q
  # logMarginalLikelihood <- fitc$LogLik+d/2*log(2*pi)+.5*fitc$logInvNegHessDet-d/2*log(N)
    
  ## pretty output ##
  out <- list()
  if ("Eta" %in% pars){
    out[["Eta"]] <- fitc$Samples
  }
  if ("Thetas_filtered" %in% pars & !ret_mean){
    out[["Thetas_filtered"]] <- fitu$Thetas_filtered_sample
  }
  if ("Thetas_smoothed" %in% pars & !ret_mean){
    out[["Thetas_smoothed"]] <- fitu$Thetas_smoothed_sample
  }
  if ("M_star" %in% pars & !ret_mean){
    out[["M_star"]] <- fitu$M_star_sample
  }
  if ("Eta_DLM" %in% pars & !ret_mean){
    out[["Eta_DLM"]] <- fitu$Eta_sample
  }
  if ("Sigma" %in% pars){
    out[["Sigma"]] <- fitu$Sigma
  }
  
  # # By default just returns all other parameters
  out$D <- D
  out$N <- N
  out$T <- T
  out$coord_system <- "alr"
  out$iter <- dim(fitc$Samples)[3]
  out$upsilonT <- fitc$upsilonT
  out$XiT <- fitc$XiT
  out$alr_base <- D
  out$Y <- Y
  out$upsilon <- upsilon
  out$gamma_scale <- gamma_scale
  out$F <- F
  out$G <- G
  out$W <- W
  out$W_scale <- W_scale
  out$M0 <- M0
  out$C0 <- C0
  out$observations <- observations
  out$init <- init
  # for other methods
  out$summary <- NULL
  out$Timer <- timer
  # out$logMarginalLikelihood <- logMarginalLikelihood
  attr(out, "class") <- c("labraduckfit", "pibblefit")
  # add names if present 

  if (use_names) out <- name(out)
  verify(out) # verify the labraduckfit object
  return(out)
}

