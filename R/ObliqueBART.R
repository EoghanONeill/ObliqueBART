

make_01_norm <- function(x) {
  a <- min(x)
  b <- max(x)
  return(function(y0) (y0 - a) / (b - a))
}


#' @title Bayesian Additive Regression Trees with Oblique Splits
#'
#' @description Bayesian Additive Regression Trees with Oblique Splits
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm' 'aggregate' 'contrasts' 'model.matrix'
#' @importFrom MCMCpack 'rdirichlet'
#' @importFrom truncnorm 'rtruncnorm'
#' @importFrom collapse fmean fsum
#' @param x Training covariate matrix. Rows correspond to observations, columns correspond to covariates used in splitting rules.
#' @param y Training outcome vector.
#' @param coef_prior Prior distribution for splitting coefficients: "univariate_normal", "multivariate_normal", or "simplex"
#' @param coef_hyperprior Hyperprior for splitting coefficients.
#' @param p_bar Parameter if coef_hyperprior is univariate_normal_fixed_binomial. Mean number of covariates included in splitting rules.
#' @return The following objects are returned:
#' \item{trees}{List (over MCMC iterations) of lists (over trees in a sum-of-trees model) of tree matrices.}
#'
#' @useDynLib ObliqueBART, .registration = TRUE
#' @export
ObliqueBART <- function(x,
                        y,
                        sparse = TRUE,
                        ntrees = 10,
                        node_min_size = 5,
                        alpha = 0.95,
                        beta = 2,
                        nu = 3,
                        lambda = 0.1,
                        mu_mu = 0,
                        sigma2 = 1,
                        sigma2_mu = 1,
                        k = 2,
                        sigquant = .90,
                        nburn = 1000,
                        npost = 1000,
                        nthin = 1,
                        penalise_num_cov = FALSE,
                        lambda_cov = 0.4,
                        nu_cov = 2,
                        coef_prior = "univariate_normal",        # Prior distribution for splitting coefficients
                        coef_hyperprior = "none",
                        coef_norm_hyperprior = "fixed", # fixed at zero, or varying
                        threshold_prior = "discrete_uniform",
                        p_bar = ceiling(ncol(x)/3),
                        norm_sigma_init = "Zellner", # if a normal prior, initialize prior variaince with Zellner's prior or ridge?
                        xi_prior = "Boojum", # only relevant if hyperprior for Dirichlet prior
                        MH_propstep_xi = 0.15,
                        norm_unit_sphere = FALSE, # only relevant for normal distribution draws
                        split_mix = FALSE, # splitting coefficients drawn from mixture distribution (mixture of vairable inclusion probabiltiies and coef values)
                        num_clust = 1 # number of clusters in the finite mixture for splitting coefficients
                        ) {



  if(!(norm_sigma_init %in% c("Zellner",
                              "ridge"))){
    stop("Currently supported norm_sigma_init options are:
         Zellner, ridge")
  }

  # check threshold prior options
  if(!(threshold_prior %in% c("discrete_uniform",
                         "continuous_min_max",
                         "continuous_minus_plus_1"))){
    stop("Currently supported threshold prior options are
         discrete_uniform, continuous_min_max, continuous_minus_plus_1")
  }

  if((threshold_prior == "continuous_minus_plus_1" ) & (coef_prior != "simplex")){
    stop("The threshold_prior of continuous_minus_plus_1 is not recommended unless
         coef_prior is simplex (and the covariates are scaled between 0 and 1),
         otherwise the linear combinations are not constrained to be
         between -1 and 1.")
  }




  # check coefficient prior options
  if(!(coef_prior %in% c("univariate_normal",
                         "multivariate_normal",
                         "simplex"))){
    stop("Currently supported splitting coefficient prior options are
         univariate_normal, multivariate_normal, simplex.")
  }


  # check coefficient hyperprior options
  # if(!is.na(coef_hyperprior)){
  if(!(coef_hyperprior %in% c("none",
                              "univariate_normal",
                              "univariate_normal_fixed_binomial",
                              "univariate_normal_betabinomial_onetheta",
                              "univariate_normal_betabinomial_theta_j",
                              "univariate_normal_betabinomial_theta_j_sigma_j",
                              "multivariate_normal",
                              "simplex_fixed_beta_binomial_theta_j",
                              "simplex_fixed_Dir_trinomial_theta_j",
                              "simplex_fixed_Dir_binomial_plusminus_theta_j",
                              "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"
  ))){
    stop("Currently supported splitting coefficient hyperprior options are:
       none, univariate_normal, multivariate_normal")
  }
  # }

  # check coefficient prior and hyperprior match
  if(coef_hyperprior != "none"){
    if(coef_prior == "univariate_normal"){
      if(!(coef_hyperprior %in% c(#"none",
                                  "univariate_normal",
                                  "univariate_normal_fixed_binomial",
                                  "univariate_normal_betabinomial_onetheta",
                                  "univariate_normal_betabinomial_theta_j",
                                  "univariate_normal_betabinomial_theta_j_sigma_j") )){
        stop("Currently only the univariate normal hyperprior is supported for a univariate normal coefficient prior")
      }
    }
    if(coef_prior == "multivariate_normal"){
      if(coef_hyperprior != "multivariate_normal"){
        stop("Currently only the multivariate normal hyperprior is supported for a multivariate normal coefficient prior")
      }
    }
    if(coef_prior == "simplex"){
      if(!(coef_hyperprior %in% c(#"none",
        "simplex_fixed_beta_binomial_theta_j",
        "simplex_fixed_Dir_trinomial_theta_j",
        "simplex_fixed_Dir_binomial_plusminus_theta_j",
        "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j") )){
        stop("Currently simplex hyperprior options are:
             none,
             simplex_fixed_beta_binomial_theta_j,
             simplex_fixed_Dir_trinomial_theta_j,
             simplex_fixed_Dir_binomial_plusminus_theta_j",
             "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j")
      }
      # stop("Currently no hyperpriors supported for simplex coefficient prior")
      # if(coef_hyperprior != "multivariate_normal"){
      #   stop("Currently only the multivariate normal hyperprior is supported for a multivariate normal coefficient prior")
      # }
    }
  }


  # Extract control parameters
  node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  # var_count = rep(0, ncol(x))
  # var_count_store = matrix(0, ncol = ncol(x), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(x), nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(x)
  s = rep(1/p, p)


  # tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))

  sigma2_mu <- (max(y_scale)-min(y_scale))/((2 * k * sqrt(ntrees))^2)

  # if(is.na(sigest)) {
    if(p < n) {
      df = data.frame(x,y_scale)
      lmf = lm(y_scale~.,df)
      sigest = summary(lmf)$sigma
    } else {
      sigest = sd(y_scale)
    }
  # }
  qchi = qchisq(1.0-sigquant,nu)
  lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior

  ##### scale covariates ############



  ecdfs   <- list()
  for(i in 1:ncol(x)) {
    ecdfs[[i]] <- ecdf(x[,i])
    if(length(unique(x[,i])) == 1) ecdfs[[i]] <- identity
    if(length(unique(x[,i])) == 2) ecdfs[[i]] <- make_01_norm(x[,i])
  }
  for(i in 1:ncol(x)) {
    x[,i] <- ecdfs[[i]](x[,i])
    # x.test[,i] <- ecdfs[[i]](x.test[,i])
  }


  ###### Initialization ################

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = y_scale,
                            X = x)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  y_hat = get_predictions(curr_trees, x, single_tree = ntrees == 1)



  # Initialize list of hyperparameters
  hyp_par_list <- list()

  if(coef_prior == "univariate_normal"){
   # sigma2_beta_scalar <- 1 # what would be a reasonable initial value?
   # maybe it would be better to have a vector of coefficient values? Scaled by a linear regression or ridge?
   # Scale covariates first?

    if(norm_sigma_init == "ridge"){
      # ridge-like prior sets prior equal to constant (another hyperpameter) times error variance
      ridge_lambda <- 1

      sigma2_beta_vec <- rep(ridge_lambda*sigma2, p)
    }
    if(norm_sigma_init == "Zellner"){
      # something like the Zellener G-prior uses the inverse sample covariance matrix of the covariates
      # multiplied by sigma
      zell_g <- 1
      zell_diag <- diag(solve( t(x) %*% x ))
      sigma2_beta_vec <- zell_g*zell_diag*sigma2
    }

    hyp_par_list$beta_bar_vec <- rep(0, p) # initialize hyperparameter mean at zero?
    hyp_par_list$sigma2_beta_vec <- sigma2_beta_vec


    if(coef_hyperprior == "univariate_normal"){
      sigma2_beta_bar <- sigma2_beta_vec # difficult to know how to set variance of hyperparameter means
    }
    if(coef_hyperprior == "univariate_normal_fixed_binomial"){
      hyp_par_list$theta <- p_bar/p
    }
    if(coef_hyperprior == "univariate_normal_betabinomial_onetheta"){
      hyp_par_list$theta <- p_bar/p #initializing, but will be learned
      a_theta <- 1
      b_theta <- 1
      hyp_par_list$num_success <- 0
      hyp_par_list$num_splits <- 0

    }
    if(coef_hyperprior == "univariate_normal_betabinomial_theta_j"){
      hyp_par_list$theta_vec <- rep(p_bar/p, p) # initializing, but will be learned
      a_theta <- 1
      b_theta <- 1

      hyp_par_list$num_success_vec <- rep(0,p)
      hyp_par_list$num_splits <- 0

    }


    if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
      hyp_par_list$theta_vec <- rep(p_bar/p, p) # initializing, but will be learned
      a_theta <- 1
      b_theta <- 1

      hyp_par_list$num_success_vec <- rep(0,p)
      hyp_par_list$num_splits <- 0

      nu_sig <- 3
      tau_sig <- 1/nu_sig

      hyp_par_list$coef_sumsq_vec <- rep(0,p)


    }

    if(coef_norm_hyperprior == "varying"){
      hyp_par_list$coef_sum_vec <- rep(0,p)
      sigma2_beta_bar <- sigma2_beta_vec # difficult to know how to set variance of hyperparameter means


    }

  }


  # if(coef_hyperprior == "univariate_normal"){
  #   hyp_par_list$beta_bar_vec <- rep(0, p) # initialize hyperparameter mean at zero?
  #   # hyp_par_list$sigma_beta_vec <- rep(1, p) # What would be a reasonable initial value?
  #
  #
  #   # ridge-like prior sets prior equal to constant (another hyperpameter) times error variance
  #   ridge_lambda <- 1
  #
  #   sigma2_beta_vec <- rep(ridge_lambda*sigma2, p)
  #
  #   # something like the Zellener G-prior uses the inverse sample covariance matrix of the covariates
  #   # multiplied by sigma
  #   zell_g <- 1
  #   zell_diag <- diag(solve( t(x) %*% x ))
  #   sigma2_beta_vec <- zell_g*zell_diag*sigma2
  #   hyp_par_list$sigma2_beta_vec <- sigma2_beta_vec
  #
  #
  #
  # }









  if(coef_prior == "simplex"){

    # still need to decide whether to constrain beta to be positive
    # or to randomly assign coefficients to be positive or negative
    # or to also put prior on the sign of the coefficients

    hyp_par_list$alpha_simplex <- p/2 # Numbers less than p encourage sparsity
    # See Fisher Simplex regression paper for hyperprior on alpha

    hyp_par_list$xi_vec <- rep(1/p, p)

    if(coef_hyperprior == "simplex_fixed_beta_binomial_theta_j"){
      hyp_par_list$theta_vec <- rep(p_bar/p, p) # initializing, but will be learned
      a_theta <- 1
      b_theta <- 1

      hyp_par_list$num_success_vec <- rep(0,p)
      hyp_par_list$num_splits <- 0

    }

    if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j"){
      hyp_par_list$theta_vec <- rep(1/2, p) # initializing, but will be learned
      a_theta <- 1
      b_theta <- 1

      hyp_par_list$num_success_vec <- rep(0,p)
      hyp_par_list$num_splits <- 0

    }



    if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
      hyp_par_list$theta_vec <- rep(1/2, p) # initializing, but will be learned
      a_theta <- 1
      b_theta <- 1

      hyp_par_list$num_success_vec <- rep(0,p)
      hyp_par_list$num_splits <- 0


      # remove alpha ?
      hyp_par_list$xi_vec <- rep(1/p, p)*hyp_par_list$alpha_simplex

      hyp_par_list$coef_logsum_vec <- rep(0,p)



      # hyperparameters for xi_vec hyperprior
      if(xi_prior == "Boojum"){
        # simplest Boojum prior settings are tau = 0 and v_j = p/alpha. Product of independent exponential distributions
        # Boojum_rate_vec <-  rep(p, p)/hyp_par_list$alpha_simplex
        Boojum_rate <-  p/hyp_par_list$alpha_simplex

      }
      if(xi_prior == "uniform"){
        # uniform from 0 to infinity
        # no parameter actually required for MH step
      }
      if(xi_prior == "gamma"){
        # independent gamma priors
        a_xi_gam <- 1
        b_xi_gam <- p / hyp_par_list$alpha_simplex

      }
      if(xi_prior == "truncnorm"){
        # independent gamma priors
        mu_xi_tnorm <- 1
        sig_xi_tnorm <- 3 # arbitrary

      }




    }


    if(coef_hyperprior == "simplex_fixed_Dir_trinomial_theta_j"){
      theta_1_bar <- rep(p_bar/(p), p) # initializing, but will be learned
      theta_0_bar <- rep(1 - 2*(p_bar/p), p) # initializing, but will be learned
      theta_min1_bar <- rep(p_bar/(p), p) # initializing, but will be learned

      hyp_par_list$theta_1_vec <- rep(p_bar/(p), p) # initializing, but will be learned
      hyp_par_list$theta_0_vec <- rep(1 - 2*(p_bar/p), p) # initializing, but will be learned
      hyp_par_list$theta_min1_vec <- rep(p_bar/(p), p) # initializing, but will be learned

      a_theta <- 1
      b_theta <- 1

      hyp_par_list$num_1_vec <- rep(0,p)
      hyp_par_list$num_0_vec <- rep(0,p)
      hyp_par_list$num_min1_vec <- rep(0,p)
      hyp_par_list$num_splits <- 0

    }

  }



  hyp_par_list$split_mix <- split_mix
  hyp_par_list$num_clust <- num_clust

  if(split_mix){
    # splits are drawn from a mixture distribution

    if(coef_hyperprior %in% c("univariate_normal_betabinomial_theta_j","univariate_normal_betabinomial_theta_j_sigma_j") ){

      hyp_par_list$theta_list <- list()
      hyp_par_list$num_success_list <- list()
      hyp_par_list$num_splits_list <- list()
      hyp_par_list$coef_sumsq_list <- list()
      hyp_par_list$coef_sum_list <- list()
      hyp_par_list$sigma2_beta_list <- list()

      hyp_par_list$beta_bar_list <- list()

      for(c_ind in 1:num_clust){
        hyp_par_list$theta_list[[c_ind]] <- rep(p_bar/p, p) # initializing, but will be learned
        hyp_par_list$num_success_list[[c_ind]] <- rep(0,p)
        hyp_par_list$num_splits_list[[c_ind]] <- 0

        hyp_par_list$coef_sumsq_list[[c_ind]] <- rep(0,p)

        hyp_par_list$coef_sum_list[[c_ind]] <- rep(0,p)
        hyp_par_list$sigma2_beta_list[[c_ind]] <- sigma2_beta_vec # difficult to know how to set variance of hyperparameter means

        hyp_par_list$beta_bar_list[[c_ind]] <- rep(0, p)

      }



      hyp_par_list$clust_counts <- rep(0, num_clust)
      # initialize cluster probabilties
      hyp_par_list$clust_probs <- rep(1/num_clust, num_clust)

      # dirichlet prior on cluster probabilities
      hyp_par_list$clust_dir_params <- rep(1, num_clust)

      hyp_par_list$clust_ind_list <- list()

      for(tree_ind in 1:ntrees){
        hyp_par_list$clust_ind_list[[tree_ind]] <- rep(NA,1)
      }

    }


  }












  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60,
                             title = 'Running rBART...')

  ##### Start the MCMC iterations loop ##########
  for (i in 1:TotIter) {

    utils::setTxtProgressBar(pb, i)

    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr = (i - nburn)/nthin
      tree_store[[curr]] = curr_trees
      sigma2_store[curr] = sigma2
      y_hat_store[curr,] = y_hat
      # var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
    }

      ###### Start looping through trees #########
      for (j in 1:ntrees) {

        current_partial_residuals = y_scale - y_hat + tree_fits_store[,j]

        # Propose a new tree via grow/change/prune/swap
        type = sample_move(curr_trees[[j]], i, 100 #nburn
                           )

        # Generate a new tree based on the current
        new_trees[[j]] = update_tree(y = y_scale,
                                     X = x,
                                     type = type,
                                     curr_tree = curr_trees[[j]],
                                     node_min_size = node_min_size,
                                     s = s,
                                     coef_prior,        # Prior distribution for splitting coefficients
                                     coef_hyperprior,   # hyperprior for splitting coefficients
                                     hyp_par_list,      # list of hyperparameters for splitting coefficients,
                                     threshold_prior,    # threshold prior
                                     coef_norm_hyperprior, # if normal hyperprior, fixed (i.e just vairable selection), or actually adapt mean of normal?
                                     norm_unit_sphere
                                     )

        # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
        l_old = tree_full_conditional(curr_trees[[j]],
                                      current_partial_residuals,
                                      sigma2,
                                      sigma2_mu) +
          get_tree_prior(curr_trees[[j]], alpha, beta)

        # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
        l_new = tree_full_conditional(new_trees[[j]],
                                      current_partial_residuals,
                                      sigma2,
                                      sigma2_mu) +
          get_tree_prior(new_trees[[j]], alpha, beta)



        a = alpha_mh(l_new,l_old, curr_trees[[j]],new_trees[[j]], type)

        # Exponentiate the results above
        # a = exp(l_new - l_old)

        if(a > runif(1)) {
          curr_trees[[j]] = new_trees[[j]]


        ### mix update counts and sums  #############
          if(split_mix){
            if(curr_trees[[j]]$var_update == 1){


              if (type =='change'){
                temp_nodeind <- new_trees[[j]]$node_updated
                c_ind_old <- hyp_par_list$clust_ind_list[[j]][temp_nodeind]
                c_ind_new <- new_trees[[j]]$c_ind_new
                hyp_par_list$clust_ind_list[[j]][temp_nodeind] <-  c_ind_new

                hyp_par_list$clust_counts[c_ind_old] <- hyp_par_list$clust_counts[c_ind_old] - 1
                hyp_par_list$clust_counts[c_ind_new] <- hyp_par_list$clust_counts[c_ind_new] + 1

              }
              if (type=='grow'){
                temp_nodeind <- new_trees[[j]]$node_updated
                # c_ind_old <- hyp_par_list$clust_ind_list[[j]][temp_nodeind]
                c_ind_new <- new_trees[[j]]$c_ind_new
                hyp_par_list$clust_ind_list[[j]][temp_nodeind] <-  c_ind_new

                hyp_par_list$clust_ind_list[[j]] <- c(hyp_par_list$clust_ind_list[[j]],
                                                      NA,NA)

                hyp_par_list$clust_counts[c_ind_new] <- hyp_par_list$clust_counts[c_ind_new] + 1


              }
              if (type=='prune'){
                par_ind <- new_trees[[j]]$parent_pruned

                c_ind_par <- hyp_par_list$clust_ind_list[[j]][par_ind]
                hyp_par_list$clust_ind_list[[j]][par_ind] <- NA
                #remove pruned nodes
                hyp_par_list$clust_ind_list[[j]] <- hyp_par_list$clust_ind_list[[j]][ - new_trees[[j]]$nodes_pruned ]

                hyp_par_list$clust_counts[c_ind_par] <- hyp_par_list$clust_counts[c_ind_par] - 1
              }
            }



            # ADD IF STATEMENT FOR RELEVANT HYPERPRIOR
            if(coef_hyperprior == "univariate_normal_betabinomial_onetheta"){
              if(curr_trees[[j]]$var_update == 1){
                if (type =='change'){
                  # hyp_par_list$num_success_list[[c_ind]] <- hyp_par_list$num_success_list[[c_ind]] + sum(curr_trees[[j]]$count_for_update)

                  # must add to count for new cluster
                  hyp_par_list$num_success_list[[c_ind_new]] <- hyp_par_list$num_success_list[[c_ind_new]] + sum(curr_trees[[j]]$count_for_update_new)

                  # and take away from count for old cluster
                  hyp_par_list$num_success_list[[c_ind_old]] <- hyp_par_list$num_success_list[[c_ind_old]] - sum(curr_trees[[j]]$count_for_update_old)


                  if(coef_norm_hyperprior == "varying"){
                    hyp_par_list$coef_sum_list[[c_ind_new]] <- hyp_par_list$coef_sum_list[[c_ind_new]] + curr_trees[[j]]$coef_for_sum_vec_new
                    hyp_par_list$coef_sum_list[[c_ind_old]] <- hyp_par_list$coef_sum_list[[c_ind_old]] - curr_trees[[j]]$coef_for_sum_vec_old
                  }
                }
                if (type=='grow'){
                  hyp_par_list$num_success_list[[c_ind_new]] <- hyp_par_list$num_success_list[[c_ind_new]] + sum(curr_trees[[j]]$count_for_update)

                  hyp_par_list$num_splits_list[[c_ind_new]] <- hyp_par_list$num_splits_list[[c_ind_new]] + 1

                  if(coef_norm_hyperprior == "varying"){
                    hyp_par_list$coef_sum_list[[c_ind_new]] <- hyp_par_list$coef_sum_list[[c_ind_new]] + curr_trees[[j]]$coef_for_sum_vec
                  }
                }
                if (type=='prune'){
                  # hyp_par_list$num_success <- hyp_par_list$num_success - sum(curr_trees[[j]]$count_for_update)
                  # hyp_par_list$num_splits <- hyp_par_list$num_splits - 1

                  hyp_par_list$num_success_list[[c_ind_par]] <- hyp_par_list$num_success_list[[c_ind_par]] -
                    sum(curr_trees[[j]]$count_for_update)

                  hyp_par_list$num_splits_list[[c_ind_par]] <- hyp_par_list$num_splits_list[[c_ind_par]] - 1

                  if(coef_norm_hyperprior == "varying"){
                    hyp_par_list$coef_sum_list[[c_ind_par]] <- hyp_par_list$coef_sum_list[[c_ind_par]] - curr_trees[[j]]$coef_for_sum_vec
                  }
                }
              }
            }

            if( coef_hyperprior %in% c( "univariate_normal_betabinomial_theta_j" ,
                                        "univariate_normal_betabinomial_theta_j_sigma_j",
                                        "simplex_fixed_beta_binomial_theta_j",
                                        "simplex_fixed_Dir_binomial_plusminus_theta_j",
                                        "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j") ){

              if(curr_trees[[j]]$var_update == 1){
                if (type =='change'){
                  # hyp_par_list$num_success_vec <- hyp_par_list$num_success_vec + curr_trees[[j]]$count_for_update
                  #
                  #must update vectors for new and old clusters
                  hyp_par_list$num_success_list[[c_ind_new]] <- hyp_par_list$num_success_list[[c_ind_new]] + curr_trees[[j]]$count_for_update_new

                  hyp_par_list$num_success_list[[c_ind_old]] <- hyp_par_list$num_success_list[[c_ind_old]] - curr_trees[[j]]$count_for_update_old


                  if(coef_norm_hyperprior == "varying"){
                    # hyp_par_list$coef_sum_vec <- hyp_par_list$coef_sum_vec + curr_trees[[j]]$coef_for_sum_vec
                    hyp_par_list$coef_sum_list[[c_ind_new]] <- hyp_par_list$coef_sum_list[[c_ind_new]] + curr_trees[[j]]$coef_for_sum_vec_new
                    hyp_par_list$coef_sum_list[[c_ind_old]] <- hyp_par_list$coef_sum_list[[c_ind_old]] - curr_trees[[j]]$coef_for_sum_vec_old
                  }
                }
                if (type=='grow'){
                  # hyp_par_list$num_success_vec <- hyp_par_list$num_success_vec + curr_trees[[j]]$count_for_update
                  # hyp_par_list$num_splits <- hyp_par_list$num_splits + 1

                  hyp_par_list$num_success_list[[c_ind_new]] <- hyp_par_list$num_success_list[[c_ind_new]] + curr_trees[[j]]$count_for_update
                  hyp_par_list$num_splits_list[[c_ind_new]] <- hyp_par_list$num_splits_list[[c_ind_new]] + 1

                  if(coef_norm_hyperprior == "varying"){
                    hyp_par_list$coef_sum_list[[c_ind_new]] <- hyp_par_list$coef_sum_list[[c_ind_new]] + curr_trees[[j]]$coef_for_sum_vec

                  }
                }
                if (type=='prune'){
                  # hyp_par_list$num_success_vec <- hyp_par_list$num_success_vec - curr_trees[[j]]$count_for_update
                  # hyp_par_list$num_splits <- hyp_par_list$num_splits - 1

                  hyp_par_list$num_success_list[[c_ind_par]] <- hyp_par_list$num_success_list[[c_ind_par]] - curr_trees[[j]]$count_for_update
                  hyp_par_list$num_splits_list[[c_ind_par]] <- hyp_par_list$num_splits_list[[c_ind_par]] - 1

                  if(coef_norm_hyperprior == "varying"){
                    hyp_par_list$coef_sum_list[[c_ind_par]] <- hyp_par_list$coef_sum_list[[c_ind_par]] - curr_trees[[j]]$coef_for_sum_vec
                  }
                }
              }
            }


            if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
              if(curr_trees[[j]]$var_update == 1){
                if (type =='change'){
                  # hyp_par_list$coef_sumsq_vec <- hyp_par_list$coef_sumsq_vec + curr_trees[[j]]$coef_for_sumsq_vec
                  hyp_par_list$coef_sumsq_list[[c_ind_new]] <- hyp_par_list$coef_sumsq_list[[c_ind_new]] + curr_trees[[j]]$coef_for_sumsq_vec_new
                  hyp_par_list$coef_sumsq_list[[c_ind_old]] <- hyp_par_list$coef_sumsq_list[[c_ind_old]] - curr_trees[[j]]$coef_for_sumsq_vec_old
                }
                if (type=='grow'){
                  # hyp_par_list$coef_sumsq_vec <- hyp_par_list$coef_sumsq_vec + curr_trees[[j]]$coef_for_sumsq_vec
                  hyp_par_list$coef_sumsq_list[[c_ind_new]] <- hyp_par_list$coef_sumsq_list[[c_ind_new]] + curr_trees[[j]]$coef_for_sumsq_vec
                }
                if (type=='prune'){
                  # hyp_par_list$coef_sumsq_vec <- hyp_par_list$coef_sumsq_vec - curr_trees[[j]]$coef_for_sumsq_vec
                  hyp_par_list$coef_sumsq_list[[c_ind_par]] <- hyp_par_list$coef_sumsq_list[[c_ind_par]] - curr_trees[[j]]$coef_for_sumsq_vec

                }
              }
            }

          }else{ # end split mix TRUE

            ### NOT mix update counts and sums  #############


            # ADD IF STATEMENT FOR RELEVANT HYPERPRIOR
            if(coef_hyperprior == "univariate_normal_betabinomial_onetheta"){
              if(curr_trees[[j]]$var_update == 1){
                if (type =='change'){
                  hyp_par_list$num_success <- hyp_par_list$num_success + sum(curr_trees[[j]]$count_for_update)

                  # print(" type=='change' ")
                  # print("hyp_par_list$num_success = ")
                  # print(hyp_par_list$num_success)
                  # print("hyp_par_list$num_splits = ")
                  # print(hyp_par_list$num_splits)
                  # print("curr_trees[[j]]$count_for_update = ")
                  # print(curr_trees[[j]]$count_for_update)
                  if(coef_norm_hyperprior == "varying"){
                    hyp_par_list$coef_sum_vec <- hyp_par_list$coef_sum_vec + curr_trees[[j]]$coef_for_sum_vec
                  }
                }
                if (type=='grow'){
                  hyp_par_list$num_success <- hyp_par_list$num_success + sum(curr_trees[[j]]$count_for_update)
                  hyp_par_list$num_splits <- hyp_par_list$num_splits + 1

                  # print(" type=='grow' ")
                  # print("hyp_par_list$num_success = ")
                  # print(hyp_par_list$num_success)
                  # print("hyp_par_list$num_splits = ")
                  # print(hyp_par_list$num_splits)
                  # print("curr_trees[[j]]$count_for_update = ")
                  # print(curr_trees[[j]]$count_for_update)
                  if(coef_norm_hyperprior == "varying"){
                    hyp_par_list$coef_sum_vec <- hyp_par_list$coef_sum_vec + curr_trees[[j]]$coef_for_sum_vec
                  }
                }
                if (type=='prune'){
                  hyp_par_list$num_success <- hyp_par_list$num_success - sum(curr_trees[[j]]$count_for_update)
                  hyp_par_list$num_splits <- hyp_par_list$num_splits - 1

                  # print(" type=='prune' ")
                  # print("hyp_par_list$num_success = ")
                  # print(hyp_par_list$num_success)
                  # print("hyp_par_list$num_splits = ")
                  # print(hyp_par_list$num_splits)
                  # print("curr_trees[[j]]$count_for_update = ")
                  # print(curr_trees[[j]]$count_for_update)
                  if(coef_norm_hyperprior == "varying"){
                    hyp_par_list$coef_sum_vec <- hyp_par_list$coef_sum_vec - curr_trees[[j]]$coef_for_sum_vec
                  }
                }
              }
            }

            if( coef_hyperprior %in% c( "univariate_normal_betabinomial_theta_j" ,
                                        "univariate_normal_betabinomial_theta_j_sigma_j",
                                        "simplex_fixed_beta_binomial_theta_j",
                                        "simplex_fixed_Dir_binomial_plusminus_theta_j",
                                        "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j") ){

              if(curr_trees[[j]]$var_update == 1){
                if (type =='change'){
                  hyp_par_list$num_success_vec <- hyp_par_list$num_success_vec + curr_trees[[j]]$count_for_update
                  if(coef_norm_hyperprior == "varying"){
                    hyp_par_list$coef_sum_vec <- hyp_par_list$coef_sum_vec + curr_trees[[j]]$coef_for_sum_vec
                  }
                }
                if (type=='grow'){
                  hyp_par_list$num_success_vec <- hyp_par_list$num_success_vec + curr_trees[[j]]$count_for_update
                  hyp_par_list$num_splits <- hyp_par_list$num_splits + 1
                  if(coef_norm_hyperprior == "varying"){
                    hyp_par_list$coef_sum_vec <- hyp_par_list$coef_sum_vec + curr_trees[[j]]$coef_for_sum_vec
                  }
                }
                if (type=='prune'){
                  hyp_par_list$num_success_vec <- hyp_par_list$num_success_vec - curr_trees[[j]]$count_for_update
                  hyp_par_list$num_splits <- hyp_par_list$num_splits - 1
                  if(coef_norm_hyperprior == "varying"){
                    hyp_par_list$coef_sum_vec <- hyp_par_list$coef_sum_vec - curr_trees[[j]]$coef_for_sum_vec
                  }
                }
              }
            }

            if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
              if(curr_trees[[j]]$var_update == 1){
                if (type =='change'){
                  hyp_par_list$coef_logsum_vec <- hyp_par_list$coef_logsum_vec + curr_trees[[j]]$coef_for_logsum_vec
                }
                if (type=='grow'){
                  hyp_par_list$coef_logsum_vec <- hyp_par_list$coef_logsum_vec + curr_trees[[j]]$coef_for_logsum_vec
                }
                if (type=='prune'){
                  hyp_par_list$coef_logsum_vec <- hyp_par_list$coef_logsum_vec - curr_trees[[j]]$coef_for_logsum_vec
                }
              }
            }

            if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
              if(curr_trees[[j]]$var_update == 1){
                if (type =='change'){
                  hyp_par_list$coef_sumsq_vec <- hyp_par_list$coef_sumsq_vec + curr_trees[[j]]$coef_for_sumsq_vec
                }
                if (type=='grow'){
                  hyp_par_list$coef_sumsq_vec <- hyp_par_list$coef_sumsq_vec + curr_trees[[j]]$coef_for_sumsq_vec
                }
                if (type=='prune'){
                  hyp_par_list$coef_sumsq_vec <- hyp_par_list$coef_sumsq_vec - curr_trees[[j]]$coef_for_sumsq_vec
                }
              }
            }


            if(coef_hyperprior == "simplex_fixed_Dir_trinomial_theta_j"){
              if(curr_trees[[j]]$var_update == 1){
                if (type =='change'){
                  hyp_par_list$num_min1_vec <- hyp_par_list$num_min1_vec + curr_trees[[j]]$count_min1_for_update
                  hyp_par_list$num_0_vec <- hyp_par_list$num_0_vec + curr_trees[[j]]$count_0_for_update
                  hyp_par_list$num_1_vec <- hyp_par_list$num_1_vec + curr_trees[[j]]$count_1_for_update
                }
                if (type=='grow'){
                  hyp_par_list$num_min1_vec <- hyp_par_list$num_min1_vec + curr_trees[[j]]$count_min1_for_update
                  hyp_par_list$num_0_vec <- hyp_par_list$num_0_vec + curr_trees[[j]]$count_0_for_update
                  hyp_par_list$num_1_vec <- hyp_par_list$num_1_vec + curr_trees[[j]]$count_1_for_update
                  hyp_par_list$num_splits <- hyp_par_list$num_splits + 1

                }
                if (type=='prune'){
                  hyp_par_list$num_min1_vec <- hyp_par_list$num_min1_vec - curr_trees[[j]]$count_min1_for_update
                  hyp_par_list$num_0_vec <- hyp_par_list$num_0_vec - curr_trees[[j]]$count_0_for_update
                  hyp_par_list$num_1_vec <- hyp_par_list$num_1_vec - curr_trees[[j]]$count_1_for_update
                  hyp_par_list$num_splits <- hyp_par_list$num_splits - 1

                }
              }
            }

          } # end else for split mix FALSE
        } # end if accept MH step code




        # Update mu whether tree accepted or not
        curr_trees[[j]] = simulate_mu(curr_trees[[j]],
                                      current_partial_residuals,
                                      sigma2,
                                      sigma2_mu)
      # Updating BART predictions
      current_fit = get_predictions(curr_trees[j], x, single_tree = TRUE)
      y_hat = y_hat - tree_fits_store[,j] # subtract the old fit
      y_hat = y_hat + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

      } # End loop through trees


    ####### draw sigma ############

    # y_hat = get_predictions(curr_trees, x, single_tree = ntrees == 1)
    sum_of_squares = sum((y_scale - y_hat)^2)

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)

    if(is.na(sigma2)){
      stop("sigma2 = NA")
    }


    if(coef_prior == "univariate_normal"){
      # sigma2_beta_scalar <- 1 # what would be a reasonable initial value?
      # maybe it would be better to have a vector of coefficient values? Scaled by a linear regression or ridge?
      # Scale covariates first?

      # ridge-like prior sets prior equal to constant (another hyperpameter) times error variance
      # ridge_lambda <- 1

      if(norm_sigma_init == "ridge"){
        sigma2_beta_vec <- rep(ridge_lambda*sigma2, p)
      }
      # something like the Zellener G-prior uses the inverse sample covariance matrix of the covariates
      # multiplied by sigma
      # zell_g <- 1
      # zell_diag <- diag(solve(t(x)%*%X))

      if(norm_sigma_init == "Zellner"){
        sigma2_beta_vec <- zell_g*zell_diag*sigma2
      }

      hyp_par_list$sigma2_beta_vec <- sigma2_beta_vec

    }

    # # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    # if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
    #   s = update_s(var_count, p, 1)
    # }

    ##### mix draw hyperparameters ###########

    if(split_mix){

      # draw new cluster probabiltiies

      hyp_par_list$clust_probs <- t(rdirichlet(n = 1,
                                               alpha = hyp_par_list$clust_dir_params  +
                                                 hyp_par_list$clust_counts))



      for(c_ind in 1:num_clust){

        if( coef_hyperprior %in% c( "univariate_normal_betabinomial_theta_j" ,
                                    "univariate_normal_betabinomial_theta_j_sigma_j" ,
                                    "simplex_fixed_beta_binomial_theta_j",
                                    "simplex_fixed_Dir_binomial_plusminus_theta_j",
                                    "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j") ){
          for(j in 1:p){
            hyp_par_list$theta_list[[c_ind]][j] <- rbeta(n = 1,
                                               shape1 = a_theta + hyp_par_list$num_success_list[[c_ind]][j] ,
                                               shape2 = b_theta + hyp_par_list$clust_counts[c_ind] - hyp_par_list$num_success_list[[c_ind]][j] )

            if(is.na(hyp_par_list$theta_list[[c_ind]][j])){
              print("a_theta = ")
              print(a_theta)

              print("hyp_par_list$num_success_list[[c_ind]][j] = ")
              print(hyp_par_list$num_success_list[[c_ind]][j])

              print("hyp_par_list$num_splits_list[[c_ind]] = ")
              print(hyp_par_list$num_splits_list[[c_ind]])
              print("hyp_par_list$clust_counts[c_ind] = ")
              print(hyp_par_list$clust_counts[c_ind])

              print("hyp_par_list$num_success_list[[c_ind]][j] = ")
              print( hyp_par_list$num_success_list[[c_ind]][j] )

              print("j = ")
              print(j)

              stop(" NA hyp_par_list$theta_list[[c_ind]][j] ")
            }

          }
        }


        if(coef_norm_hyperprior == "varying"){
          for(j in 1:p){
            tempvar <- 1/( (1/sigma2_beta_bar[j]) + (hyp_par_list$num_success_list[[c_ind]][j]/hyp_par_list$sigma2_beta_list[[c_ind]][j])   )
            # print("sigma2_beta_bar[j] = ")
            # print(sigma2_beta_bar[j])
            hyp_par_list$beta_bar_list[[c_ind]][j] <- rnorm(1,
                                                  mean = tempvar* hyp_par_list$coef_sum_list[[c_ind]][j] /hyp_par_list$sigma2_beta_list[[c_ind]][j],
                                                  sd = sqrt(tempvar)
            )

          }
        }

        if(coef_hyperprior  == "univariate_normal_betabinomial_theta_j_sigma_j"){

          for(j in 1:p){
            hyp_par_list$sigma2_beta_list[[c_ind]][j] <- 1/rgamma(n = 1,
                                                                  shape = (nu_sig +  hyp_par_list$num_success_list[[c_ind]][j])/2,
                                                                  rate = (nu_sig*tau_sig +
                                                                            hyp_par_list$coef_sumsq_list[[c_ind]][j] -
                                                                            2*hyp_par_list$beta_bar_list[[c_ind]][j]*hyp_par_list$coef_sum_list[[c_ind]][j] +
                                                                            hyp_par_list$num_success_list[[c_ind]][j]*(hyp_par_list$beta_bar_list[[c_ind]][j]^2))/2)

            if(is.na(hyp_par_list$sigma2_beta_list[[c_ind]][j])){
              print("nu_sig = ")
              print(nu_sig)

              print("hyp_par_list$coef_sum_list[[c_ind]][j] = ")
              print(hyp_par_list$coef_sum_list[[c_ind]][j])

              print("tau_sig = ")
              print(tau_sig)


              print("hyp_par_list$coef_sumsq_list[[c_ind]][j] = ")
              print(hyp_par_list$coef_sumsq_list[[c_ind]][j])

              print("hyp_par_list$num_success_list[[c_ind]][j] = ")
              print(hyp_par_list$num_success_list[[c_ind]][j])

              print("hyp_par_list$beta_bar_list[[c_ind]][j]^2 = ")
              print(hyp_par_list$beta_bar_list[[c_ind]][j]^2)


              print("hyp_par_list$coef_sumsq_list[[c_ind]][j] -
              2*hyp_par_list$beta_bar_list[[c_ind]][j]*hyp_par_list$coef_sum_list[[c_ind]][j] +
                    hyp_par_list$num_success_list[[c_ind]][j]*(hyp_par_list$beta_bar_list[[c_ind]][j]^2) = ")

              print(hyp_par_list$coef_sumsq_list[[c_ind]][j] -
                      2*hyp_par_list$beta_bar_list[[c_ind]][j]*hyp_par_list$coef_sum_list[[c_ind]][j] +
                      hyp_par_list$num_success_list[[c_ind]][j]*(hyp_par_list$beta_bar_list[[c_ind]][j]^2))

              print("j = ")
              print(j)


              stop("NA hyp_par_list$sigma2_beta_list[[c_ind]][j]")
            }
          }
        }



      }


    }else{ ################  not mix draw hyperparameters ########################

      if(coef_hyperprior == "univariate_normal_betabinomial_onetheta"){
        hyp_par_list$theta <- rbeta(n = 1,
                                    shape1 = a_theta + hyp_par_list$num_success ,
                                    shape2 = b_theta + p* hyp_par_list$num_splits - hyp_par_list$num_success )

        if(is.na(hyp_par_list$theta)){
          print("a_theta = ")
          print(a_theta)

          print("b_theta = ")
          print(b_theta)

          print("hyp_par_list$num_success = ")
          print(hyp_par_list$num_success)

          print("hyp_par_list$num_splits = ")
          print(hyp_par_list$num_splits)

        }

      }
      if( coef_hyperprior %in% c( "univariate_normal_betabinomial_theta_j" ,
                                  "univariate_normal_betabinomial_theta_j_sigma_j" ,
                                  "simplex_fixed_beta_binomial_theta_j",
                                  "simplex_fixed_Dir_binomial_plusminus_theta_j",
                                  "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j") ){
        for(j in 1:p){
          hyp_par_list$theta_vec[j] <- rbeta(n = 1,
                                       shape1 = a_theta + hyp_par_list$num_success_vec[j] ,
                                       shape2 = b_theta + hyp_par_list$num_splits - hyp_par_list$num_success[j] )
        }
      }


      if(coef_norm_hyperprior == "varying"){
        for(j in 1:p){
          tempvar <- 1/( (1/sigma2_beta_bar[j]) + (hyp_par_list$num_success[j]/sigma2_beta_vec[j])   )
          # print("sigma2_beta_bar[j] = ")
          # print(sigma2_beta_bar[j])
          hyp_par_list$beta_bar_vec[j] <- rnorm(1,
                                                mean = tempvar* hyp_par_list$coef_sum_vec[j] /sigma2_beta_vec[j],
                                                sd = sqrt(tempvar)
                                                )

        }
      }

      if(coef_hyperprior  == "univariate_normal_betabinomial_theta_j_sigma_j"){

        for(j in 1:p){
          sigma2_beta_vec[j] <- 1/rgamma(n = 1,
                                         shape = (nu_sig +  hyp_par_list$num_success_vec[j])/2,
                                         rate = (nu_sig*tau_sig +
                                                   hyp_par_list$coef_sumsq_vec[j] -
                                                   2*hyp_par_list$beta_bar_vec[j]*hyp_par_list$coef_sum_vec[j] +
                                                   hyp_par_list$num_success_vec[j]*(hyp_par_list$beta_bar_vec[j]^2))/2)

          if(is.na(sigma2_beta_vec[j])){
            print("nu_sig = ")
            print(nu_sig)

            print("hyp_par_list$coef_sum_vec[j] = ")
            print(hyp_par_list$coef_sum_vec[j])

            print("tau_sig = ")
            print(tau_sig)


            print("hyp_par_list$coef_sumsq_vec[j] = ")
            print(hyp_par_list$coef_sumsq_vec[j])

            print("j = ")
            print(j)


            stop("NA sigma2_beta_vec[j]")
          }
        }
      }


      if(coef_hyperprior == "simplex_fixed_Dir_trinomial_theta_j"){
        for(j in 1:p){
          temp_theta_vec <- rdirichlet(n = 1,
                                       alpha = c(theta_min1_bar +  hyp_par_list$num_min1_vec,
                                                 theta_0_bar +  hyp_par_list$num_0_vec,
                                                 theta_1_bar +  hyp_par_list$num_1_vec))

          hyp_par_list$theta_min1_vec[j] <- temp_theta_vec[1]
          hyp_par_list$theta_0_vec[j] <- temp_theta_vec[2]
          hyp_par_list$theta_1_vec[j] <- temp_theta_vec[3]
        }
      }


      if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){


        orig_xi <- hyp_par_list$xi_vec

        prop_xi <- rep(NA,p)

        for(j in 1:p){
          if(orig_xi[j] - MH_propstep_xi <= 0){
            prop_xi[j] <- runif(n = 1,
                             min = 0,
                             max = orig_xi[j] + MH_propstep_xi )
          }else{
            prop_xi[j] <- runif(n = 1,
                                min = orig_xi[j] - MH_propstep_xi,
                                max = orig_xi[j] + MH_propstep_xi )
          }
        }

        templogbeta_orig <- sum(lgamma(orig_xi)) - lgamma(sum(orig_xi))
        templogbeta_prop <- sum(lgamma(prop_xi)) - lgamma(sum(prop_xi))


        log_MH_lik_orig <- -1*hyp_par_list$num_splits*templogbeta_orig  +  sum( (orig_xi - 1)*hyp_par_list$coef_logsum_vec )
        log_MH_lik_prop <- -1*hyp_par_list$num_splits*templogbeta_prop  +  sum( (prop_xi - 1)*hyp_par_list$coef_logsum_vec )


        # hyperparameters for xi_vec hyperprior
        if(xi_prior == "Boojum"){
          # simplest Boojum prior settings are tau = 0 and v_j = p/alpha. Product of independent exponential distributions
          # Boojum_rate_vec <-  rep(p, p)/hyp_par_list$alpha_simplex

          log_prior_ratio <- sum(dexp(prop_xi,  rate = Boojum_rate, log = TRUE ) -
                                   dexp(orig_xi,  rate = Boojum_rate, log = TRUE ))

        }
        if(xi_prior == "uniform"){
          # uniform from 0 to infinity
          # no parameter actually required for MH step
          log_prior_ratio <- 0
        }
        if(xi_prior == "gamma"){
          # independent gamma priors
          log_prior_ratio <- sum(dgamma(prop_xi, shape = a_xi_gam, rate = b_xi_gam, log = TRUE ) -
            dgamma(orig_xi, shape = a_xi_gam, rate = b_xi_gam, log = TRUE ))
        }
        if(xi_prior == "truncnorm"){
          # independent gamma priors
          log_prior_ratio <- sum(dnorm(prop_xi, mean = mu_xi_tnorm, sd = sig_xi_tnorm, log = TRUE ) -
            dnorm(orig_xi, mean = mu_xi_tnorm, sd = sig_xi_tnorm, log = TRUE ))
        }

        MH_accept_prob <- exp(log_prior_ratio + log_MH_lik_prop - log_MH_lik_orig)

        if(is.na(MH_accept_prob)){
          print("log_prior_ratio = ")
          print(log_prior_ratio)

          print("log_MH_lik_prop = ")
          print(log_MH_lik_prop)

          print("log_MH_lik_orig = ")
          print(log_MH_lik_orig)

          print("MH_accept_prob = ")
          print(MH_accept_prob)

          stop("NA MH_accept_prob")

        }


        # print("MH_accept_prob = ")
        # print(MH_accept_prob)


        if(MH_accept_prob > runif(1)){
          hyp_par_list$xi_vec <- prop_xi
        }

      } # end xi MH step

    } # end else statement (split mix FALSE)


  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  ###### return results #############
  return(list(trees = tree_store,
              sigma2 = sigma2_store*y_sd^2,
              y_hat = y_hat_store*y_sd + y_mean,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd,
              # var_count_store = var_count_store,
              s = s_prob_store,
              ecdfs = ecdfs,
              hyp_par_list = hyp_par_list
              ))

} # End main function

