

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
ObliqueBART <- function(   x,
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
                   nburn = 1000,
                   npost = 1000,
                   nthin = 1,
                   penalise_num_cov = TRUE,
                   lambda_cov = 0.4,
                   nu_cov = 2,
                   coef_prior = "univariate_normal",        # Prior distribution for splitting coefficients
                   coef_hyperprior = "none",
                   coef_norm_hyperprior = "fixed", # fixed at zero, or varying
                   threshold_prior = "discrete_uniform",
                   p_bar = ceiling(ncol(x)/3)
                   ) {


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
                              "multivariate_normal"#,
                         # "simplex"
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
                                  "univariate_normal_betabinomial_theta_j" ) )){
        stop("Currently only the univariate normal hyperprior is supported for a univariate normal coefficient prior")
      }
    }
    if(coef_prior == "multivariate_normal"){
      if(coef_hyperprior != "multivariate_normal"){
        stop("Currently only the multivariate normal hyperprior is supported for a multivariate normal coefficient prior")
      }
    }
    if(coef_prior == "simplex"){
      stop("Currently no hyperpriors supported for simplex coefficient prior")
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

    # ridge-like prior sets prior equal to constant (another hyperpameter) times error variance
    ridge_lambda <- 1

    sigma2_beta_vec <- rep(ridge_lambda*sigma2, p)

    # something like the Zellener G-prior uses the inverse sample covariance matrix of the covariates
    # multiplied by sigma
    zell_g <- 1
    zell_diag <- diag(solve( t(x) %*% x ))
    sigma2_beta_vec <- zell_g*zell_diag*sigma2

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
        type = sample_move(curr_trees[[j]], i, nburn)

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
                                     threshold_prior    # threshold prior
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
              }
            }
          }
          if(coef_hyperprior == "univariate_normal_betabinomial_theta_j"){
            if(curr_trees[[j]]$var_update == 1){
              if (type =='change'){
                hyp_par_list$num_success_vec <- hyp_par_list$num_success_vec + curr_trees[[j]]$count_for_update
              }
              if (type=='grow'){
                hyp_par_list$num_success_vec <- hyp_par_list$num_success_vec + curr_trees[[j]]$count_for_update
                hyp_par_list$num_splits <- hyp_par_list$num_splits + 1
              }
              if (type=='prune'){
                hyp_par_list$num_success_vec <- hyp_par_list$num_success_vec - curr_trees[[j]]$count_for_update
                hyp_par_list$num_splits <- hyp_par_list$num_splits - 1
              }
            }
          }
        }




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

      sigma2_beta_vec <- rep(ridge_lambda*sigma2, p)

      # something like the Zellener G-prior uses the inverse sample covariance matrix of the covariates
      # multiplied by sigma
      # zell_g <- 1
      # zell_diag <- diag(solve(t(x)%*%X))
      sigma2_beta_vec <- zell_g*zell_diag*sigma2

      hyp_par_list$sigma2_beta_vec <- sigma2_beta_vec

    }

    # # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    # if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
    #   s = update_s(var_count, p, 1)
    # }

    ##### draw hyperparameters ###########

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
    if(coef_hyperprior == "univariate_normal_betabinomial_theta_j"){
      for(j in 1:p){
        hyp_par_list$theta_vec[j] <- rbeta(n = 1,
                                     shape1 = a_theta + hyp_par_list$num_success_vec[j] ,
                                     shape2 = b_theta + hyp_par_list$num_splits - hyp_par_list$num_success[j] )
      }
    }



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

