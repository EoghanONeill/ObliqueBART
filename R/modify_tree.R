# -------------------------------------------------------------------------#
# Description: this script contains functions that are used to perform the #
#              growing, pruning, changing, and swapping moves. It also has #
#              a function to initialise the trees to stumps                #
# -------------------------------------------------------------------------#

# 01. create_stump: initialises the trees to a stump
# 02. update_tree: calls the corresponding function associated to the move grow, prune, change, or swap.
# 03. grow_tree: grows a tree
# 04. prune_tree: prunes a tree
# 05. change_tree: changes the splitting rule that defines a pair of terminal nodes
# 06. swap_tree: exchanges the splitting rules that define two pair of terminal nodes

# Function to create stump ------------------------------------------------

create_stump = function(num_trees,
                        y,
                        X) {

  # Each tree is a list of 2 elements
  # The 2 elements are the tree matrix (8 columns), and the node indices
  # The columns of the tree matrix are:
  # Terminal (0 = no, 1 = yes)
  # Child left
  # Child right
  # Node parents
  # Split variable
  # Split value
  # mu
  # Node size

  # Create holder for trees
  all_trees = vector('list', length = num_trees)
  # Loop through trees
  for (j in 1:num_trees) {
    # Set up each tree to have two elements in the list as described above
    all_trees[[j]] = vector('list', length = 2)
    # Give the elements names
    names(all_trees[[j]]) = c('tree_matrix',
                              'node_indices')
    # Create the two elements: first is a matrix
    # all_trees[[j]][[1]] = matrix(NA, ncol = 8, nrow = 1)
    all_trees[[j]][[1]] = matrix(NA, ncol = 7 + ncol(x), nrow = 1)

    # Second is the assignment to node indices
    all_trees[[j]][[2]] = rep(1, length(y))

    # Create column names
    colnames(all_trees[[j]][[1]]) = c('terminal',
                                      'child_left',
                                      'child_right',
                                      'parent',
                                      # 'split_variable',
                                      'split_value',
                                      'mu',
                                      'node_size',
                                      paste0("coef_", 1:ncol(x)))

    # Set values for stump
    all_trees[[j]][[1]][1,] = c(1, NA, NA, NA, #NA,
                                NA, 0 , length(y),
                                rep(NA, ncol(x)))

  } # End of loop through trees

  return(all_trees)

} # End of function

# Function to update trees ------------------------------------------------

update_tree = function(y, # Target variable
                       X, # Feature matrix
                       type = c('grow',   # Grow existing tree
                                'prune',  # Prune existing tree
                                'change', # Change existing tree - change split variable and value for an internal node
                                'swap'),  # Swap existing tree - swap splitting rules for two pairs of terminal nodes
                       curr_tree,         # The current set of trees (not required if type is stump)
                       node_min_size,     # The minimum size of a node to grow
                       s,                 # probability vector to be used during the growing process
                       coef_prior,        # Prior distribution for splitting coefficients
                       coef_hyperprior,   # hyperprior for splitting coefficients
                       hyp_par_list,       # list of hyperparameters for splitting coefficients
                       threshold_prior,    # splitting value  prior
                       coef_norm_hyperprior,
                       norm_unit_sphere
                       )
  {

  # Call the appropriate function to get the new tree
  new_tree = switch(type,
                    grow = grow_tree(X, y, curr_tree, node_min_size, s,
                                     coef_prior, coef_hyperprior, hyp_par_list, threshold_prior, coef_norm_hyperprior, norm_unit_sphere),
                    prune = prune_tree(X, y, curr_tree,
                                       coef_prior, coef_hyperprior, hyp_par_list, threshold_prior, coef_norm_hyperprior),
                    change = change_tree(X, y, curr_tree, node_min_size,
                                         coef_prior, coef_hyperprior, hyp_par_list, threshold_prior, coef_norm_hyperprior, norm_unit_sphere),
                    swap = swap_tree(X, y, curr_tree, node_min_size))

  # Return the new tree
  return(new_tree)

} # End of update_tree function

# Grow_tree function ------------------------------------------------------

grow_tree = function(X, y, curr_tree, node_min_size, s,
                     coef_prior,
                     coef_hyperprior,
                     hyp_par_list,
                     threshold_prior,
                     coef_norm_hyperprior,
                     norm_unit_sphere) {

  # Set up holder for new tree
  new_tree = curr_tree

  # Get the list of terminal nodes
  terminal_nodes = as.numeric(which(new_tree$tree_matrix[,'terminal'] == 1))

  # Find terminal node sizes
  terminal_node_size = as.numeric(new_tree$tree_matrix[terminal_nodes,'node_size'])

  if (all(terminal_node_size < 2 * node_min_size)) {
    # curr_tree$var = 0
    curr_tree$var_update = 0
    if(hyp_par_list$split_mix){
      new_tree$node_updated <- NA
      new_tree$c_ind_new <- NA
    }

    curr_tree$count_for_update <- rep(0, ncol(X))
    if(coef_norm_hyperprior == "varying"){
      new_tree$coef_for_sum_vec <- rep(0, ncol(X))
    }
    if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
      new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
    }
    if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
      new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
    }
    return(curr_tree)
  }

  available_values = NULL
  max_bad_trees = 10
  count_bad_trees = 0
  bad_trees = TRUE

  while (bad_trees ){

    # Set up holder for new tree
    new_tree = curr_tree

    # Add two extra rows to the tree in question
    new_tree$tree_matrix = rbind(new_tree$tree_matrix,
                                 c(1, NA, NA, NA, # NA,
                                   NA, NA, NA,
                                   rep(NA,ncol(X))), # Make sure they're both terminal
                                 c(1, NA, NA, NA, # NA,
                                   NA, NA, NA,
                                   rep(NA,ncol(X))))


    # Choose a random terminal node to split
    node_to_split = sample(terminal_nodes, 1,
                           prob = as.integer(terminal_node_size >= 2*node_min_size)) # Choose which node to split, set prob to zero for any nodes that are too small

    # # Choose a split variable uniformly from all columns (the first one is the intercept)
    # split_variable = sample(1:ncol(X), 1, prob = s)

    # PROPOSE COEFFICIENTS FOR OBLIQUE SPLIT
    # THIS DEPENDS ON THE PRIOR
    split_coefs <- rep(NA, ncol(X))



    if(hyp_par_list$split_mix){

      if(coef_hyperprior %in% c("univariate_normal_betabinomial_theta_j",
                                "univariate_normal_betabinomial_theta_j_sigma_j")){

        # first sample the cluster

        # c_ind <- sample(1:hyp_par_list$num_clust, size = 1, replace = FALSE,
        #                     prob = hyp_par_list$clust_probs )

        c_ind <- sample(hyp_par_list$num_clust, size = 1, replace = FALSE,
                        prob = hyp_par_list$clust_probs )


        gamma_vec <- rep(0, ncol(X))
        while(all(gamma_vec ==0)){ # cannot propose a split without nonzero coefficients

          gamma_vec <- rep(NA, ncol(X))
          for(j in 1:ncol(X)){
            gamma_vec[j] <- rbinom(1,
                                   size = 1,
                                   prob =  hyp_par_list$theta_list[[c_ind]][j])
            if(gamma_vec[j] ==0){
              split_coefs[j] <- 0
            }else{
              split_coefs[j] <- rnorm( n = 1,
                                       mean = hyp_par_list$beta_bar_list[[c_ind]][j],
                                       sd =  sqrt(hyp_par_list$sigma2_beta_list[[c_ind]][j]))
            }
          }
        }
      }

    }else{
      # can vectorize this to speed up
      if(coef_prior == "univariate_normal"){
        if( coef_hyperprior == "none"){
          for(j in 1:ncol(X)){
            split_coefs[j] <- rnorm( n = 1,
                                     mean = 0,
                                     sd =  sqrt(hyp_par_list$sigma2_beta_vec[j]))
          }
        }
        if(coef_hyperprior == "univariate_normal"){
          for(j in 1:ncol(X)){
            split_coefs[j] <- rnorm( n = 1,
                                     mean = hyp_par_list$beta_bar_vec[j],
                                     sd =  sqrt(hyp_par_list$sigma2_beta_vec[j]))
          }

        }

        if((coef_hyperprior == "univariate_normal_fixed_binomial") |
           (coef_hyperprior == "univariate_normal_betabinomial_onetheta")  ){
          # sample variables to include

          gamma_vec <- rep(0, ncol(X))
          while(all(gamma_vec ==0)){ # cannot propose a split without nonzero coefficients

            gamma_vec <- rbinom(ncol(X), size = 1, prob =  hyp_par_list$theta)
            if(any(is.na(gamma_vec))){
              print("hyp_par_list$theta = ")
              print(hyp_par_list$theta)
              print("ncol(X) = ")
              print(ncol(X))
              stop("gamma_vec NA")
            }

            split_coefs[gamma_vec == 0] <- 0

            include_inds <- which(gamma_vec == 1)
            for(j in include_inds){
              split_coefs[j] <- rnorm( n = 1,
                                       mean = hyp_par_list$beta_bar_vec[j],
                                       sd =  sqrt(hyp_par_list$sigma2_beta_vec[j]))
            }
          }
        }

        if(coef_hyperprior %in% c("univariate_normal_betabinomial_theta_j",
                                  "univariate_normal_betabinomial_theta_j_sigma_j")){
          gamma_vec <- rep(0, ncol(X))
          while(all(gamma_vec ==0)){ # cannot propose a split without nonzero coefficients

            gamma_vec <- rep(NA, ncol(X))
            for(j in 1:ncol(X)){
              gamma_vec[j] <- rbinom(1,
                                     size = 1,
                                     prob =  hyp_par_list$theta_vec[j])
              if(gamma_vec[j] ==0){
                split_coefs[j] <- 0
              }else{
                split_coefs[j] <- rnorm( n = 1,
                                         mean = hyp_par_list$beta_bar_vec[j],
                                         sd =  sqrt(hyp_par_list$sigma2_beta_vec[j]))
              }
            }
          }
        }
      } # end coef_prior == "univariate_normal" if statement




      if(coef_prior == "simplex"){
        if( coef_hyperprior == "none"){
          # no hyperprior, just draw from the Dirichlet prior

          split_coefs <- t(rdirichlet(n = 1,
                                    alpha = hyp_par_list$alpha_simplex * hyp_par_list$xi_vec))

        }
        if( coef_hyperprior == "simplex_fixed_beta_binomial_theta_j"){
          # no hyperprior, just draw from the Dirichlet prior
          gamma_vec <- rep(0, ncol(X))
          while(all(gamma_vec ==0)){ # cannot propose a split without nonzero coefficients

            for(j in 1:ncol(X)){
              gamma_vec[j] <- rbinom(1,
                                     size = 1,
                                     prob =  hyp_par_list$theta_vec[j])
            }

            temp_alphs <- hyp_par_list$alpha_simplex * sum(gamma_vec)/ ncol(X)
            temp_xi_vec <-  hyp_par_list$xi_vec[ gamma_vec == 1]/ sum(hyp_par_list$xi_vec[ gamma_vec == 1])

            # if(gamma_vec[j] ==0){
              split_coefs[ gamma_vec == 0] <- 0
            # }else{
              split_coefs[ gamma_vec == 1] <- t(rdirichlet(n = 1,
                                          alpha = temp_alphs * temp_xi_vec))
            # }
          }
        }


        if( coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j"){
          # no hyperprior, just draw from the Dirichlet prior
          gamma_vec <- rep(0, ncol(X))
          while(all(gamma_vec ==0)){ # cannot propose a split without nonzero coefficients

            for(j in 1:ncol(X)){
              gamma_vec[j] <- rbinom(1,
                                     size = 1,
                                     prob =  hyp_par_list$theta_vec[j])
            }

            split_coefs <- t(rdirichlet(n = 1,
                                        alpha = hyp_par_list$alpha_simplex * hyp_par_list$xi_vec))

            split_coefs[ gamma_vec == 0] <- -1* split_coefs[ gamma_vec == 0]

            # }
          }
        }


        if( coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
          # no hyperprior, just draw from the Dirichlet prior
          gamma_vec <- rep(0, ncol(X))
          while(all(gamma_vec ==0)){ # cannot propose a split without nonzero coefficients

            for(j in 1:ncol(X)){
              gamma_vec[j] <- rbinom(1,
                                     size = 1,
                                     prob =  hyp_par_list$theta_vec[j])
            }

            split_coefs <- t(rdirichlet(n = 1,
                                        alpha = hyp_par_list$xi_vec))

            split_coefs[ gamma_vec == 0] <- -1* split_coefs[ gamma_vec == 0]

            # }
          }
        }





        if(coef_hyperprior == "simplex_fixed_Dir_trinomial_theta_j"){
          # no hyperprior, just draw from the Dirichlet prior
          delta_vec <- rep(0, ncol(X))
          while(all(delta_vec ==0)){ # cannot propose a split without nonzero coefficients

            for(j in 1:ncol(X)){
              delta_vec[j] <- sample(c(-1,0,1),
                                     size = 1,
                                     prob =  c(hyp_par_list$theta_min1_vec[j],
                                               hyp_par_list$theta_0_vec[j],
                                               hyp_par_list$theta_1_vec[j]))
            }

            temp_alphs <- hyp_par_list$alpha_simplex * sum( delta_vec != 0)/ ncol(X)
            temp_xi_vec <-  hyp_par_list$xi_vec[ delta_vec != 0]/ sum(hyp_par_list$xi_vec[ delta_vec != 0])

            # if(gamma_vec[j] ==0){
            split_coefs[ delta_vec == 0] <- 0
            # }else{
            split_coefs[ delta_vec != 0] <- t(rdirichlet(n = 1,
                                                         alpha = temp_alphs * temp_xi_vec))

            split_coefs[ delta_vec == -1] <- -1* split_coefs[ delta_vec == -1]
            # }
          }

        }

      } # end of simplex coef draw if statement

    } # end else statement not mixture



    if(norm_unit_sphere){
      split_coefs <- split_coefs/sum(split_coefs^2)
    }


    if(any(is.na(split_coefs))){
      stop("NA split_coefs values")
    }



    # # Alternatively follow BARTMachine and choose a split value using sample on the internal values of the available
    # available_values = sort(unique(X[new_tree$node_indices == node_to_split,
    #                                  split_variable]))
    #
    # if(length(available_values) == 1){
    #   split_value = available_values[1]
    # } else if (length(available_values) == 2){
    #   split_value = available_values[2]
    # }  else {
    #   # split_value = sample(available_values[-c(1,length(available_values))], 1)
    #   split_value = resample(available_values[-c(1,length(available_values))])
    # }



    if(threshold_prior != "continuous_minus_plus_1"){
      # calculate all linear combination values in the relevant node
      # if( (ncol(X[new_tree$node_indices == node_to_split, ]) != length(split_coefs) ) | (nrow(X[new_tree$node_indices == node_to_split, ]) ==0)  ){
      #   print("split_coefs = ")
      #   print(split_coefs)
      #
      #   print("X[new_tree$node_indices == node_to_split, ] = ")
      #   print(X[new_tree$node_indices == node_to_split, ])
      #
      # }
      # print("split_coefs = ")
      # print(split_coefs)
      #
      # print("X[new_tree$node_indices == node_to_split, ] = ")
      # print(X[new_tree$node_indices == node_to_split, ])


      # lincomb_values <- unique(as.vector(X[new_tree$node_indices == node_to_split, ] %*% split_coefs))

      countout <- collapse::fcount(as.vector(X[new_tree$node_indices == node_to_split, ] %*% split_coefs), sort = TRUE)
      lincomb_values <- countout$x


    }

    if(threshold_prior == "discrete_uniform"){
      # Alternatively follow BARTMachine and choose a split value using sample on the internal values of the available
      # available_values = sort(lincomb_values, na.last = TRUE)

      # if(length(available_values) == 0){
      #   print(" lincomb_values = ")
      #   print(lincomb_values)
      #   print(" X[new_tree$node_indices == node_to_split, ] = ")
      #   print(X[new_tree$node_indices == node_to_split, ])
      #
      #   print(" split_coefs = ")
      #   print(split_coefs)
      #
      #   print("hyp_par_list$beta_bar_vec = ")
      #   print(hyp_par_list$beta_bar_vec)
      #
      #   print("hyp_par_list$sigma2_beta_vec = ")
      #   print(hyp_par_list$sigma2_beta_vec)
      #
      #   stop("Line 302 (length(available_values) == 0")
      #
      # }


      # if(length(lincomb_values) == 1){
      #   split_value = lincomb_values[1]
      # } else if (length(lincomb_values) == 2){
      #   split_value = max(lincomb_values)
      # }  else {
      #   # split_value = sample(available_values[-c(1,length(available_values))], 1)
      #   # split_value = resample(available_values[-c(1,length(available_values))])
      #   # split_value = sample(x = available_values[2:(length(available_values)-1)],size = 1)
      #   split_value = sample(x = lincomb_values[-c(which.min(lincomb_values), which.max(lincomb_values))],size = 1)
      #
      # }
      if(length(lincomb_values) == 1){
        # no point in splitting if only one value
        n_bad_trees = n_bad_trees + 1
        if(n_bad_trees >= max_bad_trees) {
          # print(" reached max_bad_trees = ")
          # curr_tree$var = 0
          curr_tree$var_update = 0
          if(hyp_par_list$split_mix){
            new_tree$node_updated <- NA
            new_tree$c_ind_new <- NA
          }

          curr_tree$count_for_update <- rep(0, ncol(X))
          if(coef_norm_hyperprior == "varying"){
            new_tree$coef_for_sum_vec <- rep(0, ncol(X))
          }
          if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
            new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
          }
          if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
            new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
          }
          return(curr_tree)
        }else{
          next
        }

        # new_split_value = available_values[1]
      } else if (length(lincomb_values) == 2){

        if(any(countout$N < node_min_size)){
          n_bad_trees = n_bad_trees + 1
          if(n_bad_trees >= max_bad_trees) {
            # print(" reached max_bad_trees = ")

            # curr_tree$var = 0
            curr_tree$var_update = 0
            if(hyp_par_list$split_mix){
              new_tree$node_updated <- NA
              new_tree$c_ind_new <- NA
            }

            curr_tree$count_for_update <- rep(0, ncol(X))
            if(coef_norm_hyperprior == "varying"){
              new_tree$coef_for_sum_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
              new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
              new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
            }
            return(curr_tree)
          }else{
            next
          }
        }

        # split_value = lincomb_values[2]
        # if(splitting_rules == "continuous"){
        #   split_value = runif(1,lincomb_values[1],lincomb_values[2])
        # }else{
          split_value = lincomb_values[2]
        # }

      }  else {

        # find smallest and largest split values with less than minimum node size left and right
        runsum <- 0
        # min_ind <- 0
        for(val_ind in 1:length(lincomb_values)){
          runsum <- runsum + countout$N[val_ind]
          if(runsum >= node_min_size){
            min_ind <- val_ind +1
            break
          }
        }
        runsum <- 0
        # max_ind <- 0
        for(val_ind in length(lincomb_values):1){
          runsum <- runsum + countout$N[val_ind]
          if(runsum >= node_min_size){
            max_ind <- val_ind
            break
          }
        }

        # if((min_ind > length(lincomb_values)) | max_ind == 0  ){
        if((min_ind > length(lincomb_values)) | (min_ind > max_ind)  ){
          n_bad_trees = n_bad_trees + 1
          if(n_bad_trees >= max_bad_trees) {
            # print(" reached max_bad_trees = ")

            # curr_tree$var = 0
            curr_tree$var_update = 0
            if(hyp_par_list$split_mix){
              new_tree$node_updated <- NA
              new_tree$c_ind_new <- NA
            }

            curr_tree$count_for_update <- rep(0, ncol(X))
            if(coef_norm_hyperprior == "varying"){
              new_tree$coef_for_sum_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
              new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
              new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
            }
            return(curr_tree)
          }else{
            next
          }
        }

        if(min_ind > max_ind){
          stop("min_ind > max_ind")
        }

        if(min_ind == max_ind){
          # split_value <- lincomb_values[min_ind]
          # if(splitting_rules == "continuous"){
          #   split_value <- runif(1,lincomb_values[min_ind-1],lincomb_values[min_ind])
          # }else{
            split_value <- lincomb_values[min_ind]
          # }
        }else{
          # split_value <- sample(lincomb_values[min_ind:max_ind],1)
          # if(splitting_rules == "continuous"){
          #   split_value = runif(1,lincomb_values[min_ind-1],lincomb_values[max_ind])
          # }else{
            split_value <- sample(lincomb_values[min_ind:max_ind],1)
          # }
          # new_split_value = runif(1,lincomb_values[min_ind],lincomb_values[max_ind])
        }
        # split_value = sample(available_values[-c(1,length(available_values))], 1)
        # split_value = resample(available_values[-c(1,length(available_values))])
        # split_value = sample(available_values[-c(1)],1)
        # split_value = resample(available_values[-c(1,length(available_values))])
        # split_value = runif(1,available_values[2],available_values[length(available_values)])
      }


    } # end discrete uniform prior

    if(threshold_prior == "continuous_min_max"){


      # tempmin <- min(lincomb_values)
      # tempmax <- max(lincomb_values)
      # split_value <- runif(n = 1, min = tempmin, max = tempmax)
      if(length(lincomb_values) == 1){
        # no point in splitting if only one value
        n_bad_trees = n_bad_trees + 1
        if(n_bad_trees >= max_bad_trees) {
          # print(" reached max_bad_trees = ")
          # curr_tree$var = 0
          curr_tree$var_update = 0
          if(hyp_par_list$split_mix){
            new_tree$node_updated <- NA
            new_tree$c_ind_new <- NA
          }

          curr_tree$count_for_update <- rep(0, ncol(X))
          if(coef_norm_hyperprior == "varying"){
            new_tree$coef_for_sum_vec <- rep(0, ncol(X))
          }
          if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
            new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
          }
          if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
            new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
          }
          return(curr_tree)
        }else{
          next
        }

        # new_split_value = available_values[1]
      } else if (length(lincomb_values) == 2){

        if(any(countout$N < node_min_size)){
          n_bad_trees = n_bad_trees + 1
          if(n_bad_trees >= max_bad_trees) {
            # print(" reached max_bad_trees = ")

            # curr_tree$var = 0
            curr_tree$var_update = 0
            if(hyp_par_list$split_mix){
              new_tree$node_updated <- NA
              new_tree$c_ind_new <- NA
            }

            curr_tree$count_for_update <- rep(0, ncol(X))
            if(coef_norm_hyperprior == "varying"){
              new_tree$coef_for_sum_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
              new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
              new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
            }
            return(curr_tree)
          }else{
            next
          }
        }

        # split_value = lincomb_values[2]
        split_value = runif(1,lincomb_values[1],lincomb_values[2])


      }  else {

        # find smallest and largest split values with less than minimum node size left and right
        runsum <- 0
        # min_ind <- 0
        for(val_ind in 1:length(lincomb_values)){
          runsum <- runsum + countout$N[val_ind]
          if(runsum >= node_min_size){
            min_ind <- val_ind +1
            break
          }
        }
        runsum <- 0
        # max_ind <- 0
        for(val_ind in length(lincomb_values):1){
          runsum <- runsum + countout$N[val_ind]
          if(runsum >= node_min_size){
            max_ind <- val_ind
            break
          }
        }

        # if((min_ind > length(lincomb_values)) | max_ind == 0  ){
        if((min_ind > length(lincomb_values)) | (min_ind > max_ind)  ){
          n_bad_trees = n_bad_trees + 1
          if(n_bad_trees >= max_bad_trees) {
            # print(" reached max_bad_trees = ")

            # curr_tree$var = 0
            curr_tree$var_update = 0
            if(hyp_par_list$split_mix){
              new_tree$node_updated <- NA
              new_tree$c_ind_new <- NA
            }

            curr_tree$count_for_update <- rep(0, ncol(X))
            if(coef_norm_hyperprior == "varying"){
              new_tree$coef_for_sum_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
              new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
              new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
            }
            return(curr_tree)
          }else{
            next
          }
        }

        if(min_ind > max_ind){
          stop("min_ind > max_ind")
        }

        if(min_ind == max_ind){
          # split_value <- lincomb_values[min_ind]
          split_value <- runif(1,lincomb_values[min_ind-1],lincomb_values[min_ind])
        }else{
          # new_split_value <- sample(lincomb_values[min_ind:max_ind],1)
          # split_value = runif(1,lincomb_values[min_ind],lincomb_values[max_ind])
          split_value = runif(1,lincomb_values[min_ind-1],lincomb_values[max_ind])
        }
        # split_value = sample(available_values[-c(1,length(available_values))], 1)
        # split_value = resample(available_values[-c(1,length(available_values))])
        # new_split_value = sample(available_values[-c(1)],1)
        # split_value = resample(available_values[-c(1,length(available_values))])
        # split_value = runif(1,available_values[2],available_values[length(available_values)])
      }




    }


    if(threshold_prior == "continuous_minus_plus_1"){
      split_value <- runif(n = 1, min = -1, max = 1)
    }


    curr_parent = new_tree$tree_matrix[node_to_split, 'parent'] # Make sure to keep the current parent in there. Will be NA if at the root node
    new_tree$tree_matrix[node_to_split,1:5] = c(0, # Now not temrinal
                                                nrow(new_tree$tree_matrix) - 1, # child_left is penultimate row
                                                nrow(new_tree$tree_matrix),  # child_right is penultimate row
                                                curr_parent,
                                                # split_variable,
                                                split_value)

    new_tree$tree_matrix[node_to_split,8:ncol(new_tree$tree_matrix)] <- split_coefs

    #  Fill in the parents of these two nodes
    new_tree$tree_matrix[nrow(new_tree$tree_matrix),'parent'] = node_to_split
    new_tree$tree_matrix[nrow(new_tree$tree_matrix)-1,'parent'] = node_to_split

    # Now call the fill function on this tree
    new_tree = fill_tree_details(new_tree, X)

    # Store the covariate name to use it to update the Dirichlet prior of Linero (2016).
    # new_tree$var = split_variable
    new_tree$var_update = 1

    if(hyp_par_list$split_mix){
      new_tree$node_updated <- node_to_split
      new_tree$c_ind_new <- c_ind
    }

    if(coef_hyperprior %in% c("univariate_normal_fixed_binomial",
                              "univariate_normal_betabinomial_onetheta",
                              "univariate_normal_betabinomial_theta_j",
                              "univariate_normal_betabinomial_theta_j_sigma_j",
                              "simplex_fixed_beta_binomial_theta_j",
                              "simplex_fixed_Dir_binomial_plusminus_theta_j",
                              "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j")  ){

      new_tree$count_for_update = gamma_vec

    }

    if(coef_hyperprior == "simplex_fixed_Dir_trinomial_theta_j"){
      new_tree$count_min1_for_update = 1*(delta_vec == -1)
      new_tree$count_0_for_update = 1*(delta_vec == 0)
      new_tree$count_1_for_update = 1*(delta_vec == 1)
    }

    if(coef_norm_hyperprior == "varying"){
      new_tree$coef_for_sum_vec <- split_coefs
    }

    if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
      new_tree$coef_for_sumsq_vec <- gamma_vec*(split_coefs)^2
    }
    if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
      new_tree$coef_for_logsum_vec <- log(abs(split_coefs ))
    }




    # Check for bad tree
    if(any(as.numeric(new_tree$tree_matrix[,'node_size']) <= node_min_size)) {
      count_bad_trees = count_bad_trees + 1
    } else {
      bad_trees = FALSE
    }

    if(count_bad_trees == max_bad_trees) {
      # curr_tree$var = 0
      curr_tree$var_update = 0
      if(hyp_par_list$split_mix){
        new_tree$node_updated <- NA
        new_tree$c_ind_new <- NA
      }

      curr_tree$count_for_update <- rep(0, ncol(X))
      if(coef_norm_hyperprior == "varying"){
        new_tree$coef_for_sum_vec <- rep(0, ncol(X))
      }
      if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
        new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
      }
      if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
        new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
      }
      return(curr_tree)
    }
  }
  # Return new_tree
  return(new_tree)

} # End of grow_tree function

# Prune_tree function -----------------------------------------------------

prune_tree = function(X, y, curr_tree,
                      coef_prior, coef_hyperprior, hyp_par_list, threshold_prior,
                      coef_norm_hyperprior) {

  # Create placeholder for new tree
  new_tree = curr_tree

  if(nrow(new_tree$tree_matrix) == 1) { # No point in pruning a stump!
    # new_tree$var = 0
    new_tree$var_update = 0

    if(hyp_par_list$split_mix){
      new_tree$nodes_pruned <- NA
      new_tree$parent_pruned <- NA

      # new_tree$c_ind_new <- NA
    }


    new_tree$count_for_update <- rep(0, ncol(X))
    if(coef_norm_hyperprior == "varying"){
      new_tree$coef_for_sum_vec <- rep(0, ncol(X))
    }
    if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
      new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
    }
    if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
      new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
    }
    return(new_tree)
  }

  # Get the list of terminal nodes
  terminal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 1)

  # Pick a random terminal node to prune
  # ONLY PICK NODES WHERE BOTH LEFT AND RIGHT CHILD ARE TERMINAL
  bad_node_to_prune = TRUE # Assume a bad node pick
  while(bad_node_to_prune) {

    # Choose a random terminal node
    node_to_prune = sample(terminal_nodes, 1)

    # Find the parent of this terminal node
    parent_pick = as.numeric(new_tree$tree_matrix[node_to_prune, 'parent'])
    # var_pruned_nodes = as.numeric(new_tree$tree_matrix[parent_pick, 'split_variable'])
    pruned_coefs <- as.numeric(new_tree$tree_matrix[parent_pick, 8:ncol(new_tree$tree_matrix)])
    gamma_pruned_nodes = 1*(pruned_coefs != 0 )

    # Get the two children of this parent
    child_left = as.numeric(new_tree$tree_matrix[parent_pick, 'child_left'])
    child_right = as.numeric(new_tree$tree_matrix[parent_pick, 'child_right'])

    # See whether either are terminal
    child_left_terminal = as.numeric(new_tree$tree_matrix[child_left, 'terminal'])
    child_right_terminal = as.numeric(new_tree$tree_matrix[child_right, 'terminal'])

    # If both are terminal then great
    if( (child_left_terminal == 1) & (child_right_terminal == 1) ) {
      bad_node_to_prune = FALSE # Have chosen a pair of terminal nodes so exist while loop
    }

  }# End of bad node to prune while loop

  # Delete these two rows from the tree matrix
  new_tree$tree_matrix = new_tree$tree_matrix[-c(child_left,child_right),,
                                              drop = FALSE]
  # Make this node terminal again with no children or split values
  new_tree$tree_matrix[parent_pick,c('terminal',
                                     'child_left',
                                     'child_right',
                                     # 'split_variable',
                                     'split_value')] = c(1, NA, NA, #NA,
                                                         NA)

  new_tree$tree_matrix[parent_pick,8:ncol(new_tree$tree_matrix)] = rep(NA, ncol(X))


  # If we're back to a stump no need to call fill_tree_details
  if(nrow(new_tree$tree_matrix) == 1) {
    # new_tree$var = var_pruned_nodes
    new_tree$var_update = 1

    if(hyp_par_list$split_mix){
      new_tree$nodes_pruned <- c(child_left,child_right)
      new_tree$parent_pruned <- parent_pick

      # new_tree$c_ind_new <- c_ind
    }

    if(coef_hyperprior %in% c("univariate_normal_fixed_binomial",
                              "univariate_normal_betabinomial_onetheta",
                              "univariate_normal_betabinomial_theta_j",
                              "univariate_normal_betabinomial_theta_j_sigma_j",
                              "simplex_fixed_beta_binomial_theta_j")  ){

      new_tree$count_for_update = gamma_pruned_nodes

    }

    if(coef_hyperprior %in% c("simplex_fixed_Dir_binomial_plusminus_theta_j",
                              "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j")){
      new_tree$count_for_update <- 1*(pruned_coefs >= 0 )
    }

    if(coef_hyperprior == "simplex_fixed_Dir_trinomial_theta_j"){
      new_tree$count_min1_for_update = 1*(pruned_coefs < 0)
      new_tree$count_0_for_update = 1*(pruned_coefs == 0)
      new_tree$count_1_for_update = 1*(pruned_coefs > 0)
    }



    if(coef_norm_hyperprior == "varying"){
      new_tree$coef_for_sum_vec <- pruned_coefs
    }
    if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
      new_tree$coef_for_sumsq_vec <- gamma_pruned_nodes*(pruned_coefs )^2
    }
    if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
      new_tree$coef_for_logsum_vec <- log(abs(pruned_coefs ))
    }

    new_tree$node_indices = rep(1, length(y))
  } else {
    # If we've removed some nodes from the middle we need to re-number all the child_left and child_right values - the parent values will still be correct
    if(node_to_prune <= nrow(new_tree$tree_matrix)) { # Only need do this if we've removed some observations from the middle of the tree matrix
      # If you're pruning any nodes which affect parent indices further down the tree then make sure to shift the parent values
      bad_parents = which(as.numeric(new_tree$tree_matrix[,'parent'])>=node_to_prune)
      # Shift them back because you have removed two rows
      new_tree$tree_matrix[bad_parents,'parent'] = as.numeric(new_tree$tree_matrix[bad_parents,'parent']) - 2

      for(j in node_to_prune:nrow(new_tree$tree_matrix)) {
        # Find the current parent
        curr_parent = as.numeric(new_tree$tree_matrix[j,'parent'])
        # Find both the children of this node
        curr_children = which(as.numeric(new_tree$tree_matrix[,'parent']) == curr_parent)
        # Input these children back into the parent
        new_tree$tree_matrix[curr_parent,c('child_left','child_right')] = sort(curr_children, na.last = TRUE)
      } # End for loop of correcting parents and children
    } # End if statement to fill in tree details

    # Call the fill function on this tree
    new_tree = fill_tree_details(new_tree, X)

    # Store the covariate name that was used in the splitting rule of the terminal nodes that were just pruned
    # new_tree$var = var_pruned_nodes
    new_tree$var_update = 1
    if(hyp_par_list$split_mix){
      new_tree$nodes_pruned <- c(child_left,child_right) # use indices from old tree matrix, so this is not curr_children

      new_tree$parent_pruned <- parent_pick
      # new_tree$c_ind_new <- c_ind
    }

    if(coef_hyperprior %in% c("univariate_normal_fixed_binomial",
                              "univariate_normal_betabinomial_onetheta",
                              "univariate_normal_betabinomial_theta_j",
                              "univariate_normal_betabinomial_theta_j_sigma_j",
                              "simplex_fixed_beta_binomial_theta_j")  ){
      new_tree$count_for_update = gamma_pruned_nodes
    }
    if(coef_hyperprior %in% c("simplex_fixed_Dir_binomial_plusminus_theta_j",
                              "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j")){
      new_tree$count_for_update <- 1*(pruned_coefs >= 0 )
    }

    if(coef_hyperprior == "simplex_fixed_Dir_trinomial_theta_j"){
      new_tree$count_min1_for_update = 1*(pruned_coefs < 0)
      new_tree$count_0_for_update = 1*(pruned_coefs == 0)
      new_tree$count_1_for_update = 1*(pruned_coefs > 0)
    }

    if(coef_norm_hyperprior == "varying"){
      new_tree$coef_for_sum_vec <- pruned_coefs
    }
    if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
      new_tree$coef_for_sumsq_vec <- gamma_pruned_nodes*(pruned_coefs )^2
    }
    if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
      new_tree$coef_for_logsum_vec <- log(abs(pruned_coefs ))
    }

  }

  # Return new_tree
  return(new_tree)

} # End of prune_tree function

# change_tree function ----------------------------------------------------

change_tree = function(X, y, curr_tree, node_min_size,
                       coef_prior,
                       coef_hyperprior,
                       hyp_par_list,
                       threshold_prior,
                       coef_norm_hyperprior,
                       norm_unit_sphere) {

  # Change a node means change out the split value and split variable of an internal node. Need to make sure that this does now produce a bad tree (i.e. zero terminal nodes)

  # If current tree is a stump nothing to change
  if(nrow(curr_tree$tree_matrix) == 1) {
    # curr_tree$var = c(0, 0)
    curr_tree$var_update = 0
    curr_tree$count_for_update <- rep(0, ncol(X))

    if(hyp_par_list$split_mix){
      new_tree$node_updated <- NA
      new_tree$c_ind_new <- NA
    }

    if(coef_norm_hyperprior == "varying"){
      new_tree$coef_for_sum_vec <- rep(0, ncol(X))
    }
    if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
      new_tree$coef_for_sumsq_vec <-  rep(0, ncol(X))
    }
    if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
      new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
    }
    return(curr_tree)
  }

  # Create a holder for the new tree
  new_tree = curr_tree

  # Need to get the internal nodes
  internal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 0)
  terminal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 1)

  # Create a while loop to get good trees
  # Create a counter to stop after a certain number of bad trees
  max_bad_trees = 2
  count_bad_trees = 0
  bad_trees = TRUE
  while(bad_trees) {
    # Re-set the tree
    new_tree = curr_tree

    # choose an internal node to change
    node_to_change = sample(internal_nodes, 1)

    # Get the covariate that will be changed
    # var_changed_node = as.numeric(new_tree$tree_matrix[node_to_change, 'split_variable'])

    # Use the get_children function to get all the children of this node
    all_children = get_children(new_tree$tree_matrix, node_to_change)

    # Now find all the nodes which match these children
    use_node_indices = !is.na(match(new_tree$node_indices, all_children))

    # Create new split variable and value based on ignorance
    # then check this doesn't give a bad tree

    # available_values = NULL
    #
    # new_split_variable = sample(1:ncol(X), 1)
    #
    # available_values = sort(unique(X[use_node_indices,
    #                                  new_split_variable]))
    #
    # if (length(available_values) == 1){
    #   new_split_value = available_values[1]
    #   # new_tree$var = c(var_changed_node, new_split_variable)
    # } else if (length(available_values) == 2){
    #   new_split_value = available_values[2]
    #   # new_tree$var = c(var_changed_node, new_split_variable)
    # } else {
    #   # new_split_value = sample(available_values[-c(1,length(available_values))], 1)
    #   new_split_value = resample(available_values[-c(1,length(available_values))])
    # }



    # PROPOSE COEFFICIENTS FOR OBLIQUE SPLIT
    # THIS DEPENDS ON THE PRIOR
    new_split_coefs <- rep(NA, ncol(X))


    if(hyp_par_list$split_mix){
      if(coef_hyperprior %in% c("univariate_normal_betabinomial_theta_j",
                                "univariate_normal_betabinomial_theta_j_sigma_j")){


        # c_ind <- sample(1:hyp_par_list$num_clust, size = 1, replace = FALSE,
        #                 prob = hyp_par_list$clust_probs )
        c_ind <- sample.int(hyp_par_list$num_clust, size = 1, replace = FALSE,
                        prob = hyp_par_list$clust_probs )

        gamma_vec <- rep(0, ncol(X))
        while(all(gamma_vec ==0)){ # cannot propose a split without nonzero coefficients
          gamma_vec <- rep(NA, ncol(X))
          for(j in 1:ncol(X)){
            gamma_vec[j] <- rbinom(1,
                                   size = 1,
                                   prob =  hyp_par_list$theta_list[[c_ind]][j])
            if(gamma_vec[j] ==0){
              new_split_coefs[j] <- 0
            }else{
              new_split_coefs[j] <- rnorm( n = 1,
                                           mean = hyp_par_list$beta_bar_list[[c_ind]][j],
                                           sd =  sqrt(hyp_par_list$sigma2_beta_list[[c_ind]][j]))
            }
          }

          coefs_changed_node <- as.numeric(new_tree$tree_matrix[node_to_change, 8:ncol(new_tree$tree_matrix)])
          gamma_changed_node <- 1*( coefs_changed_node != 0)
        }
      }

    }else{
      # can vectorize this to speed up
      if(coef_prior == "univariate_normal"){
        if( coef_hyperprior == "none"){
          for(j in 1:ncol(X)){
            new_split_coefs[j] <- rnorm( n = 1,
                                     mean = 0,
                                     sd =  sqrt(hyp_par_list$sigma2_beta_vec[j]))
          }
        }
        if(coef_hyperprior == "univariate_normal"){
          for(j in 1:ncol(X)){
            new_split_coefs[j] <- rnorm( n = 1,
                                     mean = hyp_par_list$beta_bar_vec[j],
                                     sd =  sqrt(hyp_par_list$sigma2_beta_vec[j]))
          }
        }

        if((coef_hyperprior == "univariate_normal_fixed_binomial") |
           (coef_hyperprior == "univariate_normal_betabinomial_onetheta")  ){
          # sample variables to include
          gamma_vec <- rep(0, ncol(X))
          while(all(gamma_vec ==0)){ # cannot propose a split without nonzero coefficients
            gamma_vec <- rbinom(ncol(X), size = 1, prob =  hyp_par_list$theta)
            if(any(is.na(gamma_vec))){
              print("hyp_par_list$theta = ")
              print(hyp_par_list$theta)
              print("ncol(X) = ")
              print(ncol(X))
              stop("gamma_vec NA")
            }

            new_split_coefs[gamma_vec == 0] <- 0
            include_inds <- which(gamma_vec == 1)
            for(j in include_inds){
              new_split_coefs[j] <- rnorm( n = 1,
                                       mean = hyp_par_list$beta_bar_vec[j],
                                       sd =  sqrt(hyp_par_list$sigma2_beta_vec[j]))
            }

            coefs_changed_node <- as.numeric(new_tree$tree_matrix[node_to_change, 8:ncol(new_tree$tree_matrix)])
            gamma_changed_node <- 1*(coefs_changed_node != 0)
          }

        }
        if(coef_hyperprior %in% c("univariate_normal_betabinomial_theta_j",
                                  "univariate_normal_betabinomial_theta_j_sigma_j")){
          gamma_vec <- rep(0, ncol(X))
          while(all(gamma_vec ==0)){ # cannot propose a split without nonzero coefficients
            gamma_vec <- rep(NA, ncol(X))
            for(j in 1:ncol(X)){
              gamma_vec[j] <- rbinom(1,
                                     size = 1,
                                     prob =  hyp_par_list$theta_vec[j])
              if(gamma_vec[j] ==0){
                new_split_coefs[j] <- 0
              }else{
                new_split_coefs[j] <- rnorm( n = 1,
                                         mean = hyp_par_list$beta_bar_vec[j],
                                         sd =  sqrt(hyp_par_list$sigma2_beta_vec[j]))
              }
            }

            coefs_changed_node <- as.numeric(new_tree$tree_matrix[node_to_change, 8:ncol(new_tree$tree_matrix)])
            gamma_changed_node <- 1*( coefs_changed_node != 0)
          }
        }
      } # end coef_prior == "univariate_normal" if statement



      if(coef_prior == "simplex"){
        if( coef_hyperprior == "none"){
          # no hyperprior, just draw from the Dirichlet prior

          new_split_coefs <- t(rdirichlet(n = 1,
                                    alpha = hyp_par_list$alpha_simplex * hyp_par_list$xi_vec))

        }
        if( coef_hyperprior == "simplex_fixed_beta_binomial_theta_j"){
          # no hyperprior, just draw from the Dirichlet prior
          gamma_vec <- rep(0, ncol(X))
          while(all(gamma_vec ==0)){ # cannot propose a split without nonzero coefficients

            for(j in 1:ncol(X)){
              gamma_vec[j] <- rbinom(1,
                                     size = 1,
                                     prob =  hyp_par_list$theta_vec[j])
            }

            temp_alphs <- hyp_par_list$alpha_simplex * sum(gamma_vec)/ ncol(X)
            temp_xi_vec <-  hyp_par_list$xi_vec[ gamma_vec == 1]/ sum(hyp_par_list$xi_vec[ gamma_vec == 1])

            # if(gamma_vec[j] ==0){
            new_split_coefs[ gamma_vec == 0] <- 0
            # }else{
            new_split_coefs[ gamma_vec == 1] <- t(rdirichlet(n = 1,
                                                         alpha = temp_alphs * temp_xi_vec))
            # }
            coefs_changed_node <- as.numeric(new_tree$tree_matrix[node_to_change, 8:ncol(new_tree$tree_matrix)])
            gamma_changed_node <- 1*( coefs_changed_node != 0)
          }

        }

        if( coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j"){
          # no hyperprior, just draw from the Dirichlet prior
          gamma_vec <- rep(0, ncol(X))
          while(all(gamma_vec ==0)){ # cannot propose a split without nonzero coefficients

            for(j in 1:ncol(X)){
              gamma_vec[j] <- rbinom(1,
                                     size = 1,
                                     prob =  hyp_par_list$theta_vec[j])
            }

            new_split_coefs <- t(rdirichlet(n = 1,
                                        alpha = hyp_par_list$alpha_simplex * hyp_par_list$xi_vec))

            new_split_coefs[ gamma_vec == 0] <- -1* new_split_coefs[ gamma_vec == 0]
            coefs_changed_node <- as.numeric(new_tree$tree_matrix[node_to_change, 8:ncol(new_tree$tree_matrix)])

            # }
          }
        }

        if( coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
          # no hyperprior, just draw from the Dirichlet prior
          gamma_vec <- rep(0, ncol(X))
          while(all(gamma_vec ==0)){ # cannot propose a split without nonzero coefficients

            for(j in 1:ncol(X)){
              gamma_vec[j] <- rbinom(1,
                                     size = 1,
                                     prob =  hyp_par_list$theta_vec[j])
            }

            new_split_coefs <- t(rdirichlet(n = 1,
                                            alpha = hyp_par_list$xi_vec))

            new_split_coefs[ gamma_vec == 0] <- -1* new_split_coefs[ gamma_vec == 0]
            coefs_changed_node <- as.numeric(new_tree$tree_matrix[node_to_change, 8:ncol(new_tree$tree_matrix)])

            # }
          }
        }



        if(coef_hyperprior == "simplex_fixed_Dir_trinomial_theta_j"){
          # no hyperprior, just draw from the Dirichlet prior
          delta_vec <- rep(0, ncol(X))
          while(all(delta_vec ==0)){ # cannot propose a split without nonzero coefficients

            for(j in 1:ncol(X)){
              delta_vec[j] <- sample(c(-1,0,1),
                                     size = 1,
                                     prob =  c(hyp_par_list$theta_min1_vec[j],
                                               hyp_par_list$theta_0_vec[j],
                                               hyp_par_list$theta_1_vec[j]))
            }

            temp_alphs <- hyp_par_list$alpha_simplex * sum( delta_vec != 0)/ ncol(X)
            temp_xi_vec <-  hyp_par_list$xi_vec[ delta_vec != 0]/ sum(hyp_par_list$xi_vec[ delta_vec != 0])

            # if(gamma_vec[j] ==0){
            new_split_coefs[ delta_vec == 0] <- 0
            # }else{
            new_split_coefs[ delta_vec != 0] <- t(rdirichlet(n = 1,
                                                         alpha = temp_alphs * temp_xi_vec))

            new_split_coefs[ delta_vec == -1] <- -1* new_split_coefs[ delta_vec == -1]
            # }

            coefs_changed_node <- as.numeric(new_tree$tree_matrix[node_to_change, 8:ncol(new_tree$tree_matrix)])
          }

        }

      }
    } # end else split mix FALSE


    if(norm_unit_sphere){
      new_split_coefs <- new_split_coefs/sum(new_split_coefs^2)
    }

    if(any(is.na(new_split_coefs))){
      stop("NA in new_split_coefs")
    }

    if(threshold_prior != "continuous_minus_plus_1"){
      # calculate all linear combination values in the relevant node
      # lincomb_values <- unique(as.vector(X[ use_node_indices, ] %*% new_split_coefs))

      countout <- collapse::fcount(as.vector(X[ use_node_indices, ] %*% new_split_coefs), sort = TRUE)
      lincomb_values <- countout$x
    }

    if(threshold_prior == "discrete_uniform"){
      # Alternatively follow BARTMachine and choose a split value using sample on the internal values of the available
      # available_values = sort(lincomb_values, na.last = TRUE)


      # if(length(available_values) == 0){
      #   stop("length(available_values) == 0")
      # }

      # if(length(lincomb_values) == 1){
      #   new_split_value = lincomb_values[1]
      # } else if (length(lincomb_values) == 2){
      #   new_split_value = max(lincomb_values) #available_values[2]
      # }  else {
      #   # new_split_value = sample(available_values[-c(1,length(available_values))], 1)
      #   # new_split_value = resample(available_values[-c(1,length(available_values))])
      #   # new_split_value = sample(x = available_values[2:(length(available_values)-1)],size = 1)
      #   new_split_value = sample(x = lincomb_values[-c(which.min(lincomb_values), which.max(lincomb_values))],size = 1)
      #
      # }

      if(length(lincomb_values) == 1){
        # no point in splitting if only one value
        count_bad_trees = count_bad_trees + 1
        if(count_bad_trees >= max_bad_trees) {
          # print(" reached max_bad_trees = ")
          curr_tree$var_update = 0
          curr_tree$count_for_update <- rep(0, ncol(X))

          if(hyp_par_list$split_mix){
            new_tree$node_updated <- NA
            new_tree$c_ind_new <- NA
          }


          if(coef_norm_hyperprior == "varying"){
            new_tree$coef_for_sum_vec <- rep(0, ncol(X))
          }
          if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
            new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
          }
          if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
            new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
          }
          return(curr_tree)
        }else{
          next
        }

        # new_split_value = available_values[1]
      } else if (length(lincomb_values) == 2){

        if(any(countout$N < node_min_size)){
          n_bad_trees = n_bad_trees + 1
          if(n_bad_trees >= max_bad_trees) {
            # print(" reached max_bad_trees = ")
            curr_tree$var_update = 0
            curr_tree$count_for_update <- rep(0, ncol(X))

            if(hyp_par_list$split_mix){
              new_tree$node_updated <- NA
              new_tree$c_ind_new <- NA
            }


            if(coef_norm_hyperprior == "varying"){
              new_tree$coef_for_sum_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
              new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
              new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
            }
            return(curr_tree)
          }else{
            next
          }
        }

        new_split_value = lincomb_values[2]
      }  else {

        # find smallest and largest split values with less than minimum node size left and right
        runsum <- 0
        # min_ind <- 0
        for(val_ind in 1:length(lincomb_values)){
          runsum <- runsum + countout$N[val_ind]
          if(runsum >= node_min_size){
            min_ind <- val_ind +1
            break
          }
        }
        runsum <- 0
        # max_ind <- 0
        for(val_ind in length(lincomb_values):1){
          runsum <- runsum + countout$N[val_ind]
          if(runsum >= node_min_size){
            max_ind <- val_ind
            break
          }
        }

        # if((min_ind > length(lincomb_values)) | max_ind == 0  ){
        if((min_ind > length(lincomb_values)) | (min_ind > max_ind)  ){
          count_bad_trees = count_bad_trees + 1
          if(count_bad_trees >= max_bad_trees) {
            # print(" reached max_bad_trees = ")
            curr_tree$var_update = 0
            curr_tree$count_for_update <- rep(0, ncol(X))

            if(hyp_par_list$split_mix){
              new_tree$node_updated <- NA
              new_tree$c_ind_new <- NA
            }


            if(coef_norm_hyperprior == "varying"){
              new_tree$coef_for_sum_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
              new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
              new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
            }
            return(curr_tree)
          }else{
            next
          }
        }

        if(min_ind > max_ind){
          stop("min_ind > max_ind")
        }

        if(min_ind == max_ind){
          new_split_value <- lincomb_values[min_ind]
        }else{
          new_split_value <- sample(lincomb_values[min_ind:max_ind],1)
          # new_split_value = runif(1,lincomb_values[min_ind],lincomb_values[max_ind])
        }

        # split_value = sample(available_values[-c(1,length(available_values))], 1)
        # split_value = resample(available_values[-c(1,length(available_values))])
        # new_split_value = sample(available_values[-c(1)],1)
        # split_value = resample(available_values[-c(1,length(available_values))])
        # split_value = runif(1,available_values[2],available_values[length(available_values)])
      }



    }

    if(threshold_prior == "continuous_min_max"){
      # tempmin <- min(lincomb_values)
      # tempmax <- max(lincomb_values)
      # new_split_value <- runif(n = 1, min = tempmin, max = tempmax)


      if(length(lincomb_values) == 1){
        # no point in splitting if only one value
        count_bad_trees = count_bad_trees + 1
        if(count_bad_trees >= max_bad_trees) {
          # print(" reached max_bad_trees = ")
          curr_tree$var_update = 0
          curr_tree$count_for_update <- rep(0, ncol(X))

          if(hyp_par_list$split_mix){
            new_tree$node_updated <- NA
            new_tree$c_ind_new <- NA
          }


          if(coef_norm_hyperprior == "varying"){
            new_tree$coef_for_sum_vec <- rep(0, ncol(X))
          }
          if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
            new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
          }
          if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
            new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
          }
          return(curr_tree)
        }else{
          next
        }

        # new_split_value = available_values[1]
      } else if (length(lincomb_values) == 2){

        if(any(countout$N < node_min_size)){
          n_bad_trees = n_bad_trees + 1
          if(n_bad_trees >= max_bad_trees) {
            # print(" reached max_bad_trees = ")
            curr_tree$var_update = 0
            curr_tree$count_for_update <- rep(0, ncol(X))

            if(hyp_par_list$split_mix){
              new_tree$node_updated <- NA
              new_tree$c_ind_new <- NA
            }


            if(coef_norm_hyperprior == "varying"){
              new_tree$coef_for_sum_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
              new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
              new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
            }
            return(curr_tree)
          }else{
            next
          }
        }

        new_split_value = runif(1,lincomb_values[1],lincomb_values[2])
        # new_split_value = lincomb_values[2]
      }  else {

        # find smallest and largest split values with less than minimum node size left and right
        runsum <- 0
        # min_ind <- 0
        for(val_ind in 1:length(lincomb_values)){
          runsum <- runsum + countout$N[val_ind]
          if(runsum >= node_min_size){
            min_ind <- val_ind +1
            break
          }
        }
        runsum <- 0
        # max_ind <- 0
        for(val_ind in length(lincomb_values):1){
          runsum <- runsum + countout$N[val_ind]
          if(runsum >= node_min_size){
            max_ind <- val_ind
            break
          }
        }

        # if((min_ind > length(lincomb_values)) | max_ind == 0  ){
        if((min_ind > length(lincomb_values)) | (min_ind > max_ind)  ){
          count_bad_trees = count_bad_trees + 1
          if(count_bad_trees >= max_bad_trees) {
            # print(" reached max_bad_trees = ")
            curr_tree$var_update = 0
            curr_tree$count_for_update <- rep(0, ncol(X))

            if(hyp_par_list$split_mix){
              new_tree$node_updated <- NA
              new_tree$c_ind_new <- NA
            }


            if(coef_norm_hyperprior == "varying"){
              new_tree$coef_for_sum_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
              new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
            }
            if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
              new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
            }
            return(curr_tree)
          }else{
            next
          }
        }

        if(min_ind > max_ind){
          stop("min_ind > max_ind")
        }

        if(min_ind == max_ind){
          # new_split_value <- lincomb_values[min_ind]
          new_split_value <- runif(1,lincomb_values[min_ind-1],lincomb_values[min_ind])

        }else{
          # new_split_value <- sample(lincomb_values[min_ind:max_ind],1)
          # new_split_value = runif(1,lincomb_values[min_ind],lincomb_values[max_ind])
          new_split_value = runif(1,lincomb_values[min_ind-1],lincomb_values[max_ind])
        }



        # split_value = sample(available_values[-c(1,length(available_values))], 1)
        # split_value = resample(available_values[-c(1,length(available_values))])
        # new_split_value = sample(available_values[-c(1)],1)
        # split_value = resample(available_values[-c(1,length(available_values))])
        # split_value = runif(1,available_values[2],available_values[length(available_values)])
      }


    }


    if(threshold_prior == "continuous_minus_plus_1"){
      new_split_value <- runif(n = 1, min = -1, max = 1)
    }





    # Update the tree details
    # new_tree$tree_matrix[node_to_change,
    #                      c(#'split_variable',
    #                        'split_value')] = c(#new_split_variable,
    #                                            new_split_value)

    new_tree$tree_matrix[node_to_change,
                           'split_value'] = new_split_value

    new_tree$tree_matrix[node_to_change, 8:ncol(new_tree$tree_matrix)] = new_split_coefs


    # Update the tree node indices
    new_tree = fill_tree_details(new_tree, X)

    # Store the covariate name that was used in the splitting rule of the terminal node that was just changed
    # new_tree$var = c(var_changed_node, new_split_variable)
    new_tree$var_update = 1

    if(hyp_par_list$split_mix){
      new_tree$node_updated <- node_to_change
      new_tree$c_ind_new <- c_ind
      new_tree$count_for_update_new <- gamma_vec
      new_tree$count_for_update_old <- gamma_changed_node

      if(coef_norm_hyperprior == "varying"){
        new_tree$coef_for_sum_vec_new <- new_split_coefs
        new_tree$coef_for_sum_vec_old <- coefs_changed_node
      }
      if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
        new_tree$coef_for_sumsq_vec_new <- gamma_vec*(new_split_coefs )^2
        new_tree$coef_for_sumsq_vec_old <- (coefs_changed_node >= 0)*(coefs_changed_node )^2
      }
      if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
        new_tree$coef_for_logsum_vec_new <- log(abs(new_split_coefs ))
        new_tree$coef_for_logsum_vec_old <- log(abs(coefs_changed_node ))
      }

    }

    if(coef_hyperprior %in% c("univariate_normal_fixed_binomial",
                               "univariate_normal_betabinomial_onetheta",
                               "univariate_normal_betabinomial_theta_j",
                              "univariate_normal_betabinomial_theta_j_sigma_j",
                              "simplex_fixed_beta_binomial_theta_j")  ){

      new_tree$count_for_update <- gamma_vec - gamma_changed_node

    }


    if(coef_hyperprior %in% c("simplex_fixed_Dir_binomial_plusminus_theta_j",
                              "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j")){
      new_tree$count_for_update <- gamma_vec - 1*(coefs_changed_node >= 0)
    }

    if(coef_hyperprior == "simplex_fixed_Dir_trinomial_theta_j"){
      new_tree$count_min1_for_update = 1*(delta_vec < 0) - 1*(coefs_changed_node < 0)
      new_tree$count_0_for_update = 1*(delta_vec == 0) - 1*(coefs_changed_node == 0)
      new_tree$count_1_for_update = 1*(delta_vec > 0) - 1*(coefs_changed_node > 0)
    }

    if(coef_norm_hyperprior == "varying"){
      new_tree$coef_for_sum_vec <- new_split_coefs - coefs_changed_node
    }
    if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
      new_tree$coef_for_sumsq_vec <- gamma_vec*(new_split_coefs )^2 -  (coefs_changed_node >= 0)*(coefs_changed_node )^2
    }

    if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
      new_tree$coef_for_logsum_vec <- log(abs(new_split_coefs )) -  log(abs(coefs_changed_node ))
    }

    # Check for bad tree
    if(any(as.numeric(new_tree$tree_matrix[terminal_nodes, 'node_size']) <= node_min_size)) {
      count_bad_trees = count_bad_trees + 1
    } else {
      bad_trees = FALSE
    }


    if(count_bad_trees == max_bad_trees){
      # curr_tree$var = c(0, 0)
      curr_tree$var_update = 0
      curr_tree$count_for_update <- rep(0, ncol(X))

      if(hyp_par_list$split_mix){
        new_tree$node_updated <- NA
        new_tree$c_ind_new <- NA
      }


      if(coef_norm_hyperprior == "varying"){
        new_tree$coef_for_sum_vec <- rep(0, ncol(X))
      }
      if(coef_hyperprior == "univariate_normal_betabinomial_theta_j_sigma_j"){
        new_tree$coef_for_sumsq_vec <- rep(0, ncol(X))
      }
      if(coef_hyperprior == "simplex_fixed_Dir_binomial_plusminus_theta_j_xi_j"){
        new_tree$coef_for_logsum_vec <- rep(0, ncol(X))
      }
      return(curr_tree)
    }

  } # end of while loop

  # Return new_tree
  return(new_tree)

} # End of change_tree function

# swap_tree function ------------------------------------------------------

swap_tree = function(X, y, curr_tree, node_min_size) {

  stop("Swap proposals currently not supported. Code not written.")

  # Swap takes two neighbouring internal nodes and swaps around their split values and variables

  # If current tree is a stump nothing to change
  if(nrow(curr_tree$tree_matrix) == 1) return(curr_tree)

  # Create a holder for the new tree
  new_tree = curr_tree

  # Need to get the internal nodes
  internal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 0)
  terminal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 1)

  # If less than 3 internal nodes return curr_tree
  if(length(internal_nodes) < 3) return(curr_tree)

  # Find pairs of neighbouring internal nodes
  parent_of_internal = as.numeric(new_tree$tree_matrix[internal_nodes,'parent'])
  pairs_of_internal = cbind(internal_nodes, parent_of_internal)[-1,]

  # Create a while loop to get good trees
  # Create a counter to stop after a certain number of bad trees
  max_bad_trees = 2
  count_bad_trees = 0
  bad_trees = TRUE
  while(bad_trees) {
    # Re-set the tree
    new_tree = curr_tree

    # Pick a random pair
    # nodes_to_swap = sample(1:nrow(pairs_of_internal), 1)
    nodes_to_swap = sample(1:nrow(pairs_of_internal), 1)

    # Get the split variables and values for this pair
    swap_1_parts = as.numeric(new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,1],
                                                   c('split_variable', 'split_value')])
    swap_2_parts = as.numeric(new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,2],
                                                   c('split_variable', 'split_value')])

    # Update the tree details - swap them over
    new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,1],
                         c('split_variable',
                           'split_value')] = swap_2_parts
    new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,2],
                         c('split_variable',
                           'split_value')] = swap_1_parts

    # Update the tree node indices
    new_tree = fill_tree_details(new_tree, X)

    # Check for bad tree
    if(any(as.numeric(new_tree$tree_matrix[terminal_nodes, 'node_size']) <= node_min_size)) {
      count_bad_trees = count_bad_trees + 1
    } else {
      bad_trees = FALSE
    }
    if(count_bad_trees == max_bad_trees) return(curr_tree)

  } # end of while loop

  # Return new_tree
  return(new_tree)

} # End of swap_tree function
