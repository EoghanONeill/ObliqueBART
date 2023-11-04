# -------------------------------------------------------------------------#
# Description: this script contains auxiliar functions needed to update    #
# the trees with details and to map the predicted values to each obs       #
# -------------------------------------------------------------------------#

# 1. fill_tree_details: takes a tree matrix and returns the number of obs in each node in it and the indices of each observation in each terminal node
# 2. get_predictions: gets the predicted values from a current set of trees
# 3. get_children: it's a function that takes a node and, if the node is terminal, returns the node. If not, returns the children and calls the function again on the children
# 4. resample: an auxiliar function
# 5. get_ancestors: get the ancestors of all terminal nodes in a tree
# 6. update_s: full conditional of the vector of splitting probability.
# 7. get_number_distinct_cov: given a tree, it returns the number of distinct covariates used to create its structure

# Fill_tree_details -------------------------------------------------------

fill_tree_details = function(curr_tree, X) {

  # Collect right bits of tree
  tree_matrix = curr_tree$tree_matrix

  # Create a new tree matrix to overwrite
  new_tree_matrix = tree_matrix

  # Start with dummy node indices
  node_indices = rep(1, nrow(X))

  # For all but the top row, find the number of observations falling into each one
  for(i in 2:nrow(tree_matrix)) {

    # Get the parent
    curr_parent = as.numeric(tree_matrix[i,'parent'])

    # Find the split variable and value of the parent
    # split_var = as.numeric(tree_matrix[curr_parent,'split_variable'])
    split_coefs = as.numeric(tree_matrix[curr_parent,8:(ncol(tree_matrix))])
    split_val = as.numeric(tree_matrix[curr_parent, 'split_value'])

    # print("split_val= ")
    # print(split_val)

    if(is.na(split_val)){
      stop("split_val == NA")
    }
    if(any(is.na(split_coefs))){
      stop("split_coefs contains NA")
    }

    # Find whether it's a left or right terminal node
    left_or_right = ifelse(tree_matrix[curr_parent,'child_left'] == i,
                           'left', 'right')

    if(left_or_right == 'left') {
      # If left use less than condition
      # new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] < split_val)
      # node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] < split_val] = i

      # print("ncol(X[node_indices == curr_parent,]")
      # print(ncol(X[node_indices == curr_parent,]))
      # print("split_coefs = " )
      # print(split_coefs)

      new_tree_matrix[i,'node_size'] <- sum( ( X[node_indices == curr_parent,] %*% split_coefs )   < split_val )

      # print("new_tree_matrix[i,'node_size'] = ")
      # print(new_tree_matrix[i,'node_size'])

      node_indices[node_indices == curr_parent][ (X[node_indices == curr_parent,] %*% split_coefs ) < split_val ] <-  i



    } else {
      # If right use greater than condition
      # new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] >= split_val)
      # node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] >= split_val] = i


      # # replaced >= with ! < because values equal or very close to equal to split vlaue not allocated to a node.
      # temp_vals <- X[node_indices == curr_parent,] %*% split_coefs
      # new_tree_matrix[i,'node_size'] <- sum(  !(( temp_vals ) < split_val)  |(temp_vals == split_val) )
      #
      # # print("new_tree_matrix[i,'node_size'] = ")
      # # print(new_tree_matrix[i,'node_size'])
      #
      # node_indices[node_indices == curr_parent][  !( ( temp_vals ) <  split_val ) |(temp_vals == split_val) ] <- i

      # right node always after left, so if not filled in as the left node, remaining indices are all the right node
      new_tree_matrix[i,'node_size'] <- sum(node_indices == curr_parent)
      node_indices[node_indices == curr_parent] <- i



    }
  } # End of loop through table

  # print("new_tree_matrix = ")
  # print(new_tree_matrix)
  # print("new_tree_matrix$terminal = ")
  # print(new_tree_matrix$terminal)

  # if(any( sort(unique(node_indices)) != which(new_tree_matrix[,"terminal"] == 1))){
  #
  #   indnoteq <- which( sort(unique(node_indices)) != which(new_tree_matrix[,"terminal"] == 1))
  #   missing_inds <- which(node_indices == indnoteq)
  #
  #
  #   print("node_indices = ")
  #   print(node_indices)
  #
  #   print("curr_parent = ")
  #   print(curr_parent)
  #
  #   print("new_tree_matrix = ")
  #   print(new_tree_matrix)
  #   print("  sort(unique(node_indices)) = ")
  #   print(  sort(unique(node_indices)) )
  #
  #   print("X[node_indices == curr_parent,] %*% split_coefs    = ")
  #   print(X[node_indices == curr_parent,] %*% split_coefs )
  #
  #   print("split_val = ")
  #   print(split_val)
  #
  #   print("!( X[node_indices == curr_parent,] %*% split_coefs   < split_val )  ")
  #   print(!(X[node_indices == curr_parent,] %*% split_coefs   < split_val ))
  #
  #   print("X[node_indices == curr_parent,] %*% split_coefs   < split_val  ")
  #   print(X[node_indices == curr_parent,] %*% split_coefs   < split_val )
  #
  #
  #   print("X[node_indices == curr_parent,] %*% split_coefs   - split_val  ")
  #   print(X[node_indices == curr_parent,] %*% split_coefs   - split_val )
  #
  #   # print("X[missing_inds,] %*% split_coefs )   = ")
  #   # print(X[missing_inds,] %*% split_coefs )
  #
  #
  #
  #   stop("any( sort(unique(node_indices)) != which(new_tree_matrix[,1] ==1))")
  # }


  return(list(tree_matrix = new_tree_matrix,
              node_indices = node_indices))

} # End of function

# Get predictions ---------------------------------------------------------

get_predictions = function(trees, X, single_tree = FALSE) {

  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]

  # Normally trees will be a list of lists but just in case
  if(single_tree) {
    # Deal with just a single tree
    if(nrow(trees$tree_matrix) == 1) {
      predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      unique_node_indices = unique(trees$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(trees, X)$node_indices
      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {
        predictions[curr_X_node_indices == unique_node_indices[i]] =
          trees$tree_matrix[unique_node_indices[i], 'mu']
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees
  } else {
    # Do a recursive call to the function
    partial_trees = trees
    partial_trees[[1]] = NULL # Blank out that element of the list
    predictions = get_predictions(trees[[1]], X, single_tree = TRUE)  +
      get_predictions(partial_trees, X,
                      single_tree = length(partial_trees) == 1)
    #single_tree = !is.null(names(partial_trees)))
    # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)
  }

  return(predictions)
}

# get_children ------------------------------------------------------------

get_children = function(tree_mat, parent) {
  # Create a holder for the children
  all_children = NULL
  if(as.numeric(tree_mat[parent,'terminal']) == 1) {
    # If the node is terminal return the list so far
    return(c(all_children, parent))
  } else {
    # If not get the current children
    curr_child_left = as.numeric(tree_mat[parent, 'child_left'])
    curr_child_right = as.numeric(tree_mat[parent, 'child_right'])
    # Return the children and also the children of the children recursively
    return(c(all_children,
             get_children(tree_mat,curr_child_left),
             get_children(tree_mat,curr_child_right)))
  }
}

# Sample function ----------------------------------------------------------

resample <- function(x, ...) x[sample.int(length(x), size=1), ...]

update_s = function(var_count, p, alpha_s){
  s_ = rdirichlet(1, alpha_s/p + var_count)
  return(s_)
}

# get_number_distinct_cov <- function(tree){
#
#   # Select the rows that correspond to internal nodes
#   which_terminal = which(tree$tree_matrix[,'terminal'] == 0)
#   # Get the covariates that are used to define the splitting rules
#   num_distinct_cov = length(unique(tree$tree_matrix[which_terminal,'split_variable']))
#
#   return(num_distinct_cov)
# }

sample_move = function(curr_tree, i, nburn){

  if (nrow(curr_tree$tree_matrix) == 1 || i < max(floor(0.1*nburn), 10)) {
    type = 'grow'
  } else {
    type = sample(c('grow', 'prune', 'change'), 1)
  }
  return(type)
}
