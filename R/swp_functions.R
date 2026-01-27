# Core functions for computing #
# Weighted SWP (Small-World Propensity) #


# Assign weights to adjacency 
# matrix based on distance
assign_weights = function(adj, epsilon = 1) {
  N = nrow(adj)
  weights = matrix(0, N, N)
  distances = matrix(0, N, N)
  
  for (i in 1:N) {
    for (j in 1:N) {
      dist_ij = min(abs(i - j), N - abs(i - j))
      distances[i, j] = dist_ij
    }
  }
  
  Dmax = max(distances) + epsilon
  
  weights[adj == 1] = Dmax - distances[adj == 1]
  return(weights)
}


# Generate a Weighted Small-World
# Network(Lattice network)
generate_weighted_sw = function(N, r) {
  
  # Initialize adjacency matrix
  adj_matrix = matrix(0, nrow = N, ncol = N)  
  
  # Assign positions to nodes
  node_positions = 1:N  
  
  # Initialize a ring lattice
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      
      # Compute the shortest distance
      dist = min(abs(node_positions[i] - node_positions[j]),
                 N - abs(node_positions[i] - node_positions[j]))
      
      # If the distance is less than or equal to r
      if (dist <= r) {
        adj_matrix[i, j] = 1  
        adj_matrix[j, i] = 1 
      }
    }
  }
  
  # Assign weights based on distance
  weight_matrix = assign_weights(adj_matrix)
  
  return(weight_matrix)
}


# Function to re-wire the lattice
# and generate random graph
rewire_graph = function(W_latt, p) {
  n = nrow(W_latt)
  W_new = W_latt
  
  # Get all upper-triangular edge indices
  edges = which(upper.tri(W_latt), arr.ind = TRUE)
  weights = W_latt[edges]
  
  # Decide which edges to rewire (based on probability p)
  n_edges = nrow(edges)
  rewire_mask = runif(n_edges) < p
  
  # Indices of edges to rewire and keep
  rewire_edges = edges[rewire_mask, , drop = FALSE]
  keep_edges = edges[!rewire_mask, , drop = FALSE]
  
  # Clear upper triangle
  W_new[edges] = 0
  
  # Set preserved weights
  W_new[keep_edges] = weights[!rewire_mask]
  
  # Set rewired weights if any
  if (nrow(rewire_edges) > 0) {
    shuffled_weights = sample(weights[rewire_mask])
    W_new[rewire_edges] = shuffled_weights
  }
  
  # We need symmetric matrix
  W_new[lower.tri(W_new)] = t(W_new)[lower.tri(W_new)]
  
  return(W_new)
}



# Function to calculate weighted
# clustering coefficient
weighted_clustering_coeff = function(W) {
  N = nrow(W) 
  
  # Initialize clustering coefficients
  Cw = numeric(N) 
  
  # Loop through each node to 
  # compute clustering coefficient
  for (i in 1:N) {
    
    # get neighbors and its degree
    neighbors = which(W[i, ] > 0) 
    k = length(neighbors)  
    
    if (k < 2) {
      
      # If there are less than two neighbors,
      # clustering coefficient is 0
      Cw[i] = 0  
    } else {
      sum_weights = 0
      
      # Find the maximum weight in the network
      max_weight = max(W)  
      
      for (j in 1:(k-1)) {
        for (l in (j+1):k) {
          if (W[neighbors[j], neighbors[l]] > 0) {
            sum_weights = sum_weights + (W[i, neighbors[j]] * W[neighbors[j], neighbors[l]] * W[neighbors[l], i])^(1/3)
          }
        }
      }
      
      Cw[i] = sum_weights / (k * (k - 1)* max_weight)  
    }
  }
  return(mean(Cw))  
}

# Function to calculate weighted 
# characteristic path length
weighted_path_length = function(W) {
  
  N = nrow(W) 
  
  # Invert weights (avoid division by zero)
  W_inv = 1 / (W + (W == 0))  
  
  dist_matrix = matrix(Inf, nrow = N, ncol = N)  
  diag(dist_matrix) = 0  
  
  # Fill distance matrix with inverse weights
  for (i in 1:N) {
    for (j in 1:N) {
      if (W[i, j] > 0) {
        dist_matrix[i, j] = W_inv[i, j] 
      }
    }
  }
  
  # Floyd-Warshall algorithm for 
  # all-pairs shortest paths
  for (k in 1:N) {
    for (i in 1:N) {
      for (j in 1:N) {
        dist_matrix[i, j] = min(dist_matrix[i, j], dist_matrix[i, k] + dist_matrix[k, j])
      }
    }
  }
  # Compute characteristic path 
  # length using formula from paper
  
  L = sum(dist_matrix) / (N * (N - 1))
  
  return(L)
}


# Function to compute average 
# degree from weighted matrix
compute_avg_degree = function(W_obs) {
  
  # Convert to binary adjacency matrix
  A = ifelse(W_obs > 0, 1, 0)
  
  # Remove self-connections (diagonal)
  diag(A) = 0
  
  # Compute degree for each node
  degrees = rowSums(A)
  
  # Compute average degree
  avg_degree = mean(degrees)
  
  return(avg_degree)
}


# Compute appropriate r for 
# lattice model
compute_radius_from_W = function(W_obs) {
  avg_deg = compute_avg_degree(W_obs)
  
  # In a ring lattice, each node 
  # connects to 2*r neighbors
  r = round(avg_deg / 2)
  return(r)
}


# Function to estimate Weighted SWP #
weighted_SWP = function(W_hat, n_iter , p = 1) {
  
  # Compute radius from observed matrix
  r = compute_radius_from_W(W_hat)
  # Number of nodes
  N = nrow(W_hat)
  
  # Clustering coefficient and path length for observed weighted network
  C = weighted_clustering_coeff(W_hat)
  L = weighted_path_length(W_hat)
  
  # Generate weighted lattice graph
  W_latt = generate_weighted_sw(N, r)
  # Normalize lattice weights to match sum of observed matrix weights
  W_latt = W_latt * (sum(W_hat) / sum(W_latt))
  
  # Clustering coefficient and path length for lattice network
  C_latt = weighted_clustering_coeff(W_latt)
  L_latt = weighted_path_length(W_latt)
  
  # Initialize vectors for random graph stats
  C_vals = numeric(n_iter)
  L_vals = numeric(n_iter)
  
  for (i in 1:n_iter) {
    W_rand = rewire_graph(W_latt, p = p)
    C_vals[i] = weighted_clustering_coeff(W_rand)
    L_vals[i] = weighted_path_length(W_rand)
  }
  
  # Mean of finite values only
  C_rand = mean(C_vals[is.finite(C_vals)])
  L_rand = mean(L_vals[is.finite(L_vals)])
  
  
  # Normalize clustering coefficient and average path length
  C_norm = (C_latt - C) / (C_latt - C_rand)
  L_norm = (L - L_rand) / (L_latt - L_rand)
  
  # Bound normalized values to [0,1]
  C_norm = min(max(C_norm, 0), 1)
  L_norm = min(max(L_norm, 0), 1)
  
  # Compute Small World Propensity
  theta_W = 1 - sqrt((C_norm^2 + L_norm^2) / 2)
  
  return(theta_W)
}



