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
