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

