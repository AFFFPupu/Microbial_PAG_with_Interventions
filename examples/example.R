# Define the number of nodes and the edge list
N <- 14
edge_list <- list(
  c(1, 4), c(2, 4), c(4, 5), c(5, 6), 
  c(6, 7), c(7, 5), c(7, 13), c(13, 14), 
  c(12, 14), c(10, 12), c(10, 7), c(9, 10), 
  c(3, 9), c(8, 9), c(10, 11), c(11, 8)
)

# Step 1: Simulate adjacency matrix based on the edge list
A_test <- simulate1_A_matrix(N, edge_list, 1.5)
print(A_test)

# Step 2: Convert the simulated adjacency matrix into binary adjacency form
A_test_adj <- matrix_to_adj(A_test)

# Step 3: Generate the true PAG (Partially Ancestral Graph) from the adjacency matrix
A_test_truepag <- Adj_matrix_to_PAG(A_test_adj)

# Step 4: Simulate x data for the system dynamics
# Here, 10000 data points are simulated using the system dynamics
x_test <- simulate_x_data_onematrix_julia(A_test, 0, 10000)

# Step 5: Learn the PAG structure from the simulated x data
A_test_learnedpag <- data_to_pag(x_test, 0.05, c()) # Alpha is set to 0.05

# Step 6: Plot the original adjacency matrix, true PAG, and learned PAG
par(mfrow = c(1, 3))   # Set up a 1x3 grid for plots
plot_adjmatrix(A_test_adj)               # Plot original adjacency matrix
plotAG(A_test_truepag@amat)              # Plot true PAG
plotAG(A_test_learnedpag@amat)           # Plot learned PAG from data

# Step 7: Identify cycles in the true PAG
cycles <- pag2cycles(A_test_truepag)

# Step 8: Convert cycles to a list of vectors for further processing
cycles_list <- lapply(cycles, as.vector)

# Step 9: Find the minimal set of nodes to break cycles (minimal set cover)
I_break <- find_minimal_set_cover(cycles_list)
print(I_break)

# Step 10: Intervention 1 - Break cycles at nodes specified in I1
I1 <- c(6, 10)   # Nodes to intervene
A_test_I1 <- A_test

# Update adjacency matrix by setting rows and self-loops for nodes in I1
for (i in I1) {
  A_test_I1[i, ] <- 0
  A_test_I1[i, i] <- -1
}

# Step 11: Generate the new adjacency matrix and true PAG after intervention I1
A_test_adj_I1 <- matrix_to_adj(A_test_I1)
A_test_I1_truepag <- Adj_matrix_to_PAG(A_test_adj_I1)

# Step 12: Simulate new x data after the intervention
x_test_I1 <- simulate_x_data_onematrix_julia(A_test_I1, 0, 10000)

# Step 13: Learn the PAG structure from the new x data
A_test_I1_learnedpag <- data_to_pag(x_test_I1, 0.05, c())

# Step 14: Plot the results of intervention I1
par(mfrow = c(1, 3))   # Set up a 1x3 grid for plots
plot_adjmatrix(A_test_adj_I1)            # Plot new adjacency matrix after intervention I1
plotAG(A_test_I1_truepag@amat)           # Plot true PAG after intervention I1
plotAG(A_test_I1_learnedpag@amat)        # Plot learned PAG after intervention I1

# Step 15: Identify any remaining cycles after the intervention
remaining_cycles <- pag2cycles(A_test_I1_truepag)
print(remaining_cycles)