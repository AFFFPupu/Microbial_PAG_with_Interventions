library(dplyr)
library(tidyr)
library("graph")
library("BiocGenerics")
library("pcalg")
library("igraph")
library("MixedGraphs")
library(MAGsearch)
library("ADMGs2")
#library("Rgraphviz")
library("deSolve")
library(devtools)
options(languageserver.diagnostics = FALSE)

options(JULIA_HOME = "/home/aofu/julia-1.8.5/bin/")
library(JuliaCall)
julia <- julia_setup()
de <- diffeqr::diffeq_setup()
JuliaCall::julia_library("LinearAlgebra")
JuliaCall::julia_eval("
function f_HollingII(x, theta)
    N = length(x)
    x./(ones(N,1) + theta.*x)
end
")
JuliaCall::julia_eval("
function system_dynamics!(dx, x, p, t)
    N = length(x)

    # Get (A, r, theta) from p = [A r theta]
    A = p[1:N, 1:N]
    r = p[1:N, N+1]
    theta = p[1:N, N+2]

    # If some species start with negative abundance, make them zero
    pos = findall(x .<= 0)
    x[pos] .= 0

    # System dynamics
    hollingII_term = f_HollingII(x, theta)
    dx .= (Diagonal(x) * (A * hollingII_term .+ r))
end
")


solve_system_dynamics <- function(A, r, theta, x0, tspan){
  #Example: A <- matrix(
  # c(-1.000000,  1.552315,  0.000000,  0.000000,
  # 0.000000, -1.000000,  0.000000,  0.000000,
  # 1.761856,  0.000000, -1.000000, -0.457249,
  # 0.000000,  0.000000, -1.399008, -1.000000), nrow = 4, byrow = TRUE)
  # r <- c(0.1, -0.1, 0.05, -0.05)
  # theta <- c(1.0, 1.0, 1.0, 1.0)
  # x0 <- c(1.0,1.0,1.0,1.0)
  tspan <- c(0.0, 100.0)
  p <- cbind(A, r, theta)
  JuliaCall::julia_assign("x0", x0)
  JuliaCall::julia_assign("p", p)
  JuliaCall::julia_assign("tspan", tspan)
  prob <- JuliaCall::julia_eval("ODEProblem(system_dynamics!, x0, tspan, p)")
  sol = de$solve(prob,de$Tsit5(),saveat=0.5)
  mat <- sapply(sol$u,identity)
  udf <- as.data.frame(t(mat))
  return(udf)
}



simulate_x_data_onematrix_julia <- function(A_matrix, theta_max, number_datapoints){
  data_x <- data.frame(Date=as.Date(character()),
                       File=character(), 
                       User=character(), 
                       stringsAsFactors=FALSE) 
  node_label <- c()
  for (i in 1:N){
    node_label[i] <- paste0('x',i)
  }
  for (i in 1:number_datapoints){
    xd <- runif(N)
    theta = runif(N) * theta_max
    r <- xd/(1+theta*xd)
    x0 <- runif(N)
    x<- solve_system_dynamics(A = A_matrix, r=r, theta = theta, x0=x0,tspan=tspan)
    data_x <- rbind(data_x, x[nrow(x),])
  }
  names(data_x) <- node_label
  rownames(data_x) <- NULL
  
  return(data_x)
}

simulate_x_data_Julia <- function(N, theta_max, number_datapoints, data_A_comp){
  data_x_comp <- list()
  node_label <- c()
  for (i in 1:N){
    node_label[i] <- paste0('x',i)
  }
  for (k in 1:20){
    A_matrix <- data_A_comp[[k]]
    data_x <- simulate_x_data_onematrix_julia(A_matrix,theta_max = theta_max,number_datapoints)
    rownames(data_x) <- NULL
    data_x_comp[[length(data_x_comp)+1]] <- data_x
    print(paste0(k, 'finish!!!'))
  }
  return(data_x_comp)
}

matrix_to_pag <- function(A_matrix){
  for (i in 1:N){
    for (j in 1:N){
      if (i==j){
        A_matrix[i,j] <- 0
      }
      if (i!=j){
        if (A_matrix[i,j]!=0){
          A_matrix[i,j]=1
        }
      }
    }
  }
  A_matrix <- t(A_matrix)
  A_matrix_acy <- acyclification(A_matrix)
  true_pag <- admg2pag(A_matrix_acy,latent_variables = c())
  return(true_pag)
}
pag2cycles <- function(pag){
  # pag: pcalg FCI pag
  g_cyclic <- as.undirected(graph_from_adjacency_matrix(pag@amat))
  cycles <- cliques(g_cyclic,min = 2)
  return(cycles)
}


#data simulation functions
GLV <- function(t, x, parameters){
  with(as.list(c(x, parameters)), {
    x[x < 10^-8] <- 0 # prevent numerical problems
    dxdt <- x * (r + A %*% (x/(1+theta*x)))
    list(dxdt)
  })
}

# function to plot output
plot_ODE_output <- function(out){
  out <- as.data.frame(out)
  colnames(out) <- c("time", paste("sp", 1:(ncol(out) -1), sep = "_"))
  out <- as_tibble(out) %>% gather(species, density, -time)
  pl <- ggplot(data = out) + 
    aes(x = time, y = density, colour = species) + 
    geom_line()
  show(pl)
  return(out)
}

# general function to integrate GLV
integrate_GLV <- function(r, A,theta, x0, maxtime = 100, steptime = 0.5){
  times <- seq(0, maxtime, by = steptime)
  parameters <- list(r = r, A = A, theta = theta)
  # solve numerically
  out <- ode(y = x0, times = times, 
             func = GLV, parms = parameters, 
             method = "ode45")
  # plot and make into tidy form
  # out <- plot_ODE_output(out)
  return(out)
}


# Simulation of A_matrix  (here we simulate 20 cases)
simulate_A_matrix <- function(N, C, sigma){
  # N: number of species(nodes)
  # C: connectivity
  # sigma: interaction strength
  data_A_comp <- list()
  for (k in 1:20){
    #generate A matrix
    A_matrix <- matrix(0, nrow = N, ncol = N)
    for (i in 1:N){
      for (j in 1:N){
        if (runif(1) < C){
          A_matrix[i,j] <- sigma * (ifelse(runif(1) < 0.5, runif(1, -0.8, -0.5), runif(1, 0.5, 0.8)))
        }
        if (i == j){
          A_matrix[i,j] <- -1
        } 
      }
    }
    colnames(A_matrix) <- rownames(A_matrix) <- 1:N
    data_A_comp[[length(data_A_comp)+1]] <- A_matrix
  }
  return (data_A_comp)
}

#simulate acyclic A_matrix (also 20 cases)
simulate_A_matrix_acyclic <- function(N, C, sigma){
  # N: number of species(nodes)
  # C: connectivity
  # sigma: interaction strength
  data_A_comp <- list()
  for (k in 1:20){
    #generate A matrix
    dag <- randomDAG(N,C,lB = -sigma, uB = sigma, V = as.character(1:N))
    A_matrix <- as(dag,"matrix")
    colnames(A_matrix) <- rownames(A_matrix) <- 1:N
    for (i in 1:N){
      A_matrix[i,i] <- -1
    }
    data_A_comp[[length(data_A_comp)+1]] <- A_matrix
  }
  return (data_A_comp)
}

simulate_x_data_onematrix_julia <- function(A_matrix, theta_max, number_datapoints){
  data_x <- data.frame(Date=as.Date(character()),
                       File=character(), 
                       User=character(), 
                       stringsAsFactors=FALSE) 
  node_label <- c()
  for (i in 1:N){
    node_label[i] <- paste0('x',i)
  }
  for (i in 1:number_datapoints){
    xd <- runif(N)
    theta = runif(N) * theta_max
    r <- xd/(1+theta*xd)
    x0 <- runif(N)
    x<- solve_system_dynamics(A = A_matrix, r=r, theta = theta, x0=x0,tspan=tspan)
    data_x <- rbind(data_x, x[nrow(x),])
  }
  names(data_x) <- node_label
  rownames(data_x) <- NULL
  
  return(data_x)
}

simulate_x_data_Julia <- function(N, theta_max, number_datapoints, data_A_comp){
  data_x_comp <- list()
  node_label <- c()
  for (i in 1:N){
    node_label[i] <- paste0('x',i)
  }
  for (k in 1:20){
    A_matrix <- data_A_comp[[k]]
    data_x <- simulate_x_data_onematrix_julia(A_matrix,theta_max = theta_max,number_datapoints)
    rownames(data_x) <- NULL
    data_x_comp[[length(data_x_comp)+1]] <- data_x
    print(paste0(k, 'finish!!!'))
  }
  return(data_x_comp)
}

# Simulation of x_data with corresponding A_matrices 
simulate_x_data <- function(N, theta_max, number_datapoints, data_A_comp){
  # N: number of species(nodes)
  # theta_max: intrinsic growth strength (r)
  # number_datapoints: number of steady data points generated for each A_matrix
  # data_A_comp: list containing 20 A_matrix
  data_x_comp <- list()
  node_label <- c()
  for (i in 1:N){
    node_label[i] <- paste0('x',i)
  }
  for (k in 1:20){
    A_matrix <- data_A_comp[[k]]
    
    #simulate x_data which consists of steady states
    x_data <- data.frame(Date=as.Date(character()),
                         File=character(), 
                         User=character(), 
                         stringsAsFactors=FALSE) 
    for (i in 1:number_datapoints){
      xd <- runif(N)
      theta = runif(N) * theta_max
      r <- xd/(1+theta*xd)
      x0 <- runif(N)
      x <- integrate_GLV(r, A_matrix, theta,x0)
      x_data <- rbind(x_data, x[201, 2:(N+1)])
      
      # y_data[[length(y_data)+1]] <- y[201, 2:(N+1)]
    }
    names(x_data) <- node_label
    
    data_x_comp[[length(data_x_comp)+1]] <- x_data
    print(paste0(k, 'finished'))
  }
  return(data_x_comp)
}


add_error <- function(data_x_comp, error_para, number_datapoints){
  for (k in 1:20){
    datax <- data_x_comp[[k]]
    N <- dim(datax)[2]
    for (i in 1:number_datapoints){
      for (j in 1:N){
        a <-abs(datax[i,j])
        if (a < 10e-2 || is.na(a)){
          error <-0
        } else{
          error <- runif(1,-a,a) * error_para
        }
        datax[i,j] <- datax[i,j]+error
      }
    }
    data_x_comp[[k]] <- datax
  }
  return(data_x_comp)
}

remove_zero_data <- function(data_x_comp){
  for (k in 1:20){
    datax <- data_x_comp[[k]]
    datax[datax<10e-4] <- NA
    data<-datax[complete.cases(datax),]
    data_x_comp[[k]] <- data
  }
  return(data_x_comp)
}

# represent Adjacency matrix in graph required matrix
# A(i,j) != 0 to A(i,j) = 1
# ignore self_cycles i.e. A(i,i) = 0
A_matrix_to_Adj <- function(data_A_comp){
  N <- dim(data_A_comp[[1]])[1]
  data_A_after <- list()
  for (k in 1:20){
    A_matrix <- data_A_comp[[k]]
    for (i in 1:N){
      for (j in 1:N){
        if (i==j){
          A_matrix[i,j] <- 0
        }
        if (i!=j){
          if (A_matrix[i,j]!=0){
            A_matrix[i,j]=1
          }
        }
      }
    }
    data_A_after[[length(data_A_after)+1]] <- t(A_matrix)
  }
  return (data_A_after)
}


ADMGtoMAG <- function(G){
  if (!is_ADMG(G)){
    return("input has to be an ADMG")
  }
  n <- nv(G)
  ht <- headsTails3(G,r=FALSE,max_head=2)
  
  singlevertexindex <- which(lengths(ht$heads) == 1)
  parentset <- ht$tails[singlevertexindex]
  headofsize2 <- ht$heads[which(lengths(ht$heads) == 2)]
  
  H <- mixedgraph(v = seq_len(n))
  directed <- adjMatrix(n=n)
  bidirected <- adjMatrix(n=n)
  for (i in seq_len(n)) {
    vertex <- ht$heads[[singlevertexindex[i]]]
    for (j in seq_len(length(parentset[[i]]))){
      directed[parentset[[i]][j],vertex] <- 1 
    }
  }
  for (i in seq_len(length(headofsize2))){
    bidirected[headofsize2[[i]][1],headofsize2[[i]][2]] <- bidirected[headofsize2[[i]][2],headofsize2[[i]][1]] <- 1
  }
  H$edges$bidirected <- bidirected
  H$edges$directed <- directed
  return(H)
  
}

# acyclification
# a function of acyclification of matrix
# input: adjacency matrix of a Directed Mixed Graph
# output: adjacency matrix of an Acyclic Directed Mixed Graph
acyclification <- function(Adj_matrix) {
  
  g_ <- graph_from_adjacency_matrix(Adj_matrix)
  A_matrix_removed <- Adj_matrix
  # for (i in 1:N){
  #   for (j in 1:N){
  #     if(A_matrix_removed[i,j]==1 & A_matrix_removed[j,i]==1){
  #       A_matrix_removed[i,j]<- A_matrix_removed[j,i]<-0
  #     }
  #   }
  # }
  g_removed <- graph_from_adjacency_matrix(A_matrix_removed)
  
  # compute strong connected components of g_removed
  components <- components(g_removed, mode= "strong")
  g_membership <- components[[1]]
  cluster_sizes <- components[[2]]
  g_cluster <- c()
  for (i in 1:length(cluster_sizes)){
    if (cluster_sizes[i] > 1){
      g_cluster[length(g_cluster)+1] <- i
    }
  }
  
  A_matrix_acy <- Adj_matrix
  if (!(is.null(g_cluster))){
    
    g_cluster_nodes <- list()
    for (i in 1:length(g_cluster)){
      cluster_no <- g_cluster[i]
      cluster_no_nodes <- c()
      for (j in 1:length(g_membership)){
        if (g_membership[[j]] == cluster_no){
          cluster_no_nodes[length(cluster_no_nodes)+1] <- j
        }
      }
      g_cluster_nodes[[length(g_cluster_nodes)+1]] <- cluster_no_nodes
    }
    for (i in 1:length(g_cluster_nodes)){
      # we use one specific acyclification that we replace all strongly connected components of G by fully connected bidirected components without any direct edges.
      connected_nodes <- g_cluster_nodes[[i]]
      for (j in connected_nodes){
        for (k in connected_nodes){
          if (j != k){
            A_matrix_acy[j,k] <- 1
            for (t in 1:N){
              if (t != j & t != k){
                if (A_matrix_acy[t,j] == 1){
                  A_matrix_acy[t,k] <- 1
                }
              }
              
            }
          }
        }
      }
    }
  }
  return(A_matrix_acy)
}

#encode ADMG as  adjcency matrix of PAG-type
# amat[i,j] = 0 iff no edge btw i,j
# amat[i,j] = 1 iff i *-o j
# amat[i,j] = 2 iff i *-> j
# amat[i,j] = 3 iff i *-- j

# admg to pag
# input: admg adjcency matrix (in our case, since no bidirected edges in admg, entries with 1 on both i,j and j,i )
admg2pag <- function(admg_matrix,latent_variables){
  if (is.null(latent_variables)){
    edges = c()
    N = dim(admg_matrix)[1]
    amat = admg_matrix
    for (i in 1:N){
      for (j in i:N ){
        if (amat[i,j] ==1 & amat[j,i]==1){
          edges[length(edges)+1] <- paste0(i,"<->",j)
        }
        if (amat[i,j] == 1 & amat[j,i] == 0){
          edges[length(edges)+1] <- paste0(i,"->",j)
        }
        if (amat[j,i] == 1 & amat[i,j] == 0){
          edges[length(edges)+1] <- paste0(j,"->",i)
        }
      }
      edges[length(edges)+1] <- toString(i)
    }
    graph_admg <- graphCr(edges)
    G <- ADMGtoMAG(graph_admg)
    Gpag <- ConstructPAG(G)
    H <- convert(Gpag,format = "PAG") 
    amat <- matrix(H,N,N)
    labels = as.character(1:N)
    rownames(amat) <- labels
    colnames(amat) <- labels
    fci.pag <- new("fciAlgo",amat = amat)
    return(fci.pag)
  } else{
    #if (is.null(latent_variables)){
    #   A_matrix_acy_admg <- admg_matrix
    #   N = dim(admg_matrix)[1]
    #   
    #   for (i in 1:N){
    #     for (j in 1:N){
    #       if (i!=j & A_matrix_acy_admg[i,j]==1 & A_matrix_acy_admg[j,i]==1){
    #         A_matrix_acy_admg[i,j] <- A_matrix_acy_admg[j,i] <- 2
    #       }
    #       if (i!=j & A_matrix_acy_admg[i,j]==1 & A_matrix_acy_admg[j,i]!=1){
    #         A_matrix_acy_admg[i,j] <- 2
    #         A_matrix_acy_admg[j,i] <- 3 
    #       }
    #     }
    #   }
    #   
    #   indepTest <- dsepAMTest
    #   suffStat<-list(g=A_matrix_acy_admg,verbose=FALSE)
    #   labels = as.character(1:N)
    #   (fci.pag <- fci(suffStat,indepTest,alpha = 0.5,verbose=TRUE,selectionBias=FALSE, labels = labels))
    #   return (fci.pag)
    # } else{
    N_latent = length(latent_variables)
    edges = c()
    N = dim(admg_matrix)[1]
    amat = admg_matrix
    for (i in 1:N){
      for (j in i:N ){
        if (amat[i,j] ==1 & amat[j,i]==1){
          edges[length(edges)+1] <- paste0(i,"<->",j)
        }
        if (amat[i,j] == 1 & amat[j,i] == 0){
          edges[length(edges)+1] <- paste0(i,"->",j)
        }
        if (amat[j,i] == 1 & amat[i,j] == 0){
          edges[length(edges)+1] <- paste0(j,"->",i)
        }
      }
      edges[length(edges)+1] <- toString(i)
    }
    graph_admg <- graphCr(edges)
    graph_latent_proj <- latentProject(graph_admg, latent = latent_variables, only_directed = FALSE, sort = 1)
    graph_latent_proj$v <- graph_admg$v
    G <- ADMGtoMAG(graph_latent_proj)
    G1 <- latentProject(G, latent = latent_variables, only_directed = FALSE, sort = 1)
    Gpag <- ConstructPAG(G)
    H <- convert(Gpag,format = "PAG") 
    amat <- matrix(H,N,N)
    L <- latent_variables
    labels = as.character(1:N)
    amat <- amat[-L,-L]
    rownames(amat) <- labels[-L]
    colnames(amat) <- labels[-L]
    # indepTest <- dsepAMTest
    # suffStat<-list(g=amat,verbose=FALSE)
    # (fci.pag <- fci(suffStat,indepTest,alpha = 0.5,labels =labels,verbose=FALSE))
    fci.pag <- new("fciAlgo",amat = amat)
    return(fci.pag)
    
  }
  
}
simulate_latent_variables_comp <- function(N,p){
  latent_variables_comp <-list()
  if (p==0){
    latent_variables_comp <- vector(mode = "list", length = 20)
    return (latent_variables_comp)
  } else{
    for (i in 1:20){
      N_ <- floor(N*p)
      latent_variables <- sample(1:N,N_,replace = F)
      latent_variables_comp[[length(latent_variables_comp)+1]] <- latent_variables
    }
    return (latent_variables_comp)
  }
}


# from ground truth directed graph to ground truth pag
A_comp_to_pag_comp <- function(data_A_comp,latent_variables_comp){
  true_pag_comp <- list()
  for (i in 1:20){
    A_matrix <- data_A_comp[[i]]
    A_matrix_acy <- acyclification(A_matrix)
    true_pag <- admg2pag(A_matrix_acy,latent_variables = latent_variables_comp[[i]])
    true_pag_comp[[length(true_pag_comp)+1]] <- true_pag
    print(paste0(i,'done'))
  }
  return (true_pag_comp)
}

# learn pag from data
# input: data_x[[i]], significance level
# output: fci.pag object
data_to_pag <- function(data_x, alpha, latent_variables){
  N <- dim(data_x)[2]
  if(length(latent_variables) == 0){
    suffstatdata_x <- list(C = cor(data_x), n=nrow(data_x))
    labels = colnames(data_x)
    fci.data_x <- fci(suffstatdata_x, indepTest = gaussCItest,
                      alpha = alpha,labels = labels)
    return (fci.data_x)
  } else{
    l<- latent_variables
    data_x <- data_x[-l]
    suffstatdata_x <- list(C = cor(data_x), n=nrow(data_x))
    labels = colnames(data_x)
    fci.data_x <- fci(suffstatdata_x, indepTest = gaussCItest,
                      alpha = alpha,labels = labels)
    return (fci.data_x)
  }
  
  
}


# learn pags from data_x_comp
# input: data_x_comp
# output: pag_comp
data_to_pag_comp <- function(data_x_comp,latent_variables_comp, alpha){
  pag_comp <- list()
  for (i in 1:20){
    data_x <- data_x_comp[[i]]
    l <- latent_variables_comp[[i]]
    fci.pag <- data_to_pag(data_x,alpha=alpha,latent_variables = l)
    pag_comp[[length(pag_comp)+1]] <- fci.pag
  }
  return (pag_comp)
}

ADMGtoMAG<-function(G){
  if (!is_ADMG(G)){
    return("input has to be an ADMG")
  }
  n <- nv(G)
  ht <- headsTails3(G,r=FALSE,max_head=2)
  
  singlevertexindex <- which(lengths(ht$heads) == 1)
  parentset <- ht$tails[singlevertexindex]
  headofsize2 <- ht$heads[which(lengths(ht$heads) == 2)]
  
  H <- mixedgraph(v = seq_len(n))
  directed <- adjMatrix(n=n)
  bidirected <- adjMatrix(n=n)
  for (i in seq_len(n)) {
    vertex <- ht$heads[[singlevertexindex[i]]]
    for (j in seq_len(length(parentset[[i]]))){
      directed[parentset[[i]][j],vertex] <- 1 
    }
  }
  for (i in seq_len(length(headofsize2))){
    bidirected[headofsize2[[i]][1],headofsize2[[i]][2]] <- bidirected[headofsize2[[i]][2],headofsize2[[i]][1]] <- 1
  }
  H$edges$bidirected <- bidirected
  H$edges$directed <- directed
  return(H)
  
}



precision_pag <- function(true_pag, learned_pag){
  true_pag_adj <- true_pag@amat
  learned_pag_adj <- learned_pag@amat
  N <- dim(learned_pag_adj)[1]
  #compute the number of edges in learned PAG
  number_of_edges <- 0
  for (i in 1:N){
    for (j in 1:N){
      if (learned_pag_adj[i,j] !=0 & learned_pag_adj[j,i] !=0){
        number_of_edges <- number_of_edges + 1
      }
    }
  }
  
  #compute the number of edges with corrent orientations
  number_of_correct_edges <- 0
  for (i in 1:N){
    for (j in 1:N){
      if (learned_pag_adj[i,j] == true_pag_adj[i,j] &
          learned_pag_adj[j,i] == true_pag_adj[j,i] &
          learned_pag_adj[i,j] != 0 &
          learned_pag_adj[j,i] != 0){
        number_of_correct_edges <- number_of_correct_edges +1 
      }
    }
  }
  if (number_of_correct_edges == 0 & number_of_edges == 0){
    precision <- 1
  }
  else {
    precision <- number_of_correct_edges/number_of_edges
  }
  
  return (precision)
}

recall_pag <- function(true_pag, learned_pag){
  true_pag_adj <- true_pag@amat
  learned_pag_adj <- learned_pag@amat
  N <- dim(learned_pag_adj)[1]
  #compute the number of edges in ground truth PAG
  number_of_edges <- 0
  for (i in 1:N){
    for (j in 1:N){
      if (true_pag_adj[i,j] !=0 & true_pag_adj[j,i] !=0){
        number_of_edges <- number_of_edges + 1
      }
    }
  }
  
  #compute the number of edges with corrent orientations
  number_of_correct_edges <- 0
  for (i in 1:N){
    for (j in 1:N){
      if (learned_pag_adj[i,j] == true_pag_adj[i,j] &
          learned_pag_adj[j,i] == true_pag_adj[j,i] &
          learned_pag_adj[i,j] != 0 &
          learned_pag_adj[j,i] != 0){
        number_of_correct_edges <- number_of_correct_edges +1 
      }
    }
  }
  if (number_of_correct_edges == 0 & number_of_edges == 0){
    recall <- 1
  }
  else {
    recall <- number_of_correct_edges/number_of_edges
  }
  return (recall)
}

adjrecall_pag <- function(true_pag, learned_pag){
  true_pag_adj <- true_pag@amat
  learned_pag_adj <- learned_pag@amat
  N <- dim(learned_pag_adj)[1]
  #compute the number of edges in ground truth PAG
  number_of_edges <- 0
  for (i in 1:N){
    for (j in 1:N){
      if (true_pag_adj[i,j] !=0 & true_pag_adj[j,i] !=0){
        number_of_edges <- number_of_edges + 1
      }
    }
  }
  
  #compute the number of edges with corrent orientations
  number_of_correct_edges <- 0
  for (i in 1:N){
    for (j in 1:N){
      if (true_pag_adj[i,j] !=0 &
          true_pag_adj[j,i] != 0 &
          learned_pag_adj[i,j] != 0 &
          learned_pag_adj[j,i] != 0){
        number_of_correct_edges <- number_of_correct_edges +1 
      }
    }
  }
  if (number_of_correct_edges == 0 & number_of_edges == 0){
    adjrecall <- 1
  }
  else {
    adjrecall <- number_of_correct_edges/number_of_edges
  }
  return (adjrecall)
}

adjprecision_pag <- function(true_pag, learned_pag){
  true_pag_adj <- true_pag@amat
  learned_pag_adj <- learned_pag@amat
  N <- dim(learned_pag_adj)[1]
  #compute the number of edges in learned PAG
  number_of_edges <- 0
  for (i in 1:N){
    for (j in 1:N){
      if (learned_pag_adj[i,j] !=0 & learned_pag_adj[j,i] !=0){
        number_of_edges <- number_of_edges + 1
      }
    }
  }
  
  #compute the number of edges with corrent orientations
  number_of_correct_edges <- 0
  for (i in 1:N){
    for (j in 1:N){
      if (true_pag_adj[i,j] !=0 &
          true_pag_adj[j,i] != 0 &
          learned_pag_adj[i,j] != 0 &
          learned_pag_adj[j,i] != 0){
        number_of_correct_edges <- number_of_correct_edges +1 
      }
    }
  }
  if (number_of_correct_edges == 0 & number_of_edges == 0){
    adjrecall <- 1
  }
  else {
    adjrecall <- number_of_correct_edges/number_of_edges
  }
  return (adjrecall)
}

compute_metrics_comp <- function(N, C, sigma, theta_max, number_datapoints, alpha){
  data_A_comp <- simulate_A_matrix(N, C, sigma)
  data_x_comp <- simulate_x_data(N, theta_max, number_datapoints, data_A_comp = data_A_comp)
  data_A_comp <- A_matrix_to_Adj(data_A_comp = data_A_comp)
  true_pag_comp <- A_comp_to_pag_comp(data_A_comp = data_A_comp)
  learned_pag_comp <- data_to_pag_comp(data_x_comp = data_x_comp, alpha)
  precision_comp <- c()
  for (i in 1:200){
    precision_comp[i] <- precision_pag(true_pag = true_pag_comp[[i]], learned_pag = learned_pag_comp[[i]])
  }
  recall_comp <- c()
  for (i in 1:200){
    recall_comp[i] <- recall_pag(true_pag = true_pag_comp[[i]], learned_pag = learned_pag_comp[[i]])
  }
  adjrecall_comp <- c()
  for (i in 1:200){
    adjrecall_comp[i] <- adjrecall_pag(true_pag = true_pag_comp[[i]], learned_pag = learned_pag_comp[[i]])
  }
  return (list(precision_comp,recall_comp,adjrecall_comp))
}


pag2cycles <- function(pag){
  # pag: pcalg FCI pag
  g_cyclic <- as.undirected(graph_from_adjacency_matrix(pag@amat))
  #cycles <- largest_cliques(g_cyclic)
  cycles <- max_cliques(min =3, g_cyclic)
  return(cycles)
}

plot_adjmatrix <- function(adjmatrix){
  N<-dim(adjmatrix)[1]
  graph <- new("graphNEL", nodes=as.character(1:N), edgemode="directed")
  for (i in 1:N){
    for (j in 1:N){
      if (adjmatrix[i,j] == 1){
        graph <- addEdge(as.character(i),as.character(j), graph)
      }
    }
  }
  plot(graph)
}


find_minimal_set_cover <- function(sets) {
  all_elements <- unique(unlist(sets))
  covered <- vector("list", length = 0)
  target_set <- c()
  
  while(length(covered) < length(all_elements)) {
    # Find the set with the most uncovered elements
    set_sizes <- sapply(sets, function(s) length(setdiff(s, covered)))
    max_index <- which.max(set_sizes)
    
    # Pick an element from this set and add to the target set
    element <- setdiff(sets[[max_index]], covered)[1]
    target_set <- c(target_set, element)
    
    # Update covered elements
    covered <- unique(c(covered, sets[[max_index]]))
  }
  
  return(target_set)
}

simulate1_A_matrix <- function(N,edge_list, sigma){
  #generate A matrix
  A_matrix <- matrix(0, nrow = N, ncol = N)
  for (i in 1:length(edge_list)){
    A_matrix[edge_list[[i]][2],edge_list[[i]][1]] <- sigma * (ifelse(runif(1) < 0.5, runif(1, -0.8, -0.5), runif(1, 0.5, 0.8)))
  }
  for (i in 1:N){
    A_matrix[i,i] <- -1
  }
  colnames(A_matrix) <- rownames(A_matrix) <- 1:N
  return (A_matrix)
}

simulate_x_data_onematrix_julia <- function(A_matrix, theta_max, number_datapoints){
  data_x <- data.frame(Date=as.Date(character()),
                       File=character(), 
                       User=character(), 
                       stringsAsFactors=FALSE) 
  node_label <- c()
  for (i in 1:N){
    node_label[i] <- paste0('x',i)
  }
  for (i in 1:number_datapoints){
    xd <- runif(N)
    theta = runif(N) * theta_max
    r <- xd/(1+theta*xd)
    x0 <- runif(N)
    x<- solve_system_dynamics(A = A_matrix, r=r, theta = theta, x0=x0,tspan=tspan)
    data_x <- rbind(data_x, x[nrow(x),])
  }
  names(data_x) <- node_label
  rownames(data_x) <- NULL
  
  return(data_x)
}

matrix_to_adj <- function(A_matrix){
  for (i in 1:N){
    for (j in 1:N){
      if (i==j){
        A_matrix[i,j] <- 0
      }
      if (i!=j){
        if (A_matrix[i,j]!=0){
          A_matrix[i,j]=1
        }
      }
    }
  }
  return(t(A_matrix))
}

Adj_matrix_to_PAG <- function(A_matrix){
  A_matrix_acy <- acyclification(A_matrix)
  true_pag <- admg2pag(A_matrix_acy,latent_variables = c())
  return (true_pag)
}





