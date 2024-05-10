library(dplyr)

algo_2DDR <- function(score, xi, y, 
                      cost_matrix,
                      k = 10,
                      pp_max = 1,
                      represent = F, 
                      regular = F){
  
  index <- (0:k+1)
  
  ### --- Grid definition --- ###
  d1 <- as.vector(quantile(score, p = seq(1, 0, -1/k)))
  d2 <- as.vector(quantile(xi, p = seq(1, 0, -1/k)))
  
  if(regular){
    d1 <- seq( max(score), min(score), by = -(max(score) - min(score))/k)
    d2 <- seq( max(xi), min(xi), by = -(max(xi) - min(xi))/k)
  }

  ### --- Loss function definition --- ###
  loss_function_grid <- function(score, xi, cost_matrix, region){
    y_hat <- rep(0, length(score))
    for (j in 1:nrow(region)){
      y_hat <- ifelse(score >= d1[region[j,1]] & xi >= d2[region[j,2]], 1, y_hat)
    }
    return( list(lf = sum(cost_matrix[y_hat == 0, 1]) + sum(cost_matrix[y_hat == 1, 2]),
                 pp = sum(y_hat)/length(score)
    )
    )
  }
  
  ### --- Starting values --- ###
  region <- data.frame(d1 = c(1), d2 = c(1))
  FD_prev <- region
  lf <- sum(cost_matrix[,1])
  t0 <- Sys.time()
  iter <- 1

  ##### --- Iterative process --- #####
  while(TRUE){
    
    ### --- New frontier points --- ###
    FD <- c()
    for (i in 1:nrow(FD_prev)){
      FD <- rbind(FD, 
                      c(FD_prev[i,1] + 1, FD_prev[i, 2]), 
                      c(FD_prev[i,1], FD_prev[i, 2] + 1)
      )
    }
    
    # --- Select only points in the decision space
    FD <- data.frame( unique(FD[FD[,1] <= k + 1 | FD[,2] <= k + 1, ]) )
    FD[FD > k + 1] <- k + 1
    
    # --- Omit points inside the previous DF and region
    for (i in 1:nrow(FD)){
      for (j in 1:nrow(FD)){
        if(FD[i,1] < FD[j,1] & FD[i,2] < FD[j,2]){
          FD[i,] <- c(0,0)
        }
      }
    }
    FD <- FD[FD[,1] > 0,]
    
    for (i in 1:nrow(FD)){
      for (j in 1:nrow(region)){
        if(FD[i,1] <= region[j,1] & FD[i,2] <= region[j,2]){
          FD[i,] <- c(0,0)
        }
      }
    }
    FD <- FD[FD[,1] > 0,]
    
    # --- Formatting
    colnames(FD) <- c("d1", "d2")

    # --- Check limit in the decision space is not reached 
    if (max(apply(FD, 1, sum)) >= 2*(k+1)){
      break
    }
    
    # --- Represent current decision region and frontier
    if(represent){
      plot(score[y==1], xi[y==1], col = 2, main = iter)
      abline(v = d1)
      abline(h = d2)
      points(data.frame(d1[region[,1]], d2[region[,2]]), pch = 16, col = 5, cex = 1.5)
      points(data.frame(d1[FD_prev[,1]], d2[FD_prev[,2]]), pch = 1, cex = 2.5, col = 4)
      points(data.frame(d1[FD[,1]], d2[FD[,2]]), pch = 16)
    }
    
    ### --- Loss function adding the frontier points --- ###
    lf_new <- c()
    for (i in 1:nrow(FD)){
      region_new <- rbind(region, FD[i,])
      lf_new <- c(lf_new, loss_function_grid(score = score, 
                                             xi = xi, 
                                             cost_matrix = cost_matrix, 
                                             region = region_new)$lf
      )
    }
    
    # --- Add the point leading to the greatest improvement
    # --- If no improvement respect previous step, consider the next surrounding points in the grid
    region_prev <- region # Save current region in case PP threshold is surpassed
    new_point <- which(lf_new == min(lf_new))
    if (lf_new[new_point[1]] < lf){
      region <- rbind(region, FD[new_point,])
      for (i in 1:nrow(region)){
        for (j in 1:nrow(region)){
          if(region[i,1] <= region[j,1] & region[i,2] < region[j,2] | 
             region[i,1] < region[j,1] & region[i,2] <= region[j,2]){
            region[i,] <- c(0,0)
          }
        }
      }

      region <- region[region[,1] > 0,]
      
      lf <- lf_new[new_point[1]]
      # FD_prev <- region
      
      # --- Consider the frontier --- #
      FD_prev <- region
      for (i in 1:nrow(FD_prev)){
        FD_prev <- rbind(
              FD_prev,
              data.frame(d1 = region[i,1], d2 = index[index <= region[i,2]]),
              data.frame(d1 = index[index <= region[i,1]], d2 = region[i,2])
            )
      }
      
      FD_prev <- unique(FD_prev)
      
      # --- Take all the frontier --- #
      # FD_prev <- data.frame()
      # for (i in 1:(dim(region)[1])){
      #   FD_prev <- rbind(
      #     FD_prev,
      #     data.frame(d1 = region[i,1], d2 = index[index <= region[i,2]]),
      #     data.frame(d1 = index[index <= region[i,1]], d2 = region[i,2])
      #   )
      # }
      
    } else{
      FD_prev <- FD
    }
    
    for (i in 1:nrow(FD_prev)){
      for (j in 1:nrow(FD_prev)){
        if(FD_prev[i,1] < FD_prev[j,1] & FD_prev[i,2] < FD_prev[j,2]){
          FD_prev[i,] <- c(0,0)
        }
      }
    }
    FD_prev <- FD_prev[FD_prev[,1] > 0,]
    
    # --- Check if PP threshold is surpassed --- #
    FD_max <- head(FD %>% 
      mutate(d3 = d1 + d2) %>%
      slice_max(d3), 1)[,1:2]
      
    
    pp_act <- loss_function_grid(score = score, 
                                 xi = xi, 
                                 cost_matrix = cost_matrix, 
                                 region = rbind(region, FD_max))$pp 
    
    # print(loss_function_grid(score = score,
    #                          xi = xi,
    #                          cost_matrix = cost_matrix,
    #                          region = region)$pp)

    if(pp_act > pp_max){
      region <- region_prev
      break
    }
    
    iter <- iter + 1
    
  } # End loop
  
  ### --- Estimated region and metrics --- ###
  region <- region[region$d1 > 1 & region$d2 > 1,]
  region <- data.frame(d1 = d1[region[,1]], d2 = d2[region[,2]])
  return(list(region = region, lf = lf, sav = 1 - lf/sum(cost_matrix[,1])))
}


# a_cost <- 0.047
# b_cost <- .01
# c_cost <- .6
# n <- 5000
# 
# score <- c(runif(n), .5)
# y <- c(as.numeric(plogis(-2 + 3*score[1:n] + rnorm(n)) > .8), 1); mean(y)
# xi <- c(rchisq(n, df = 2), 18)
# 
# # --- First/second column: Cost when y_hat = 0 / 1
# cost_matrix <- matrix(nrow = length(score), ncol = 2)
# cost_matrix[, 1] <- ifelse(y == 1, c_cost*xi, 0)
# cost_matrix[, 2] <- ifelse(y == 1, b_cost, b_cost + a_cost*xi)

# algo_2DDR(score = score, xi = xi, y = y, cost_matrix = cost_matrix, k = 10, 
#           pp_max = 1, represent = T, regular = T)
# algo_2DDR(score = score, xi = xi, y = y, cost_matrix = cost_matrix, k = 10, 
#           pp_max = .05, represent = F, regular = F)
  

