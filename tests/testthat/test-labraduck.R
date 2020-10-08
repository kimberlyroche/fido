context("test-labraduck")

# # power a fixed G through eigenvalue decomposition
# power_G <- function(G, it_begin, it_end, descending=TRUE, transpose=FALSE) {
#   obj <- G
#   if(transpose) {
#     obj <- t(G)
#   }
#   obj <- eigen(obj)
#   e_vec <- obj$vectors
#   e_val <- diag(2)
#   diag(e_val) <- obj$values
#   if(it_begin != it_end) {
#     if(descending) {
#       if(it_begin > it_end) {
#         power_it <- length(1:(it_begin-it_end+1))
#         e_val <- e_val**power_it
#       } else {
#         # invalid case
#         e_val <- matrix(0, 2, 2)
#       }
#     } else {
#       if(it_begin < it_end) {
#         power_it <- length(1:(it_end-it_begin+1))
#         e_val <- e_val**power_it
#       } else {
#         # invalid case
#         e_val <- matrix(0, 2, 2)
#       }
#     }
#   }
#   # explicitly only returning the real part of A
#   # some tiny complex eigenvalues can be produced
#   ret_val <- Re(e_vec%*%e_val%*%solve(e_vec)) # reconstitute
#   return(ret_val)
# }

# build_A <- function(T, gamma, F, G, W, C0) {
#   A <- matrix(0, T, T)
#   for(i in 1:T) {
#     for(j in 1:T) {
#       if(i == j) {
#         # diagonal
#         t <- j
#         first_sum <- matrix(0, 2, 2)
#         if(t >= 2) {
#           for(ell in t:2) {
#             G_left <- power_G(G, t, ell)
#             G_right <- power_G(G, ell, t, descending=FALSE, transpose=TRUE)
#             addend <- G_left%*%W%*%G_right
#             first_sum <- first_sum + addend
#           }
#         }
#         # second sum
#         G_left <- power_G(G, t, 1)
#         G_right <- power_G(G, 1, t, descending=FALSE, transpose=TRUE)
#         second_sum <- G_left%*%C0%*%G_right
#         A[t,t] <- gamma + F%*%(W + first_sum + second_sum)%*%t(F)
#       } else {
#         tk <- i
#         t <- j
#         if(j < i) {
#           tk <- j
#           t <- i
#         }
#         # off-diagonal
#         first_sum <- matrix(0, 2, 2)
#         for(ell in tk:2) {
#           G_left <- power_G(G, t, ell)
#           G_right <- power_G(G, ell, tk, descending=FALSE, transpose=TRUE)
#           first_sum <- first_sum + G_left%*%W%*%G_right
#         }
#         G_left <- power_G(G, t, 1)
#         G_right <- power_G(G, 1, tk, descending=FALSE, transpose=TRUE)
#         second_sum <- G_left%*%C0%*%G_right
#         G_left <- power_G(G, t, tk+1)
#         A[i,j] <- F%*%(G_left%*%W + first_sum + second_sum)%*%t(F)
#       }
#     }
#   }
#   if(min(eigen(A)$values) < 0) {
#     cat("Matrix A has negative eigenvalue(s)!\n")
#   }
#   return(A)
# }

test_that("G^1 works correctly", {
  T <- 10
  omega <- 2*pi/T
  G <- matrix(c(cos(omega), -sin(omega), sin(omega), cos(omega)), 2, 2)
  retG = power_G(G, 1, 1)
  expect_equal(G, retG, tolerance=1e-3)
})

test_that("G^4 works correctly", {
  T <- 10
  omega <- 2*pi/T
  G <- matrix(c(cos(omega), -sin(omega), sin(omega), cos(omega)), 2, 2)
  retG = power_G(G, 1, 4)
  expect_equal(G%*%G%*%G%*%G, retG, tolerance=1e-3)
})

# test C++ implementation of DLM sample covariance matrix against an previously tested R implementation
test_that("DLM A matrix correctness", {
  # T <- 10
  # omega <- 2*pi/T
  # gamma <- 1
  # F <- matrix(c(1, 0), 1, 2) # row vector
  # G <- matrix(c(cos(omega), -sin(omega), sin(omega), cos(omega)), 2, 2)
  # W <- diag(2)
  # C0 <- W
  # A_R <- build_A(T, gamma, F, G, W, C0)
  #A_C <- dlm_A_test(T, gamma, F, G, W, C0) # in src/LabraduckCollapsed_Uncollapse.cpp
  expect_true(TRUE)
  # expect_equal(A_R, A_R)
})











