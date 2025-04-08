# generate <- function(n, alpha = 0.5, missing = TRUE, A = NULL, V1_1 = NULL, V1_2 = NULL, V2_1 = NULL, V2_2 = NULL) {
#   S <- stats::rbinom(n, 1, alpha)
#   W <- stats::rbinom(n, 1, 0.2 + 0.7*S)
#   
#   if (is.null(V1_1) & is.null(V1_2)) {
#     V1_1 <- stats::rbinom(n, 1, 0.5 - S*0.2 + 0.15*W)
#     V1_2 <- stats::rbinom(n, 1, 0.25 + S*0.2 - 0.15*W)
#   }
#   
#   if (is.null(V2_1) & is.null(V2_2)) {
#     V2_1 <- stats::rbinom(n, 1, 0.3 + 0.2*V1_1 + 0.2*W)
#     V2_2 <- stats::rbinom(n, 1, 0.7 - 0.2*V1_1 - 0.1*V1_2 + 0.05*W)
#   }
#   
#   if (is.null(A)) {
#     A <- stats::rbinom(n, 1, 0.5)
#   }
#   
#   # V1 and V2 modify the effect of A on Y
#   logit_Y <- -1 - 0.8*W + 0.5*V1_1 - 0.02*V1_2 - 0.5*V2_2 + 0.7*V2_1 + (0.3*V1_1 + 0.5*V2_1 - 1.3*V1_2 - 0.75*V2_2) * A
#   Y <- stats::rbinom(n, 1, stats::plogis(logit_Y))
#   
#   if (missing) {
#     V2_1[S == 0] <- NA_real_
#     V2_2[S == 0] <- NA_real_
#   }
#   
#   data.frame(S, W, V1_1, V1_2, V2_1, V2_2, A, Y)
# }

# generate <- function(n, alpha = 0.5, missing = TRUE, A = NULL, V1 = NULL, V2 = NULL) {
#   S <- stats::rbinom(n, 1, alpha)
#   W <- stats::rbinom(n, 1, 0.4 + 0.2*S)
#   
#   if (is.null(V1)) V1 <- stats::rbinom(n, 1, 0.5 - S*0.2 + 0.15*W)
#   if (is.null(V2)) V2 <- stats::rbinom(n, 1, 0.3 + 0.2*V1 + 0.2*W)
#   if (is.null(A)) A <- stats::rbinom(n, 1, 0.5)
#   
#   # V1 and V2 modify the effect of A on Y
#   logit_Y <- -1 - 0.8*W + 0.5*V1 + 0.7*V2 + (0.6*V1 + 0.2*V2) * A
#   Y <- stats::rbinom(n, 1, stats::plogis(logit_Y))
#   
#   if (missing) {
#     V2[S == 0] <- NA_real_
#   }
#   
#   data.frame(S, W, V1, V2, A, Y)
# }

generate <- function(n, alpha, missing_V2 = TRUE, A) {
  # S: Binary covariate with P(S=1) = alpha
  S <- rbinom(n, 1, alpha)
  
  # A: Binary covariate with P(A=1) = 0.5
  if (missing(A)) A <- rbinom(n, 1, 0.5)
  
  # --- 2. Generate W ~ S ---
  # W consists of 2 variables: W1 (numeric), W2 (binary)
  W1 <- rnorm(n, mean = 1 + 0.5 * S, sd = 1.5)
  W2 <- rbinom(n, 1, plogis(-0.5 + 1.5 * S))
  
  # --- 3. Generate V1 ~ W, S ---
  # V1 consists of 3 binary variables: V1_1, V1_2, V1_3
  prob_V1_1 <- plogis(0.2 + 0.3*S - 0.2*W1 + 0.5*W2)
  V1_1 <- rbinom(n, 1, prob_V1_1)
  
  prob_V1_2 <- plogis(-0.1 - 0.4*S + 0.1*W1 - 0.6*W2)
  V1_2 <- rbinom(n, 1, prob_V1_2)
  
  prob_V1_3 <- plogis(0.0 + 0.1*S + 0.3*W1 + 0.2*W2 - 0.1*W1*S) # Example interaction
  V1_3 <- rbinom(n, 1, prob_V1_3)
  
  # --- 4. Generate V2 ~ V1, W (Missing if S=1) ---
  # V2 is a single binary variable
  prob_V2_raw <- plogis(-0.5 + 0.6*V1_1 - 0.4*V1_2 + 0.3*V1_3 + 0.1*W1 - 0.2*W2)
  
  # Generate the potential V2 value
  V2 <- rbinom(n, 1, prob_V2_raw)
  
  # --- 5. Generate Y ~ W, V1, V2, A (with interactions V1:A, V2:A) ---
  eta_Y <- (-1.5                # Intercept
            + 0.3 * W1          # Effect of W1
            - 0.4 * W2          # Effect of W2
            + 0.1 * A           # Effect of A
            + 0.5 * V1_1        # Effect of V1_1
            - 0.8 * V1_2        # Effect of V1_2
            + 0.2 * V1_3        # Effect of V1_3
            + 1.0 * V1_1 * A    # Interaction V1_1:A
            - 1.2 * V1_2 * A    # Interaction V1_2:A
            + 0.5 * V1_3 * A    # Interaction V1_3:A
            + 0.9 * V2          # Effect of V2
            + 1.2 * V2 * A)     # Interaction V2:A (to be added only if S == 0
  
  # Generate Y
  Y <- rbinom(n, 1, plogis(eta_Y))
  # Set V2 to NA if S == 1
  if (missing_V2) V2 <- ifelse(S == 0, NA, V2) 
  
  # --- 6. Combine into a data.frame ---
  data.frame(
    S = S,
    A = A,
    W1 = W1,
    W2 = W2,
    V1_1 = V1_1,
    V1_2 = V1_2,
    V1_3 = V1_3,
    V2 = V2,
    Y = Y
  )
}
