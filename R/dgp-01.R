generate <- function(n, missing_V2 = TRUE, A = NULL) {
  W1 <- rbinom(n, 1, 0.33)
  W2 <- rbeta(n, 2, 2, ncp = 0)
  
  # V1 consists of 3 binary variables: V1_1, V1_2, V1_3
  V1_1 <- rbinom(n, 1, plogis(0.5 - 0.2*W1 + 0.15*W2))
  V1_2 <- rbinom(n, 1, plogis(-0.3 + 0.1*W1 - 0.6*W2))
  V1_3 <- rbinom(n, 1, plogis(0.1 + 0.3*W1 + 0.2*W2))
  
  # V2 is a single binary variable
  V2 <- rbinom(n, 1, plogis(-0.5 + 0.6*V1_1 - 0.4*V1_2 + 0.3*V1_3 + 0.1*W1 - 0.2*W2))
  
  S <- rbinom(n, 1, plogis(0 + 0.5*W1 - 0.3*W2 + 0.2*V1_1 - 0.4*V1_2 + 0.3*V1_3))
  
  if (is.null(A)) {
    A <- stats::rbinom(n, 1, 0.5)
  }
  
  # Generate Y ~ W, V1, V2, A (with interactions V1:A, V2:A) ---
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
            + 1.2 * V2 * A)     # Interaction V2:A
  
  # Generate Y
  Y <- rbinom(n, 1, plogis(eta_Y))
  # Set V2 to NA if S == 1
  if (missing_V2) V2 <- ifelse(S == 0, NA, V2) 
  
  # Combine into a data.frame ---
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
