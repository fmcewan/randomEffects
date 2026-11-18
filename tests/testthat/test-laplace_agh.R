# Testing both laplace and agh methods simultanously

# Test data
integration_cd <- list(
    y = c(0, 0, 0, 0,
            1, 1, 1, 1,
            0, 1, 0, 1),
    X = cbind(1, rep(c(-1, -0.5, 0.5, 1), 3)),
    group = rep(1:3, each = 4),
    J = 3
)
integration_theta <- c(-0.5, 1, log(0.8))

# TEST 1) Check the returned structure and finite values.
test_that("both methods return one (named) finite value per practice", {
    results <- list(
        laplace = laplace_loglik(integration_theta, integration_cd),
        agh = agh_loglik(integration_theta, integration_cd, K = 7)
    )

    for (method in names(results)) {
        result <- results[[method]]
        # Standard type, length, finiteness check
        expect_type(result, "double")
        expect_length(result, integration_cd$J)
        expect_true(all(is.finite(result)), info = method)
        # Checking named output matching the group numbers
        expect_equal(names(result), as.character(seq_len(integration_cd$J)), info = method)
    }
})

# TEST 2) Adaptive Gauss-Hermite with one node = Laplace
test_that("one-node adaptive quadrature equals Laplace", {
    laplace <- laplace_loglik(integration_theta, integration_cd)
    agh_one <- agh_loglik(integration_theta, integration_cd, K = 1)
    expect_equal(agh_one, laplace, tolerance = 1e-8)
})

# TEST 3) Comparing AGH with independent numerical integration
test_that("adaptive quadrature agrees with integrate for every practice", {
        
    theta <- integration_theta; cd <- integration_cd
    beta <- theta[1:2]; sigma <- exp(theta[3])
    agh <- agh_loglik(theta, cd, K = 40)

    # Independent integrations per cluster to provide an independent reference 
    for (j in seq_len(cd$J)) {  # (again, looping over clusters just for testing)

        # Filtering only patients being currently evaluated and storing fixed predictions and outcomes
        patients <- cd$group == j 
        fixed <- drop(cd$X[patients, , drop = FALSE] %*% beta) 
        outcomes <- cd$y[patients]

        h_j <- function(u) {
        sum(dbinom(outcomes, size = 1, prob = plogis(fixed + u), log = TRUE)) +
        dnorm(u, mean = 0, sd = sigma, log = TRUE)
        }

        # General logic adapted from practical 7 L1
        # To pass to integrate, vapply() evaluates h_j() at each supplied value
        integrand <- function(u) { vapply(u, function(value) exp(h_j(value)), numeric(1)) }
        # Reference integral for comparison 
        reference <- integrate(integrand,lower = -Inf, upper = Inf, rel.tol = 1e-10)$value

        # agh returns log-integrals, so log-ing the reference too
        # Removing the cluster label for this test, just comparing value output
        expect_equal(unname(agh[j]), log(reference), tolerance = 1e-6, info = paste("practice", j))
    }
})

# TEST 4) Checking numerical stability with larger sigma 
test_that("both methods give finite values with larger sigma", {
    theta <- c(-0.5, 1, log(5)) # Theta with larger sigma
    laplace <- laplace_loglik(theta, integration_cd)
    agh <- agh_loglik(theta, integration_cd, K = 7)
    expect_true(all(is.finite(laplace)))
    expect_true(all(is.finite(agh)))
})

# TEST 5) Checking that inputting an invalid numbers of nodes gives an informative error
test_that("adaptive quadrature rejects invalid K", {
    # Listing a range of invalid K and checking each
    invalid_K <- list(0, -1, 1.5, Inf, NA_real_, c(2, 3), "seven")
    for (K in invalid_K) {
        expect_error(
        agh_loglik(integration_theta, integration_cd, K = K),
        regexp = "K must be a positive integer"
        )
    }
})