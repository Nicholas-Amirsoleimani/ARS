############ Unit tests of individual functions ############

## Load testthat
library(testthat)
library(pracma)

## grad_construct() tests

# test with normal input
test_that("grad_construct() works with normal input", {
  m <- 1
  h <- function(x) {
    log(dexp(x, m))
  }
  
  expect_type(grad_construct(h), 'closure')
  expect_is(grad_construct(h), 'function')
})


# test with derivative that doesn't exist
test_that("grad_construct() doesn't work with unevaluatable derivative", {
  x <- -1
  h <- function(x) {
    log(x)
  }
  
  expected_msg <- "!is.nan(derivative) is not TRUE"
  expect_error(supresssWarnings(grad_construct(h)(x), regexp = expected_msg, fixed = TRUE))
})

## get_abscissae() tests

# test with normalized normal log density and inf bounds
test_that("get_abscissae() works with normal input inf bounds", {
  m <- 1
  v <- 5
  ld <- function(x) {
    log(dnorm(x, m, v))
  }
  der <- grad_construct(ld)
  bounds <- c(-Inf,Inf)
  
  # we expect one positive and one negative deriv if -inf,inf bounds 
  expect_gt(pracma::grad(ld, get_abcissae(ld, der, bounds)[1]), 0)
  expect_lt(pracma::grad(ld, get_abcissae(ld, der, bounds)[3]), 0)
  expect_length(get_abcissae(ld, der, bounds), 3)
  expect_type(get_abcissae(ld, der, bounds), 'double')
})

# test with normalized normal log density and -inf bound
test_that("get_abscissae() works with normal input and -inf bound", {
  m <- 1
  v <- 5
  library(truncnorm)
  ld <- function(x) {
    return(log(dtruncnorm(x, b=4)))
  }
  der <- grad_construct(ld)
  bounds <- c(-Inf,4)
  
  # if lower bound is negative inf there must be one positive deriv 
  expect_gt(pracma::grad(ld, get_abcissae(ld, der, bounds)[1]), 0)
  expect_length(get_abcissae(ld, der, bounds), 3)
  expect_type(get_abcissae(ld, der, bounds), 'double')
})

# test with normalized normal log density and +inf bound
test_that("get_abscissae() works with normal input and +inf bound", {
  m <- 1
  v <- 5
  library(truncnorm)
  ld <- function(x) {
    return(log(dtruncnorm(x, a=-1)))
  }
  der <- grad_construct(ld)
  bounds <- c(-1, Inf)
  
  # if upper bound is inf there must be one negative deriv
  expect_lt(pracma::grad(ld, get_abcissae(ld, der, bounds)[3]), 0)
  expect_length(get_abcissae(ld, der, bounds), 3)
  expect_type(get_abcissae(ld, der, bounds), 'double')
})

# test with unnormalized log density 
test_that("get_abscissae() works with an unnormalized density", {
  m <- 10
  v <- 4
  n <- 5
  ld <- function(x) {
    d <- log(dnorm(x, m, v))
    density <- d * n
    return(density)
  }
  der <- grad_construct(ld)
  bounds <- c(-Inf,Inf)
  
  expect_gt(pracma::grad(ld, get_abcissae(ld, der, bounds)[1]), 0)
  expect_lt(pracma::grad(ld, get_abcissae(ld, der, bounds)[3]), 0)
  expect_length(get_abcissae(ld, der, bounds), 3)
  expect_type(get_abcissae(ld, der, bounds), 'double')
})

# test with non-infinite bounds
test_that("get_abscissae() works with non-infinite bounds", {
  m <- 1
  v <- 3
  n <- 5
  ld <- function(x) {
    d <- log(dnorm(x, m, v))
    density <- d * n
    return(density)
  }
  bounds <- c(-10,10)
  der <- grad_construct(ld)
  
  expect_length(get_abcissae(ld, der, bounds), 3)
  expect_type(get_abcissae(ld, der, bounds), 'double')
})

## initialize_data() tests

# test with normal inputs 
test_that("initialize_data() works with normal input", {
  m <- 1
  v <- 3
  ld <- function(x) {
    log(dnorm(x, m, v))
  }
  der <- grad_construct(ld)
  bounds <- c(-Inf,Inf)
  ab <- get_abcissae(ld, der, bounds)
  
  expect_length(initialize_data(ab, ld, der, bounds), 9)
  expect_type(initialize_data(ab, ld, der, bounds), 'list')
})


## k_samples() tests

# test with normal inputs 
test_that("k_samples() works with normal input", {
  m <- 1
  v <- 3
  ld <- function(x) {
    log(dnorm(x, m, v))
  }
  der <- grad_construct(ld)
  bounds <- c(-Inf,Inf)
  ab <- get_abcissae(ld, der, bounds)
  num_samples <- 100
  accepted_points <- rep(0, num_samples)
  accepted_count <- 0
  init <- initialize_data(ab, ld, der, bounds)
  upper_x <- init[[1]]
  upper_y <- init[[2]]
  upper_slopes <- init[[3]]
  upper_ints <- init[[4]]
  upper_areas <- init[[5]]
  lower_x <- init[[6]]
  lower_y <- init[[7]]
  lower_slopes <- init[[8]]
  lower_ints <- init[[9]]
  n <- length(ab)
  k <- 1000
  
  expect_length(k_samples(k, accepted_points, accepted_count, upper_areas, upper_x, upper_slopes, upper_ints, n, bounds), k)
  expect_type(k_samples(k, accepted_points, accepted_count, upper_areas, upper_x, upper_slopes, upper_ints, n, bounds), 'double')
})

## find_first_points() tests

# test with normal inputs
test_that("find_first_points() works with normal input", {
  m <- 1
  v <- 3
  dens <- function(x) {
    return(dnorm(x, 5, 3))
  }
  ld <- function(x) {
    log(dnorm(x, m, v))
  }
  der <- grad_construct(ld)
  bounds <- c(-Inf,Inf)
  ab <- get_abcissae(ld, der, bounds)
  num_samples <- 1000
  accepted_points <- rep(0, num_samples)
  accepted_count <- 0
  init <- initialize_data(ab, ld, der, bounds)
  upper_x <- init[[1]]
  upper_y <- init[[2]]
  upper_slopes <- init[[3]]
  upper_ints <- init[[4]]
  upper_areas <- init[[5]]
  lower_x <- init[[6]]
  lower_y <- init[[7]]
  lower_slopes <- init[[8]]
  lower_ints <- init[[9]]
  n <- length(ab)
  k <- 1000
  x_samples <- k_samples(k, accepted_points, accepted_count, upper_areas, upper_x, upper_slopes, upper_ints, n, bounds)
  
  # expect correct length 
  expect_length(find_first_points(x_samples, k, dens, num_samples, accepted_points, accepted_count, lower_x, lower_y, lower_slopes, upper_x, upper_y, upper_slopes)[[1]], k)
  expect_length(find_first_points(x_samples, k, dens, num_samples, accepted_points, accepted_count, lower_x, lower_y, lower_slopes, upper_x, upper_y, upper_slopes)[[2]], 1)
  expect_length(find_first_points(x_samples, k, dens, num_samples, accepted_points, accepted_count, lower_x, lower_y, lower_slopes, upper_x, upper_y, upper_slopes)[[3]], 1)
  expect_length(find_first_points(x_samples, k, dens, num_samples, accepted_points, accepted_count, lower_x, lower_y, lower_slopes, upper_x, upper_y, upper_slopes)[[4]], 1)
  
  # expect list of vectors 
  expect_type(find_first_points(x_samples, k, dens, num_samples, accepted_points, accepted_count, lower_x, lower_y, lower_slopes, upper_x, upper_y, upper_slopes), 'list')
  expect_is(find_first_points(x_samples, k, dens, num_samples, accepted_points, accepted_count, lower_x, lower_y, lower_slopes, upper_x, upper_y, upper_slopes)[[1]], "numeric")
  expect_is(find_first_points(x_samples, k, dens, num_samples, accepted_points, accepted_count, lower_x, lower_y, lower_slopes, upper_x, upper_y, upper_slopes)[[2]], "numeric")
  expect_is(find_first_points(x_samples, k, dens, num_samples, accepted_points, accepted_count, lower_x, lower_y, lower_slopes, upper_x, upper_y, upper_slopes)[[3]], "numeric")
  expect_is(find_first_points(x_samples, k, dens, num_samples, accepted_points, accepted_count, lower_x, lower_y, lower_slopes, upper_x, upper_y, upper_slopes)[[4]], "numeric")
})


############ Tests of primary function (using Kolmogorov-Smirnov) ############

# Test that ARS doesn't take invalid inputs

test_that("Error for non-numeric num_samples", {
          expect_error(ARS("hello", dunif, c(-Inf, Inf)),
                       "`num_samples` must be a double")})

test_that("Error for Inf if abcissae", {
  expect_error(ARS(10, dunif, c(-Inf, Inf), abcissae = c(-Inf, 5)),
               "`abcissae` cannot contain `Inf` of `-Inf`")})

# beta density
set.seed(1)
f <- function(x) {
  return(dbeta(x,2,2))
}
# ARS Function
samples <- ARS(10000, dens = f, bounds = c(0,1))
hist(samples, breaks = 20, freq = FALSE)
# Real Function
d <- dbeta(seq(0,1,length.out=1000),2,2)
lines(seq(0,1,length.out=1000),d)

# Kolmogorov-Smirnov Test
test_that("ARS() works with beta density", {
  set.seed(1)
  f <- function(x) {
    return(dbeta(x,2,2))
  }
  samples <- ARS(10000, dens = f, bounds = c(0,1))
  
  expect_gt(ks.test(samples,rbeta(10000,2,2))$p.value, 0.05)
})

# normal density
set.seed(1)
f <- function(x) {
  return(dnorm(x, 5, 3))
}
# ARS Function
samples <- ARS(10000, dens = f, bounds = c(-Inf,Inf))
hist(samples, breaks = 20,freq = FALSE)
# Real Function
d <- dnorm(seq(-10,20,length.out=1000),5,3)
lines(seq(-10,20,length.out=1000),d)

# Kolmogorov-Smirnov Test
test_that("ARS() works with normal density", {
  set.seed(1)
  f <- function(x) {
    return(dnorm(x, 5, 3))
  }
  samples <- ARS(10000, dens = f, bounds = c(-Inf,Inf))
  
  expect_gt(ks.test(samples,rnorm(10000,5,3))$p.value, 0.05)
})

# uniform density
set.seed(1)
f <- function(x) {
  return(dunif(x, 5, 14))
}
# ARS Function
samples <- ARS(10000, dens = f, bounds = c(5,14))
hist(samples, breaks = 20,freq = FALSE)
# Real Function
d <- dunif(seq(5,14,length.out=1000),5,14)
lines(seq(5,14,length.out=1000),d)

# Kolmogorov-Smirnov Test
test_that("ARS() works with uniform density", {
  set.seed(1)
  f <- function(x) {
    return(dunif(x, 5, 14))
  }
  samples <- ARS(10000, dens = f, bounds = c(5,14))
  
  expect_gt(ks.test(samples,runif(10000,5,14))$p.value, 0.05)
})


# exponential density
set.seed(1)
f <- function(x) {
  return(dexp(x, 2))
}
# ARS Function
samples <- ARS(10000, dens = f, bounds = c(0,Inf))
hist(samples, breaks = 20,freq = FALSE)
# Real Function
d <- dexp(seq(0,10,length.out=1000),2)
lines(seq(0,10,length.out=1000),d)

# Kolmogorov-Smirnov Test
test_that("ARS() works with exponential density", {
  set.seed(1)
  f <- function(x) {
    return(dexp(x, 2))
  }
  samples <- ARS(10000, dens = f, bounds = c(0,Inf))
  
  expect_gt(ks.test(samples,rexp(10000,2))$p.value, 0.05)
})


# gamma density
set.seed(1)
f <- function(x) {
  return(dgamma(x, 10))
}
# ARS Function
samples <- ARS(10000, dens = f, bounds = c(0,Inf))
hist(samples, breaks = 20,freq = FALSE)
# Real Function
d <- dgamma(seq(0,30,length.out=1000),10)
lines(seq(0,30,length.out=1000),d)

# Kolmogorov-Smirnov Test
test_that("ARS() works with gamma density", {
  set.seed(1)
  f <- function(x) {
    return(dgamma(x, 10))
  }
  samples <- ARS(10000, dens = f, bounds = c(0,Inf))
  
  expect_gt(ks.test(samples,rgamma(10000,10))$p.value, 0.05)
})


# truncated normal density
library(truncnorm)
set.seed(1)
f <- function(x) {
  return(truncnorm::dtruncnorm(x, a=-1, b=4))
}
# ARS Function
samples <- ARS(10000, dens = f, bounds = c(-1,4))
hist(samples, breaks = 20,freq = FALSE)
# Real Function
d <- dtruncnorm(seq(-1,4,length.out=1000), a=-1, b=4)
lines(seq(-1,4,length.out=1000),d)

# Kolmogorov-Smirnov Test
test_that("ARS() works with truncated normal density", {
  set.seed(1)
  f <- function(x) {
    return(truncnorm::dtruncnorm(x, a=-1, b=4))
  }
  samples <- ARS(10000, dens = f, bounds = c(-1,4))
  
  expect_gt(ks.test(samples,truncnorm::rtruncnorm(10000, a=-1, b=4))$p.value, 0.05)
})



## Adding custom derivative function

set.seed(1)
f <- function(x) {
  return(dnorm(x, 5, 3))
}
# derivative of log(dnorm(x,5,3))
d_func <- function(x){
  return(-(x-5)/3^2)
}
# ARS Function
samples <- ARS(10000, dens = f, bounds = c(-Inf,Inf), deriv_func = d_func)
hist(samples, breaks = 20,freq = FALSE)
# Real Function
d <- dnorm(seq(-10,20,length.out=1000),5,3)
lines(seq(-10,20,length.out=1000),d)

# Kolmogorov-Smirnov Test
test_that("ARS() works with custom derivative function", {
  set.seed(1)
  f <- function(x) {
    return(dnorm(x, 5, 3))
  }
  # derivative of log(dnorm(x,5,3))
  d_func <- function(x){
    return(-(x-5)/3^2)
  }
  samples <- ARS(10000, dens = f, bounds = c(-Inf,Inf), deriv_func = d_func)
  
  expect_gt(ks.test(samples,rnorm(10000, 5, 3))$p.value, 0.05)
})


## Adding custom starting abcissae (requires at least 3 points)

set.seed(1)
f <- function(x) {
  return(dexp(x, 2))
}
# ARS Function
samples <- ARS(10000, dens = f, bounds = c(0,Inf), abcissae = c(1,2,5))
hist(samples, breaks = 20,freq = FALSE)
# Real Function
d <- dexp(seq(0,10,length.out=1000),2)
lines(seq(0,10,length.out=1000),d)

# Kolmogorov-Smirnov Test
test_that("ARS() works with custom abscissae", {
  set.seed(1)
  f <- function(x) {
    return(dexp(x, 2))
  }
  samples <- ARS(10000, dens = f, bounds = c(0,Inf), abcissae = c(1,2,5))
  
  expect_gt(ks.test(samples,rexp(10000,2))$p.value, 0.05)
})



# Test that main function throws correct error when bounds are -inf,inf 
# and abscissae are all positive 
test_that("ARS() throws expected error when abscissae all positive and bounds -inf to inf", {
  set.seed(1)
  f <- function(x) {
    return(dnorm(x))
  }
  
  expected_msg <- "No elements of upper_slopes > 0 are true"
  expect_error(ARS(10000, dens = f, bounds = c(-Inf,Inf), abcissae = c(1,2,3)), regexp = expected_msg, fixed = TRUE)
})

# show that this works when abscissae set to {-1, 2, 3}
set.seed(1)
f <- function(x) {
  return(dnorm(x))
}
# ARS Function
samples <- ARS(10000, dens = f, bounds = c(-Inf,Inf), abcissae = c(-1,2,3))
hist(samples, breaks = 20,freq = FALSE)
# Real Function
d <- dnorm(seq(-10,20,length.out=1000))
lines(seq(-10,20,length.out=1000),d)

# Kolmogorov-Smirnov Test
test_that("ARS() works with adjusted abscissae", {
  set.seed(1)
  f <- function(x) {
    return(dnorm(x))
  }
  samples <- ARS(10000, dens = f, bounds = c(-Inf,Inf), abcissae = c(-1,2,3))
  
  expect_gt(ks.test(samples,rnorm(10000))$p.value, 0.05)
})
          


