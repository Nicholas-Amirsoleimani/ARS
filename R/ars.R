library(assertthat)
library(pracma)

####################################################################################

get_z_yint_left <- function(x, h_x, dh_x, left_x, left_y, left_slope) {
  # gets the coordinates of the upper polygon given x and the point left of x
  
  if (abs(dh_x - left_slope) < 1e-8) {
    zx <- x
  } else {
    zx <- (h_x - left_y - x * dh_x + left_x * left_slope) / (left_slope - dh_x)
  }
  zy <- left_y + (zx - left_x) * left_slope
  y_int <- h_x - dh_x * x

  # assert that it exists
  assertthat::assert_that(!is.nan(zx))

  # returns two new z and y-int
  return(list(zx, zy, y_int))
}

####################################################################################

get_z_yint_right <- function(x, h_x, dh_x, right_x, right_y, right_slope) {
  # gets the coordinates of the upper polygon given x and the point right of x
  
  if (abs(dh_x - right_slope) < 1e-8) {
    zx <- x
  } else {
    zx <- (right_y - h_x - right_x * right_slope + x * dh_x) / (dh_x - right_slope)
  }

  zy <- h_x + (zx - x) * dh_x
  y_int <- h_x - dh_x * x

  assertthat::assert_that(!is.nan(zx))

  # returns two new z and y-int
  return(list(zx, zy, y_int))
}

####################################################################################

get_slope_yint <- function(x, h_x, point2_x, point2_y) {
  # gets the slope of the lower polygon given two points
  
  slope <- (h_x - point2_y) / (x - point2_x)
  y_int <- h_x - slope * x

  # ensures slope is not 0
  if (abs(slope) <= 1e-8) {
    if (sign(slope) == -1){
      slope <- -1e-8
    } else {
      slope <- 1e-8
    }
    
  }

  # return slope and y-int
  return(list(slope, y_int))
}

####################################################################################

update_upper <- function(x, position, h_x, dh_x, lower_x, lower_y, upper_x, upper_y, upper_slopes, upper_ints, n, deriv_func) {
  # updates upper polygon
  
  # cycle through different possibilities of location of x
  # it could be the leftmost point, rightmost point, or somewhere in the middle
  if (position == -Inf) {
    z_yint <- get_z_yint_right(x, h_x, dh_x, lower_x[1], lower_y[1], upper_slopes[1])
    
    # update
    upper_x <- append(upper_x, z_yint[[1]], 0)
    upper_y <- append(upper_y, z_yint[[2]], 0)
    upper_slopes <- append(upper_slopes, dh_x, 0)
    upper_ints <- append(upper_ints, z_yint[[3]], 0)
  } else if (position == n) {
    z_yint <- get_z_yint_left(x, h_x, dh_x, lower_x[n], lower_y[n], upper_slopes[n])

    # update
    upper_x <- append(upper_x, z_yint[[1]], n)
    upper_y <- append(upper_y, z_yint[[2]], n)
    upper_slopes <- append(upper_slopes, dh_x, n)
    upper_ints <- append(upper_ints, z_yint[[3]], n)
  } else {
    # left point
    z_yint_left <- get_z_yint_left(x, h_x, dh_x, lower_x[position], lower_y[position], upper_slopes[position])

    # right point
    z_yint_right <- get_z_yint_right(x, h_x, dh_x, lower_x[position + 1], lower_y[position + 1], upper_slopes[position + 1])

    # update
    upper_x[position] <- z_yint_left[[1]]
    upper_x <- append(upper_x, z_yint_right[[1]], position)
    upper_y[position] <- z_yint_left[[2]]
    upper_y <- append(upper_y, z_yint_right[[2]], position)
    upper_slopes <- append(upper_slopes, dh_x, position)
    upper_ints <- append(upper_ints, z_yint_left[[3]], position)
  }

  return(list(upper_x, upper_y, upper_slopes, upper_ints))
}

####################################################################################

update_lower <- function(x, h_x, position, lower_x, lower_y, lower_slopes, lower_ints, n) {
  # updates lower polygon
  
  # cycle through different possibilities
  # it could be the leftmost point, rightmost point, or somewhere in the middle
  if (position == -Inf) {
    slopes_int <- get_slope_yint(x, h_x, lower_x[1], lower_y[1])

    # update
    lower_x <- append(lower_x, x, 0)
    lower_y <- append(lower_y, h_x, 0)
    lower_slopes <- append(lower_slopes, slopes_int[[1]], 0)
    lower_ints <- append(lower_ints, slopes_int[[2]], 0)
  } else if (position == n) {
    slopes_int <- get_slope_yint(x, h_x, lower_x[n], lower_y[n])

    # update
    lower_x <- append(lower_x, x, n)
    lower_y <- append(lower_y, h_x, n)
    lower_slopes <- append(lower_slopes, slopes_int[[1]], n)
    lower_ints <- append(lower_ints, slopes_int[[2]], n)
  } else {
    ## left point
    slopes_int_left <- get_slope_yint(x, h_x, lower_x[position], lower_y[position])

    ## right point
    slopes_int_right <- get_slope_yint(x, h_x, lower_x[position + 1], lower_y[position + 1])

    # update
    lower_x <- append(lower_x, x, position)
    lower_y <- append(lower_y, h_x, position)
    lower_slopes[position] <- slopes_int_left[[1]]
    lower_slopes <- append(lower_slopes, slopes_int_right[[1]], position)
    lower_ints[position] <- slopes_int_left[[2]]
    lower_ints <- append(lower_ints, slopes_int_right[[2]], position)
  }

  return(list(lower_x, lower_y, lower_slopes, lower_ints))
}

####################################################################################

update_areas <- function(n, position, upper_x, upper_slopes, upper_ints, upper_areas, bounds) {
  # updates areas of upper polygon
  
  # cycle through different possibilities
  # it could be the leftmost point, rightmost point, or somewhere in the middle
  if (position == -Inf) {
    area_left <- get_area(bounds[1], upper_x[1], upper_slopes[1], upper_ints[1])
    area_right <- get_area(upper_x[1], upper_x[2], upper_slopes[2], upper_ints[2])

    # update
    upper_areas[1] <- area_right
    upper_areas <- append(upper_areas, area_left, 0)
  } else if (position == n) {
    area_left <- get_area(upper_x[position - 1], upper_x[position], upper_slopes[position], upper_ints[position])
    area_right <- get_area(upper_x[position], bounds[2], upper_slopes[position + 1], upper_ints[position + 1])

    # update
    upper_areas[position] <- area_left
    upper_areas <- append(upper_areas, area_right, position)
  } else {
    if (position == 1) {
      area_left <- get_area(bounds[1], upper_x[1], upper_slopes[1], upper_ints[1])
      area_middle <- get_area(upper_x[1], upper_x[2], upper_slopes[2], upper_ints[2])
      area_right <- get_area(upper_x[2], upper_x[3], upper_slopes[3], upper_ints[3])
    } else if (position == (n - 1)) {
      area_left <- get_area(upper_x[position - 1], upper_x[position], upper_slopes[position], upper_ints[position])
      area_middle <- get_area(upper_x[position], upper_x[position + 1], upper_slopes[position + 1], upper_ints[position + 1])
      area_right <- get_area(upper_x[position + 1], bounds[2], upper_slopes[position + 2], upper_ints[position + 2])
    } else {
      area_left <- get_area(upper_x[position - 1], upper_x[position], upper_slopes[position], upper_ints[position])
      area_middle <- get_area(upper_x[position], upper_x[position + 1], upper_slopes[position + 1], upper_ints[position + 1])
      area_right <- get_area(upper_x[position + 1], upper_x[position + 2], upper_slopes[position + 2], upper_ints[position + 2])
    }

    # update
    upper_areas[position] <- area_left
    upper_areas[position + 1] <- area_right
    upper_areas <- append(upper_areas, area_middle, position)
  }

  return(upper_areas)
}

####################################################################################

get_area <- function(bound1, bound2, slope, y_int) {
  # gets the area of the region under the exponentiated polygon given two bounds
  
  if (bound1 == -Inf) {
    area <- exp(slope * bound2 + y_int) / slope
  } else if (bound2 == Inf) {
    area <- -exp(slope * bound1 + y_int) / slope
  } else {
    area <- (exp(slope * bound2 + y_int) - exp(slope * bound1 + y_int)) / slope
  }

  # assert that area is positive
  assertthat::assert_that(area >= 0)

  return(area)
}

####################################################################################

k_samples <- function(k, accepted_points, accepted_count, upper_areas, upper_x, upper_slopes, upper_ints, n, bounds) {
  # sample k points from upper polygon

  # first, sample the pieces in the piecewise function
  pieces <- sample(n, k, prob = upper_areas, replace = TRUE)
  quantiles <- runif(k) * upper_areas[pieces]
  xvals <- rep(0, k)

  # then, sample from the pieces
  first_inds <- (pieces == 1)
  other_inds <- (!first_inds)

  # find what value of x is required to get the area obtained by the quantiles variable
  if (sum(first_inds) > 0) {
    if (bounds[1] == -Inf) {
      xvals[first_inds] <- (log(upper_slopes[1] * quantiles[first_inds]) - upper_ints[1]) / upper_slopes[1]
    } else {
      xvals[first_inds] <- log(upper_slopes[pieces[first_inds]] * exp(-upper_ints[pieces[first_inds]]) * quantiles[first_inds] + exp(upper_slopes[pieces[first_inds]] * bounds[1])) / upper_slopes[pieces[first_inds]]
    }
  }

  # find what value of x is required to get the area obtained by the quantiles variable
  if (sum(other_inds) > 0) {
    xvals[other_inds] <- log(upper_slopes[pieces[other_inds]] * exp(-upper_ints[pieces[other_inds]]) * quantiles[other_inds] + exp(upper_slopes[pieces[other_inds]] * upper_x[pieces[other_inds] - 1])) / upper_slopes[pieces[other_inds]]
  }
  
  xvals = xvals[is.finite(xvals)]
  
  return(xvals)
}

####################################################################################

get_uk <- function(x, upper_x, upper_y, upper_slopes) {
  # evaluates x at the upper polygon
  
  left_upper <- suppressWarnings(max(which(x > upper_x)))

  if (left_upper == -Inf) {
    uk <- upper_y[1] - (upper_x[1] - x) * upper_slopes[1]
  } else {
    uk <- upper_y[left_upper] + (x - upper_x[left_upper]) * upper_slopes[left_upper + 1]
  }
  return(uk)
}

####################################################################################

get_lk <- function(x, lower_x, lower_y, lower_slopes) {
  # evaluates x at the lower polygon
  
  left_lower <- suppressWarnings(max(which(x > lower_x)))

  if ((left_lower == -Inf) | (left_lower == length(lower_x))) {
    lk <- -Inf
  } else {
    lk <- lower_y[left_lower] + (x - lower_x[left_lower]) * lower_slopes[left_lower]
  }
  return(lk)
}

####################################################################################

find_first_points <- function(x_samples, k, dens, num_samples, accepted_points, accepted_count, lower_x, lower_y, lower_slopes, upper_x, upper_y, upper_slopes) {
  # finds first points that don't require evaluation of f(x) (points under lower polygon)
    
  # set default x as nan in case the lower polygon covers all cases
  x_to_update <- NaN
  h_x <- NaN

  # vectorize random unifs
  ws <- runif(k)

  # cycle through each point
  for (i in 1:length(x_samples)) {
    x <- x_samples[i]
    lk <- get_lk(x, lower_x, lower_y, lower_slopes)
    uk <- get_uk(x, upper_x, upper_y, upper_slopes)
    lower_upper_diff <- lk - uk

    # assert lower_upper_diff is negative
    assertthat::assert_that(lower_upper_diff <= 1e-9)

    if (ws[i] <= exp(lower_upper_diff)) {
      accepted_count <- accepted_count + 1
      accepted_points[accepted_count] <- x
    } else {
      h_x <- h(dens, x)

      # assert it exists
      assertthat::assert_that(!is.na(h_x))

      if (ws[i] <= exp(h_x - uk)) {
        accepted_count <- accepted_count + 1
        accepted_points[accepted_count] <- x
        x_to_update <- x
        break
      } else {
        x_to_update <- x
        break
      }
    }

    # if finished, return and exit
    if (accepted_count == num_samples) {
      return(list(accepted_points, accepted_count, x_to_update, h_x))
    }
  }

  return(list(accepted_points, accepted_count, x_to_update, h_x))
}

####################################################################################

g <- function(dens) {
  # function that finds log of the density
  
  function(x) {
    val <- log(dens(x))
    if (is.infinite(val)) {
      return(1e100)
    } else if (is.infinite(-val)) {
      return(-1e100)
    } else {
      return(val)
    }
  }
}

h <- function(dens, x) {
  # function that returns log of density
  return(log(dens(x)))
}

####################################################################################

grad_construct <- function(h) {
  # finds the derivative of h
  
  function(x) {
    derivative <- pracma::grad(h, x)

    # assert derivative must exist
    assertthat::assert_that(!is.nan(derivative))
    assertthat::assert_that(!is.infinite(derivative))

    # make sure not 0
    if (abs(derivative) <= 1e-8) {
      if (sign(derivative) == -1) {
        derivative <- -1e-8
      } else {
        derivative <- 1e-8
      }
    }

    return(derivative)
  }
}

####################################################################################

get_abcissae <- function(h, deriv_func, bounds) {
  # gets initial abcissae if not given by user
  
  # cycles through different options of bounds since they require different starting requirements
  
  # if bounds are -Inf, Inf
  if ((bounds[1] == -Inf) & (bounds[2] == Inf)) {
    opt <- suppressWarnings(optim(1, h, control = list(fnscale = -1)))
    x_opt <- opt$par
    lower_abs <- suppressWarnings(pracma::newtonRaphson((function(x) deriv_func(x) - 0.8), x_opt))
    upper_abs <- suppressWarnings(pracma::newtonRaphson((function(x) deriv_func(x) + 0.8), x_opt))
    
    # only lower bound is -Inf
  } else if (bounds[1] == -Inf) {
    opt <- suppressWarnings(optim(1, h, control = list(fnscale = -1), upper = bounds[2], method = "L-BFGS-B"))
    x_opt <- opt$par

    if (abs(x_opt - bounds[1]) < 1e-5) {
      return(c(bounds[2] - 3, bounds[2] - 2, bounds[2] - 1))
    } else {
      lower_abs <- suppressWarnings(pracma::newtonRaphson((function(x) deriv_func(x) - 0.8), x_opt))
      upper_abs <- suppressWarnings(pracma::newtonRaphson((function(x) deriv_func(x) + 0.8), x_opt))
    }

    lower_abs <- suppressWarnings(pracma::newtonRaphson((function(x) deriv_func(x) - 0.8), x_opt))
    upper_abs <- suppressWarnings(pracma::newtonRaphson((function(x) deriv_func(x) + 0.8), x_opt))
    
    # only upper bound is Inf
  } else if (bounds[2] == Inf) {
    opt <- suppressWarnings(optim(1, h, control = list(fnscale = -1), lower = bounds[1], method = "L-BFGS-B"))
    x_opt <- opt$par

    if (abs(x_opt - bounds[1]) < 1e-5) {
      return(c(bounds[1] + 1, bounds[1] + 2, bounds[1] + 3))
    } else {
      lower_abs <- suppressWarnings(pracma::newtonRaphson((function(x) deriv_func(x) - 0.8), x_opt))
      upper_abs <- suppressWarnings(pracma::newtonRaphson((function(x) deriv_func(x) + 0.8), x_opt))
    }
    
    # both bounds are not infinite
  } else {
    opt <- suppressWarnings(optim((bounds[1] + bounds[2]) / 2, h, control = list(fnscale = -1), lower = bounds[1], upper = bounds[2], method = "L-BFGS-B"))
    x_opt <- opt$par
    delta <- (bounds[2] - bounds[1]) / 10

    if (abs(x_opt - bounds[1]) < 1e-5) {
      return(c(bounds[1] + delta, bounds[1] + 2 * delta, bounds[1] + 3 * delta))
    } else if (abs(x_opt - bounds[2]) < 1e-5) {
      return(c(bounds[2] - 3 * delta, bounds[2] - 2 * delta, bounds[2] - delta))
    } else {
      return(c(x_opt - delta, x_opt, x_opt + delta))
    }
  }

  return(c(lower_abs$root, x_opt, upper_abs$root))
}

####################################################################################

initialize_data <- function(abcissae, h, deriv_func, bounds) {
  # initializes all data structures that will be used (polygons & areas)
  
  l <- length(abcissae)

  upper_x <- rep(0, l - 1)
  upper_y <- rep(0, l - 1)
  upper_slopes <- rep(0, l)
  upper_ints <- rep(0, l)
  upper_areas <- rep(0, l)
  lower_x <- rep(0, l)
  lower_y <- rep(0, l)
  lower_slopes <- rep(0, l - 1)
  lower_ints <- rep(0, l - 1)

  for (i in 1:l) {
    x <- abcissae[i]
    h_x <- h(x)
    dh_x <- deriv_func(x)

    # assert that they exist
    assertthat::assert_that(!is.na(h_x))

    upper_slopes[i] <- dh_x
    upper_ints[i] <- h_x - dh_x * x
    lower_x[i] <- x
    lower_y[i] <- h_x
  }

  # assert slopes must follow bounds
  if (bounds[1] == -Inf) {
    assertthat::assert_that(any(upper_slopes > 0))
  }

  if (bounds[2] == Inf) {
    assertthat::assert_that(any(upper_slopes < 0))
  }

  # cycle through abcissae
  for (i in 1:l) {
    x <- abcissae[i]

    # initialize points in data structures
    if (i == 1) {
      z_yint <- get_z_yint_right(lower_x[1], lower_y[1], upper_slopes[1], lower_x[2], lower_y[2], upper_slopes[2])
      slope_yint <- get_slope_yint(lower_x[1], lower_y[1], lower_x[2], lower_y[2])

      upper_x[i] <- z_yint[[1]]
      upper_y[i] <- z_yint[[2]]
      upper_areas[i] <- get_area(bounds[1], upper_x[1], upper_slopes[1], upper_ints[1])
      lower_slopes[i] <- slope_yint[[1]]
      lower_ints[i] <- slope_yint[[2]]
    } else if (i == l) {
      upper_areas[i] <- get_area(upper_x[i - 1], bounds[2], upper_slopes[i], upper_ints[i])
    } else {
      z_yint <- get_z_yint_right(lower_x[i], lower_y[i], upper_slopes[i], lower_x[i + 1], lower_y[i + 1], upper_slopes[i + 1])
      slope_yint <- get_slope_yint(lower_x[i], lower_y[i], lower_x[i + 1], lower_y[i + 1])

      upper_x[i] <- z_yint[[1]]
      upper_y[i] <- z_yint[[2]]
      upper_areas[i] <- get_area(upper_x[i - 1], upper_x[i], upper_slopes[i], upper_ints[i])
      lower_slopes[i] <- slope_yint[[1]]
      lower_ints[i] <- slope_yint[[2]]
    }
  }
  return(list(upper_x, upper_y, upper_slopes, upper_ints, upper_areas, lower_x, lower_y, lower_slopes, lower_ints))
}

####################################################################################

ARS <- function(num_samples, dens, bounds, deriv_func = NULL, abcissae = NULL) {
  #' ARS
  #'
  #' @description
  #' This function generates samples from a custom distribution.
  #'
  #' @param num_samples integer -> The sample size.
  #' @param dens function -> The target density. Must be a continuous log-concave function
  #' that takes one input as a float and returns the density as a float. The density does
  #' not need to be normalized.
  #' @param bounds vector -> Bounds of target density. Must be a two-element vector
  #' of the lower and upper bounds of the target density describing where the density
  #' is non-zero.
  #' @param deriv_func (Optional) function -> Custom derivative function. Must be a
  #' function that takes one input as a float and returns the derivative of the log
  #' of the target density.
  #' @param abcissae (Optional) vector -> Starting abcissae locations. Can give starting
  #' abcissae locations. This can be helpful if you want the help the function find where
  #' to start searching. If the lower bound is -Inf, you need at least one abcissae location
  #' where the derivative of the logarithm of the density is positive. If the upper bound is
  #' Inf, then you need at least one location where the derivative of the logarithm of the
  #' density is negative. Requires at least 3 points.
  #'
  #' @return A vector of floats
  #' 
  #' @references
  #' Journal of the Royal Statistical Society. Series C (Applied Statistics), 1992, Vol.
  #' 41, No. 2 (1992), pp. 337-348
  #'
  #' @examples
  #'
  #' norm_samp <- ARS(10000, dens = dnorm, bounds = c(-Inf,Inf))
  #' 
  #' 
  #'
  #' 
  #' set.seed(1)
  #' f <- function(x) {
  #' return(dnorm(x, 5, 3))
  #' }
  #' samples <- ARS(10000, dens = f, bounds = c(-Inf,Inf))
  #' hist(samples, breaks = 20,freq = FALSE)
  #' d <- dnorm(seq(-10,20,length.out=1000),5,3)
  #' lines(seq(-10,20,length.out=1000),d)
  #' ks.test(samples,rnorm(10000,5,3))
  #' 
  #' set.seed(1)
  #' f <- function(x) {
  #'     return(1000*dbeta(x,2,2))
  #' }
  #' samples <- ARS(10000, dens = f, bounds = c(0,1))
  #'
  #' set.seed(1)
  #' f <- function(x) {
  #'     return(dunif(x, 5, 14))
  #' }
  #' samples <- ARS(10000, dens = f, bounds = c(5,14))
  #'
  #' set.seed(1)
  #' f <- function(x) {
  #'     return(dexp(x, 2))
  #' }
  #' samples <- ARS(10000, dens = f, bounds = c(0,Inf))
  #'
  #' set.seed(1)
  #' f <- function(x) {
  #'     return(dgamma(x, 10))
  #' }
  #' samples <- ARS(10000, dens = f, bounds = c(0,Inf))
  #' 
  #' 
  #'
  #' set.seed(1)
  #' f <- function(x) {
  #'     return(dnorm(x, 5, 3))
  #' }
  #' d_func <- function(x){
  #'     return(-(x-5)/3^2)
  #' }
  #' samples <- ARS(10000, dens = f, bounds = c(-Inf,Inf), deriv_func = d_func)
  #'
  #'
  #' 
  #' 
  #' set.seed(1)
  #' f <- function(x) {
  #' return(dexp(x, 2))
  #' }
  #' samples <- ARS(10000, dens = f, bounds = c(0,Inf), abcissae = c(1,2,5))
  #'
  #' @export
  
  # Checks that all arguments are correct type
  
  assertthat::assert_that(length(num_samples) == 1, msg = "`num_samples` must be a single value")
  
  assertthat::assert_that(num_samples > 0, msg = "`num_samples` must be postive")
  
  assertthat::assert_that(is.double(num_samples), msg = "`num_samples` must be a double")
  
  assertthat::assert_that(class(dens) == "function", msg = "`dens` must be of class function")
  
  assertthat::assert_that(is.numeric(bounds), msg = "`bounds` must be a vector of numerics")
  
  if(!is.null(deriv_func)) {
  
    assertthat::assert_that(class(deriv_func) == "function", msg = "`deriv_func` must be of class function")  
  }
  
  
  if(!is.null(abcissae)) {

    assertthat::assert_that(is.numeric(abcissae), msg = "`abcissae` must be a vector of numerics")
    
    assertthat::assert_that(sum(is.infinite(abcissae)) == 0, msg = "`abcissae` cannot contain `Inf` of `-Inf`")
  }
  

  # Rounds num_samples to integer
  num_samples <- round(num_samples)
  
  
  # defining h
  h <- g(dens)

  # defining derivative function if not given by user
  if (is.null(deriv_func)) {
    deriv_func <- grad_construct(h)
  }

  # initialize vector and count of accepted samples
  accepted_points <- rep(0, num_samples)
  accepted_count <- 0

  # intializing abcissae if not given by user
  if (is.null(abcissae)) {
    abcissae <- get_abcissae(h, deriv_func, bounds)
  }
  abcissae <- unique(sort(abcissae))

  # initialize storage vectors
  init <- initialize_data(abcissae, h, deriv_func, bounds)
  upper_x <- init[[1]]
  upper_y <- init[[2]]
  upper_slopes <- init[[3]]
  upper_ints <- init[[4]]
  upper_areas <- init[[5]]
  lower_x <- init[[6]]
  lower_y <- init[[7]]
  lower_slopes <- init[[8]]
  lower_ints <- init[[9]]

  # assert that slopes are decreasing
  assertthat::assert_that(is.unsorted(rev(upper_slopes), FALSE) == FALSE)
  assertthat::assert_that(is.unsorted(rev(lower_slopes), FALSE) == FALSE)

  # initialize number of polygon points
  n <- length(abcissae)

  while (accepted_count < num_samples) {
    ## sample xs from upper polygon k times
    k <- 1000
    x_samples <- k_samples(k, accepted_points, accepted_count, upper_areas, upper_x, upper_slopes, upper_ints, n, bounds)

    if (length(x_samples) > 0) {
      # perform squeezing test
      squeeze_test <- find_first_points(x_samples, k, dens, num_samples, accepted_points, accepted_count, lower_x, lower_y, lower_slopes, upper_x, upper_y, upper_slopes)
      accepted_points <- squeeze_test[[1]]
      accepted_count <- squeeze_test[[2]]
      x_to_update <- squeeze_test[[3]]
      h_x <- squeeze_test[[4]]

      ## update polygons
      if (!is.na(x_to_update)) {

        # find position
        position <- suppressWarnings(max(which(x_to_update > lower_x)))

        # find derivative
        dh_x <- deriv_func(x_to_update)

        # update upper polygon:
        upper_update <- update_upper(x_to_update, position, h_x, dh_x, lower_x, lower_y, upper_x, upper_y, upper_slopes, upper_ints, n, deriv_func)
        upper_x <- upper_update[[1]]
        upper_y <- upper_update[[2]]
        upper_slopes <- upper_update[[3]]
        upper_ints <- upper_update[[4]]

        # update lower polygon:
        lower_update <- update_lower(x_to_update, h_x, position, lower_x, lower_y, lower_slopes, lower_ints, n)
        lower_x <- lower_update[[1]]
        lower_y <- lower_update[[2]]
        lower_slopes <- lower_update[[3]]
        lower_ints <- lower_update[[4]]

        # assert that z values are increasing
        assertthat::assert_that(is.unsorted(upper_x, FALSE) == FALSE)

        # assert that slopes are not increasing
        assertthat::assert_that(is.unsorted(rev(round(upper_slopes, 10)), FALSE) == FALSE)
        assertthat::assert_that(is.unsorted(rev(round(lower_slopes, 10)), FALSE) == FALSE)

        # update upper areas
        upper_areas <- update_areas(n, position, upper_x, upper_slopes, upper_ints, upper_areas, bounds)

        # update number of points in lower polygon
        n <- n + 1
      }
    }
  }
  return(accepted_points)
}
