## The functions used in the hybrid zone analysis

# Summarize phenotype counts ----------------------------------------------
sumPhenoCounts <- function(df, factor_vec) {

  # Take in the data frame of individual capture records,
  # and summarize the number of phenotypes captured.
  # Can summarize by site, subsite, or both: simply provide a 
  # vector of quoted column names to summarize over,
  # E.g.: factor_vec = "subsite"  or 
  # factor_vec = c("subsite", "site.collected")
  
  # Make sure the necessary columns are present
  assert_that("A.melanized" %in% colnames(df), msg = "Input data must contain a column named \"A.melanized\"")
  assert_that("B.hetero" %in% colnames(df), msg = "Input data must contain a column named \"B.hetero\"")
  assert_that("C.west.col" %in% colnames(df), msg = "Input data must contain a column named \"C.west.col\"")
  assert_that("D.postman" %in% colnames(df), msg = "Input data must contain a column named \"D.postman\"")
  # Make sure the grouping columns are all present
  assert_that(min(unique(factor_vec %in% colnames(df)))== 1, msg = "The columns given in factor_vec must be present in the input data frame")
  # do the actual counting of phenotypes
  phenoCounts<- df %>%
    group_by_(.dots = factor_vec) %>%
    summarize(A.melanized = sum(A.melanized, na.rm = T),
              B.hetero = sum(B.hetero, na.rm = T),
              C.west.col = sum(C.west.col, na.rm = T),
              D.postman = sum(D.postman, na.rm = T)) %>%
    ungroup()
  return(phenoCounts)
}


# Convert pheno counts to allele counts -----------------------------------
phenoToAlleleCount <- function(df, remove = T) {
  # Takes a data frame with numbers of individuals in the 4 phenotype classes
  # for each site/population, and converts those to counts of alleles.
  # Default behavior is to lump the west-colombian yellow hindwing allele
  # with the central-american yellow hindwing allele
  
  # Make sure the necessary columns are present
  assert_that("A.melanized" %in% colnames(df), msg = "Input data must contain a column named \"A.melanized\"")
  assert_that("B.hetero" %in% colnames(df), msg = "Input data must contain a column named \"B.hetero\"")
  assert_that("C.west.col" %in% colnames(df), msg = "Input data must contain a column named \"C.west.col\"")
  assert_that("D.postman" %in% colnames(df), msg = "Input data must contain a column named \"D.postman\"")
  
  
  #First, make allele columns for allele counts
  alleleCounts <- df %>%
    mutate(HYDallele = 2*A.melanized + B.hetero,
           YELallele = 2*D.postman + B.hetero + 2*C.west.col) %>%
    mutate(totalAlleles = HYDallele + YELallele)
  
  #Now, some checks
  if(unique(alleleCounts$HYDallele + alleleCounts$YELallele == alleleCounts$totalAlleles) != T) {
    warning("Allele counts don't add up properly")
  }
  if(unique((df$A.melanized + df$B.hetero + df$C.west.col + df$D.postman)*2 == alleleCounts$totalAlleles) != T) {
    warning("Mismatch between # of individuals and # of alleles")
  }
  
  # if we pass the checks, remove the phenotype columns if remove == T
  if (remove == T) {
    alleleCounts <- alleleCounts %>%
    dplyr::select(-A.melanized, -B.hetero, -C.west.col, -D.postman)
  return(alleleCounts)
  }
  else {
    return(alleleCounts)
  }
}




# Do all pairwise fisher’s exact test on a df of subsites -----------------
pairwiseFisherInSite <- function(df) {
  # Takes in a data frame with allele counts for each of the subsets within a site
  # Does a fisher's exact test to compare allele counts for each possible pairwise comparison
  # returns those p-values as a data frame
  subsites <- dim(df)[1]
  #determine # of pairwise comparisons
  pw_comps <- (subsites * (subsites-1))/2 
  # set up a results data frame, so this will play nice with dplyr
  results <- data.frame(subsite1 = rep("NA", times = pw_comps),
                        subsite2 = rep("NA", times = pw_comps),
                        p.value = rep(99, times = pw_comps), stringsAsFactors=F)
  # do the pairwise comparisons across the data frame
  results_row <- 1
  for (row1 in 1:subsites) {
    for (row2 in 1:subsites) {
      if (row1 == row2) {
        next # skip iteration if comparing same row
      } 
      else if (row1 > row2) {
        next # skip iteration if row1 > row2, so we only do one direction of each pairwise comparison
      }
      else {
        results[results_row, 1] <- df$subsite[row1]
        results[results_row, 2] <- df$subsite[row2]
        # do the fisher exact test
        testVal <- fisher.test(matrix(c(df$HYDallele[row1], df$YELallele[row1],
                                        df$HYDallele[row2], df$YELallele[row2]),
                                      nrow = 2, byrow = T))$p.value
        results[results_row, 3] <- testVal
        
        # iterate the results row
        results_row <- results_row + 1
      }
    }
  }
  return(results)
}


# Test subsite differentiation for all sites ------------------------------
testSubsiteDiff <- function(df) {
  ## A function to test differentiation bewtween a (variable) number of sites. 
  ## Does pairwise Fisher's exact tests for all subsites within a site.
  ## Can't do a chi-square contingency test: some cells have 0 or <5.
  ## Relies on the custom pairwiseFisherInSite funciton to do the pairwise
  ## tests for each site
  results <- df %>%
    group_by(site.collected) %>%
    do(pairwiseFisherInSite(.)) %>%
    ungroup()
  return(results)
}



# Add coordinates for point on transect -----------------------------------
addPointOnTransect <- function(df) {
  # this function reads in a dataframe, which must contain columns:
  # coord.W.decdeg (GPS coords in decimal degrees)
  # coord.N.decdeg.
  # For each row in the data frame, if calls the findClosestPointonCubic
  # function, which finds the closest point on the fitted transect.
  # It then adds new columns giving the site's new coordinates on the cubic transect
  for (row in 1:dim(df)[1]) {
    point_on_transect <- findClosestPointonCubic(df[row,])
    df$tran.Coord.N[row] <- as.numeric(point_on_transect[1])
    df$tran.Coord.W[row] <- as.numeric(point_on_transect[2]) 
  }
  return(df)
}


# Find closest point on transect ------------------------------------------
findClosestPointonCubic <- function(df_row) {
  # This function takes in a row of a data frame containing the N and W decimal degree coordinates
  # for a collection site. 
  # It then examines the range of points on the cubic transect that are within 0.5 degrees longitude
  # (30 minutes) of the site, and calculates the distance between the each of those points on the transect,
  # and the focal site. This distance calculation accounts for the curvature of the earth. 
  # it then returns a vector of GPS coordinates for the closest site on the transect
  results <- rep(NA, 30001)
  site.p1 <- c(df_row$coord.W.decdeg, df_row$coord.N.decdeg)
  # don't test all points, just the ones within 0.5 degree longitude of our focal site
  range <- which(cubic.trans$coord.W.decdeg < (site.p1[1] + 0.5) & cubic.trans$coord.W.decdeg > (site.p1[1] - 0.5))
  for (point in range) {
    results[point] <- distVincentyEllipsoid(p1 = site.p1, p2 = c(cubic.trans$coord.W.decdeg[point],
                                                                 cubic.trans$coord.N.decdeg[point]))
  }
  closest <- cubic.trans[which(results == min(results, na.rm = T)),]
  return(as.vector(closest))
}


# Add transect distance columns -----------------------------------------------------
addDistanceColsArc <- function(site_data, coefficients) {
  # This function takes in 1) a data frame, which must contain the GPS data
  # for the site on the transect in columns named tran.Coord.N and
  # tran.Coord.W.
  # and 2) a vector of coefficients for the cubic function describing the transect
  #
  # It adds two columns: dist.to.last, which contains the distance, in km,
  # from a site to the site before it.
  # and transect.dist, which contains the distance along the transect
  # This distance is calculated as the distance along the cubic transect between
  # the two points, accounting for curvature of the earth.
  

  # First, tests
  
  # Make sure the input data frame is sorted properly.
  test <- unique(arrange(site_data, tran.Coord.W)$tran.Coord.W == site_data$tran.Coord.W)
  if (F %in% test) {# sort it if it isn't
    site_data <- arrange(site_data, tran.Coord.W)
  }
  # Make sure the coefficients are numeric and of length 4
  if (length(coefficients) != 4) {
    stop("coefficients not of proper length for a cubic function")
  }
  if (is.numeric(coefficients) != T) {
    stop("coefficients not of type numeric")
  }
  
  # Then, add the distance to last column
  for (row in 1:dim(site_data)[1]) {
    if (row == 1) {# is this the first site?
      # if so, dist.to.last = 0
      site_data$dist.to.last[row] <- 0
    } else { # for all other sites, set p2 and p1 as the current row and row before it
      # then calculate the distance between the points
      low <- site_data$tran.Coord.W[row-1]
      up <- site_data$tran.Coord.W[row]
      site_data$dist.to.last[row] <- arc_length(coefficients, low, up)
    }
  }
  ## Then, create the transect.dist column
  site_data <- mutate(site_data, transect.dist = cumsum(dist.to.last))
  return(site_data)
}

# Calculus function ------------------------------------------------------
# Written by Joe Thurman
arc_length <- function(coefficients, lower, upper) {
  ## Arc length on cubic between points lower and upper, accounting for Earth's
  ## curvature
  ##
  ## Args:
  ## coefficients - numerical vector of length 4, giving the coefficients of a cubic polynomial in ascending degree (first entry is constant term, last entry is coefficient on cubic)
  ## lower - number, longitude of left point on cubic
  ## upper - number, longitude of right poing on cubic
  ##
  ## Returns:
  ## Arc length, in km, of cubic described by coefficients between lower and upper.
  
  R <- 6378.1370 #Radius of earth, in km, at equator
  
  cubic <- function(x) {coefficients[1] + coefficients[2]*x + coefficients[3] * x^2 + coefficients[4]*x^3}
  derivative <- function(x) {coefficients[2] + 2*coefficients[3]*x + 3* coefficients[4] * x^2}
  integrand <- function(x) {sqrt((cos(cubic(x)*pi/180))^2 + (derivative(x))^2)}
  
 return(R* pi/180 * integrate(integrand, lower, upper)$value)
}



# Generate points along a transect at a specified interval -----------------
# Given a cubic transect and a starting point, return 
# a data frame containing points along that transect that are a specified
# distance from each other (interval), and that extend to the desired point
# (end.dist) on the transect. 

gen_points_on_trans <- function(transect.coeffs, start.long, end.dist, interval) {

  # Given the endpoint of the transect and the desired interval between points
  # Find the number of rows needed to make it to the end (go 1 past if necessary)
  # and generate a results dataframe
  res.rows <- ceiling(end.dist/interval)
  
  # Internally define the cubic function that defines our transect. 
  # Cubic function to describe the transect. 
  cubic <- function(x) {transect.coeffs[1] + transect.coeffs[2]*x + transect.coeffs[3] * x^2 + transect.coeffs[4]*x^3}
  
  # internal function that returns the absolute value of the difference
  # between the calculated arc length and the desired arc length.
  # Will use this in R's optimizer to find the minimum, with the first 
  # argument (upper limit) the one to be optimized. 
  arc_diff <- function(upper, low, target.dist, coefficients) {
    z <- arc_length(coefficients = coefficients, lower = low, upper = upper)
    result <- abs(target.dist - z)
  }
  
  # Make results data frame outside the loop
  results <- data.frame(tran.Coord.N = rep(NA, times = res.rows), 
                        tran.Coord.W = rep(NA, times = res.rows),
                        dist.to.last = rep(NA, times = res.rows),
                        transect.dist = rep(NA, times = res.rows))
  
  # Fill it with results
  for (row in 1:dim(results)[1]) {
    if (row == 1) { # in first row, want to fill in the coordinates with the starting point
      results$tran.Coord.W[row] <- start.long
      results$tran.Coord.N[row] <- cubic(results$tran.Coord.W[row]) 
      results$dist.to.last[row] <- 0
      results$transect.dist[row] <- 0
    }
    else {
      results$tran.Coord.W[row] <- optimize(arc_diff, 
                                            interval = c(results$tran.Coord.W[row-1], results$tran.Coord.W[row-1] + 1), 
                                            maximum = F, 
                                            tol = 0.000001, 
                                            target.dist = interval,
                                            coefficients = transect.coeffs,
                                            low = results$tran.Coord.W[row-1])$minimum
      results$tran.Coord.N[row] <- cubic(results$tran.Coord.W[row])
      results$dist.to.last[row] <- round(arc_length(coefficients = transect.coeffs, 
                                                    lower = results$tran.Coord.W[row-1], 
                                                    upper = results$tran.Coord.W[row]), digits = 3)
      results$transect.dist[row] <- round(arc_length(coefficients = transect.coeffs, 
                                                     lower = results$tran.Coord.W[1], 
                                                     upper = results$tran.Coord.W[row]), digits = 3)
    }
  }
  return(results)
}


# correct FIS, for effective samples in hzar ------------------------------
# a quick little function to turn NaNs and negative numbers into 0. 
correct_Fis <- function(df) {
  for (element in 1:length(df$Fis)) {
    if (is.nan(df$Fis[element]) == T) {
      df$Fis[element] <- 0
    }
    else if (df$Fis[element] < 0) {
      df$Fis[element] <- 0
    }
  }
  return(df)
}


# General cline equation --------------------------------------------------
# A flexible function that generates values from clines

general_cline_eqn <- function(transectDist, decrease,
                              center, width,
                              pmin = 0, pmax = 1,
                              deltaL = NULL, tauL = NULL,
                              deltaR = NULL, tauR = NULL) {
  
  # Start with an ungodly amount of argument checking
  # check transect distance is a numeric vector
  assertthat::assert_that(is.vector(transectDist) == T, msg = "transect_distances must be a vector")
  assertthat::assert_that(is.numeric(transectDist) == T, msg = "transect_distances must be of type numeric")
  
  # Check decrease is T/F
  assertthat::assert_that(is.logical(decrease) == T, msg = "decrease must be either TRUE (T) or FALSE (F)")
  
  # Check the cline parameters are numeric vectors of length 1
  for (num.arg in alist(center, width, pmin, pmax, deltaL, deltaR, tauL, tauR)) {
    if (is.null(eval(num.arg)) == F) {
      assertthat::assert_that(is.vector(eval(num.arg)) == T, msg = paste(num.arg, "must be a vector", sep = " "))
      assertthat::assert_that(is.numeric(eval(num.arg)) == T, msg = paste(num.arg, "must be numeric", sep = " "))
      assertthat::assert_that(length(eval(num.arg)) == 1, msg = paste(num.arg, "must be of length 1", sep = " "))
    }
  }
  # Center and width must be greater than 0
  assertthat::assert_that(center >= 0, msg = "center must be greater than 0")
  assertthat::assert_that(width >= 0, msg = "width must be greater than 0")
  
  # Pmin, pmax, tauL, and tauR must be between 0 and 1
  for (num.arg in alist(pmin, pmax, tauL, tauR)) {
    if (is.null(eval(num.arg)) == F) {
      assertthat::assert_that(eval(num.arg) >= 0, msg = paste(num.arg, " must be between 0 and 1 (inclusive)", sep = ""))
      assertthat::assert_that(eval(num.arg) <= 1, msg = paste(num.arg, " must be between 0 and 1 (inclusive)", sep = ""))
    }
  }
  
  # Check to make sure both delta and tau are supplied for a given set of tails
  assertthat::assert_that(sum(is.null(deltaL), is.null(tauL)) %in% c(0,2),
                          msg = "If using deltaL or tauL, must supply both of them")
  assertthat::assert_that(sum(is.null(deltaR), is.null(tauR)) %in% c(0,2),
                          msg = "If using deltaR or tauR, must supply both of them")
  
  # The cline equations are written for increasing clines
  # With two alleles, a decreasing cline is simply
  # 1- increasing
  # Will do that at the end, after computng p.
  
  # Calculate p, using the proper equation based on the parameters prodived
  if (sum(is.null(deltaL), is.null(deltaR)) == 2) { # if no tails
    prop <- (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
  }
  else if (sum(is.null(deltaL), is.null(deltaR)) == 0) {# If both tails
    if (transectDist <= center - deltaL) { # and we're in the left tail
      # use the left tail equation
      prop <- (1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist - center + deltaL)/width)/(1 + exp(-4*deltaL/width)))
    }
    else if (transectDist >= center + deltaR) { # and we're in the right tail
      # use the right tail equation
      prop <- (1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist - center - deltaR)/width)/(1 + exp(-4*deltaR/width))))
    }
    else {# Otherwise, we're in the sigmoid center
      prop <- (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
    }
  }
  else if (is.null(deltaL) == F) { # If left tail only
    if (transectDist <= center - deltaL) { # and we're in the left tail
      # use the left tail equation
      prop <- (1/(1 + exp(4*deltaL/width)))*exp((4*tauL*(transectDist - center + deltaL)/width)/(1 + exp(-4*deltaL/width)))
    }
    else {# Otherwise, we're in the sigmoid center
      prop <- (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
    }
  }
  else { # must be right tail cline
    if(transectDist >= center + deltaR) { # and we're in the right tail
      # use the right tail equation
      prop <- (1-(1/(1 + exp(4*deltaR/width)))*exp((-4*tauR*(transectDist - center - deltaR)/width)/(1 + exp(-4*deltaR/width))))
    }
    else {# Otherwise, we're in the sigmoid center or there are no tails
      prop <- (exp(4*(transectDist - center)/width)/(1 + exp(4 * (transectDist - center)/width)))
    }
  }
  
  # The exponential function in R will got to Inf once it hits 710.
  # So, if (distance-center)/width gets to be 177.5, the exp() function goes to Inf
  # This turns prop into NaN, when it should be 1.
  # Will just add a check for this.
  if (is.nan(prop) == T) {
    prop <- 1
  }
  
  
  if (decrease == T) {
    p <- pmin + (pmax - pmin)*(1 - prop)
  }
  if (decrease == F) {
    p <- pmin + (pmax - pmin)*prop
  }
  p
}


# Sim data from cline -----------------------------------------------------
# Simulate data from a cline

sim_data_from_cline <- function(transect_distances, n_ind,
                                Fis, decrease,
                                center, width,
                                pmin = 0, pmax = 1,
                                deltaL = NULL, tauL = NULL,
                                deltaR = NULL, tauR = NULL) {
  
  
  # Check the sampling and inbreeding options
  for (vec.arg in alist(n_ind, Fis)) {
    assertthat::assert_that(is.vector(eval(vec.arg)) == T,
                            msg = paste(vec.arg, "must be a vector", sep = " "))
    assertthat::assert_that(is.numeric(eval(vec.arg)) == T,
                            msg = paste(vec.arg, "must be numeric", sep = " "))
    assertthat::assert_that((length(eval(vec.arg)) %in% c(1, length(transect_distances))) == T,
                            msg = paste(vec.arg, " must be either be of length 1, for constant ", vec.arg,
                                        ", or must match the length of transect_distances (",
                                        length(transect_distances), sep = ""))
  }
  assertthat::assert_that(min(Fis) >=0, msg = "Fis values cannot be less than 0")
  assertthat::assert_that(min(Fis) <=1, msg = "Fis values cannot be greater than 1")
  # All other args will get checked in the cline equation.
  
  # Get number of sites from the vector of transect data.
  sites <- length(transect_distances)
  # Get the vector of f values for each site
  if (length(n_ind) == 1) {
    Ns <- as.integer(rep(n_ind, times = sites))
  } else {
    Ns <- as.integer(n_ind)
  }
  if (length(Fis) == 1) {
    fs <- rep(Fis, times = sites)
  } else {
    fs <- Fis
  }
  
  # Make the empty results data frame
  fk.dt <- data.frame(site = 1:sites,
                      transectDist = transect_distances,
                      cline.p = rep(NA, times = sites),
                      cline.f = fs,
                      AA = rep(NA, times = sites),
                      Aa = rep(NA, times = sites),
                      aa = rep(NA, times = sites),
                      N = Ns)
  
  # Then add the simulated genotypes to each row
  for (row in 1:sites) {
    fk.dt$cline.p[row] <- general_cline_eqn(transectDist = fk.dt$transectDist[row], center = center, width = width,
                                            pmin = pmin, pmax = pmax, deltaL = deltaL, deltaR = deltaR, tauL = tauL,
                                            tauR = tauR, decrease = decrease)
    AA <- fk.dt$cline.p[row]^2 + fk.dt$cline.f[row]*fk.dt$cline.p[row]*(1-fk.dt$cline.p[row])
    Aa<- 2*fk.dt$cline.p[row]*(1-fk.dt$cline.p[row])*(1-fk.dt$cline.f[row])
    aa <- (1-fk.dt$cline.p[row])^2 +fk.dt$cline.f[row]*fk.dt$cline.p[row]*(1-fk.dt$cline.p[row])
    genotypes <- rowSums(stats::rmultinom(n = 1, size = fk.dt$N[row], prob = c(AA, Aa, aa)))
    fk.dt$AA[row] <- as.integer(genotypes[1])
    fk.dt$Aa[row] <- as.integer(genotypes[2])
    fk.dt$aa[row] <- as.integer(genotypes[3])
  }
  

  # Calculate empirical p and Fis value from the simulated data
  fk.dt <- fk.dt %>%
    dplyr::mutate(emp.p = (2*.data$AA + .data$Aa)/(2*.data$N)) %>%
    dplyr::mutate(Hexp = 2*.data$emp.p*(1-.data$emp.p),
                  Hobs = .data$Aa/.data$N) %>%
    dplyr::mutate(Fis = (.data$Hexp - .data$Hobs)/.data$Hexp) %>%
    correct_fis(.) %>%
    dplyr::rename(emp.f = .data$Fis) %>%
    dplyr::select(-.data$Hexp, -.data$Hobs)
  
  # Do some rounding
  fk.dt$cline.p <- round(fk.dt$cline.p, digits = 3)
  fk.dt$emp.p <- round(fk.dt$emp.p, digits = 3)
  fk.dt$emp.f <- round(fk.dt$emp.f, digits = 3)
  fk.dt
}

# Summarize cline data ----------------------------------------------------



cline_summary <- function(stanfit, prob = .95, method = "HPDI", show.all = F) {
  assertthat::assert_that(class(stanfit)[1] == "stanfit",
                          msg = "Object to be summarized must be of class stanfit")
  assertthat::assert_that(is.numeric(prob) == T, msg = "prob must be numeric")
  assertthat::assert_that(length(prob) == 1, msg = "prob must be of length 1")
  assertthat::assert_that(prob <= 1, msg = "prob must be between 0 and 1")
  assertthat::assert_that(prob > 0, msg = "prob must be between 0 and 1")
  assertthat::assert_that((method %in% c("HPDI", "ET")) == T,
                          msg = "method must be either 'HPDI' or 'ET")
  assertthat::assert_that(is.logical(show.all) == T, msg = "show.all must be either TRUE or FALSE")
  
  # Give a method to do HPDI vs. eqaul-tail interval
  tail <- (1-prob)/2
  low.name <- paste("low", prob, method, sep = "_")
  up.name <- paste("up", prob, method, sep = "_")
  
  # could add [abcdeghijklmnopqrstuvwxyz] in the reg expression below to
  # also keep the column with inbreeding values
  if (show.all == F) {
    keep <- grep("\\[|_", names(stanfit), invert = T, value = T)
  } else {
    keep <- names(stanfit)
  }
  
  res <- rstan::summary(stanfit, probs = c(0 + tail, 1 - tail), pars = keep)$summary %>%
    as.data.frame(.) %>%
    round(., digits = 4) %>%
    dplyr::mutate(n_eff = as.integer(.data$n_eff)) %>%
    cbind(keep, .)
  
  if (method == "HPDI") {
    hpd_cols <- as.matrix(stanfit, pars = keep) %>%
      coda::as.mcmc(.) %>%
      coda::HPDinterval(obj = ., prob = prob) %>%
      round(., digits = 4) %>%
      as.matrix(.)
    res[,5:6] <- hpd_cols
  }
  
  names(res)[c(1,5,6)] <- c("param", low.name, up.name)
  
  res
}

