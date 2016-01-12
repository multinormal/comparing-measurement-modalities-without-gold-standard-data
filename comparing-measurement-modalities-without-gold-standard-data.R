# This code and simulation are based on Hoppin et al., IEEE TMI 21(5) 2002.
# The file is in four parts. The first sets up the R environment.
# The second sets up the environment further to replicates the first experiment
# in the Hoppin et al. paper, but using JAGS. The simulation uses invented gold
# standard data that would of course not be available if using this method for a
# real analysis, since assessing accuracy and precision of a number of modalities
# in the absence of gold standard data is the problem being addressed! The third
# part applies the method to the invented data. This part could be reused to
# actually perform such an analysis on real data. Note that you would need to edit
# the JAGS model file, in particular the specification of the distribution of the
# quantity being estimated. The fourth part compares the modalities in terms of
# the slope, intercept, and standard deviation parameters (defined by Hoppin et al.)
# by computing the posterior probabilities that the magnitude of a parameter for
# a given modality is lower than the magnitude of a parameter for another modality.
# This fourth part is specific to the replication of the Hoppin et al. analysis
# since it assumes three modalities; it also assumes the use of two chains in the
# MCMC run, and as such it is not sufficiently generic to run as-is if these
# assumptions are violated. However, it should give you an idea as to how to code
# such an analysis.

#############################################
# PART 1: Set up the environment
#############################################

library(runjags) # You need the JAGS binary installed, the runjags package, and its dependencies.

# Set the random number generator seed, so we (should) get repeatable results each run.
set.seed(1234)

# Set the working directory to the directory containing this .R file, so that
# the path to the JAGS model file can be specified simply.
current.directory <- dirname(sys.frame(1)$ofile)
setwd(current.directory)

#############################################
# PART 2: Set up the Hoppin et al. simulation
#############################################

# Specify the number of observations (e.g., patients).
n.obs <- 100

# Specify and sample from the distribution of the gold standard data.
gold.standard.data <- rbeta(n = n.obs, shape1 = 1.5, shape2 = 2.0) # Values as for Exp A in Hoppin et al.

# The relationships between the estimates made by the modalities and the gold standard
# data are modeled as linear with normally-distributed residuals. Specify the slopes,
# intercepts, and standard deviations on the distributions of the residuals for the simulation.
# These are the parameters that are of interest and will be estimated. Compare the result
# of the analysis to these values.
a.1 <- 0.6
a.2 <- 0.7
a.3 <- 0.8 # The slopes.
b.1 <- -0.1
b.2 <- 0.0
b.3 <- 0.1 # The intercepts.
s.1 <- 0.05
s.2 <- 0.03
s.3 <- 0.08 # The standard deviations.

# Generate the estimates for the three modalities.
estimates.1 <- (a.1 * gold.standard.data) + b.1 + rnorm(n = n.obs, mean = 0, sd = s.1)
estimates.2 <- (a.2 * gold.standard.data) + b.2 + rnorm(n = n.obs, mean = 0, sd = s.2)
estimates.3 <- (a.3 * gold.standard.data) + b.3 + rnorm(n = n.obs, mean = 0, sd = s.3)

# Plot the gold standard data against the estimates.
op <- par(mfrow = c(1,3))
plot(gold.standard.data, estimates.1)
plot(gold.standard.data, estimates.2)
plot(gold.standard.data, estimates.3)
par(op)

# Bundle the estimates into a matrix; this would be the raw data in an actual analysis,
# but here we use the simulated values. Note that we do not include the simulated gold
# standard data in this matrix, as that must not be used in the analysis.
# The rows of observations matrix are the individuals (e.g., patients) and the columns
# are the modalities.
observations <- matrix(c(estimates.1, estimates.2, estimates.3), nrow = n.obs)

#############################################
# PART 3: Perform the analysis
#############################################

# A real analysis would start at this point, with an n.obs x n.mods matrix called observations
# containing the n.obs estimates of the gold standard values made using the n.mods modalities.

# Specify the variables that we want to estimate; these must have the same names as the corresponding
# variables in the JAGS file.
model.parameters <- c('a', 'b', 's') # This must be a vector (i.e., made using c()).
  # Note that a, b, and s are vectors, each with the same number of elements as there are modalities.

# Specify the path to the JAGS file (here, relative to the current directory).
model.file <- 'comparing-measurement-modalities-without-gold-standard-data.jags'

# Specify initialization parameters for JAGS. In particular, specify the pseudo-random
# number generators and their seeds, which should guarantee reproducible analyses.
# Note that this also implicitly specifies the number of chains.
inits <- list(
              list(.RNG.name = 'base::Super-Duper', .RNG.seed = 12345),
              list(.RNG.name = 'base::Mersenne-Twister', .RNG.seed = 67890)
              )

# Run JAGS.
result <- run.jags(data = list(observations = observations),
                      model = model.file,
                      monitor = model.parameters,
                      thin = 2,
                      inits = inits,
                      modules = 'glm'
                      )

# Print a tabular summary of the estimates.
print(result)

# Plot summaries of the MCMC sampling and posterior distributions for the estimated parameters.
plot(result)


#############################################
# PART 4: Posterior comparative probabilities
#############################################

# NOTE: This assumes you're using two chains in run.jags, which is the default.
mcmc.samples <- rbind(result$mcmc[,][[1]], result$mcmc[,][[2]])

# Define a function that estimates the posterior probability that the
# value of a linear model parameter is closer to the ideal value of
# that parameter, w, than another linear model parameter. This allows
# us to calculate quantities such as P(|a[1] - w| < |a[2] - w|).
probability.left.better.than.right <- function(left.colname, right.colname, w = 0) {
  samples <- mcmc.samples[, c(left.colname, right.colname)]
  stopifnot(nrow(mcmc.samples) == nrow(samples))
  n.better <- length(which( abs(samples[,1] - w) < abs(samples[,2] - w) ))
  n.better / nrow(mcmc.samples)
}

# Calculate the posterior probabilities.
p.a1.better.than.a2 <- probability.left.better.than.right('a[1]', 'a[2]', w = 1)
p.a1.better.than.a3 <- probability.left.better.than.right('a[1]', 'a[3]', w = 1)
p.a2.better.than.a3 <- probability.left.better.than.right('a[2]', 'a[3]', w = 1)

p.b1.better.than.b2 <- probability.left.better.than.right('b[1]', 'b[2]')
p.b1.better.than.b3 <- probability.left.better.than.right('b[1]', 'b[3]')
p.b2.better.than.b3 <- probability.left.better.than.right('b[2]', 'b[3]')

p.s1.better.than.s2 <- probability.left.better.than.right('s[1]', 's[2]')
p.s1.better.than.s3 <- probability.left.better.than.right('s[1]', 's[3]')
p.s2.better.than.s3 <- probability.left.better.than.right('s[2]', 's[3]')

# Print the posterior probabilities.

print('Posterior inferences for the slope parameter:')
print(paste('Posterior P(|a[1] - 1| < |a[2] - 1|) = ', p.a1.better.than.a2, sep = ' '))
print(paste('Posterior P(|a[1] - 1| < |a[3] - 1|) = ', p.a1.better.than.a3, sep = ' '))
print(paste('Posterior P(|a[2] - 1| < |a[3] - 1|) = ', p.a2.better.than.a3, sep = ' '))

print('Posterior inferences for the intercept parameter:')
print(paste('Posterior P(|b[1]| < |b[2]|) = ', p.b1.better.than.b2, sep = ' '))
print(paste('Posterior P(|b[1]| < |b[3]|) = ', p.b1.better.than.b3, sep = ' '))
print(paste('Posterior P(|b[2]| < |b[3]|) = ', p.b2.better.than.b3, sep = ' '))

print('Posterior inferences for the standard deviation parameter:')
print(paste('Posterior P(|s[1]| < |s[2]|) = ', p.s1.better.than.s2, sep = ' '))
print(paste('Posterior P(|s[1]| < |s[3]|) = ', p.s1.better.than.s3, sep = ' '))
print(paste('Posterior P(|s[2]| < |s[3]|) = ', p.s2.better.than.s3, sep = ' '))