
#######################
# Demo from Baptiste Alglave
#######################

## Construct the seasonnal variable
p <- ncol(S) # number of time steps
ts <- seq(1,p,1) # sequence of time steps
component.strength = 1 # amplitude of the signal
component.freqs = 1/12 # frequency of signal (month^-1)
f.0 = 2 * pi # fundamental frequency (month^-1)
component.delay = 1 # delay of signal components (radians)

# Create the seasonal signal
signal_vec = - component.strength *
sin(component.freqs * f.0 * ts + component.delay)

## Filter EOFs dimensions and apply CCA
ndim <- 2 # number of dimension to filter in the EOFs
EOFset1 <- E$v[,1:ndim] # filter the loadings
cca_res = stats::cancor(EOFset1, signal_vec) # perform CCA

## Rotate the temporal loadings and spatial factors
## with CCA results
# Rotate the temporal loadings of the EOFs
ccavar1_t <- (as.matrix(EOFset1) %*% cca_res$xcoef)[,1]
# Rotate the spatial factor of the EOFs
ccavar1_x <- as.matrix(E$u[,1:ndim) %*% cca_res$xcoef[,1]


#####################
# Modified for VAST using a model with four quarters per year
#####################

# Assuming fit is output from `fit_model`

# Loop through species c \in {1,2,...C}
c_index = 1
E = list(
          "v" = fit$Report$Ltime_epsilon1_tf,
          "u" = fit$Report$Epsiloninput1_gff[,c_index,]
)

## Construct the seasonnal variable
p <- ncol(S) # number of time steps
ts <- seq(1,p,1) # sequence of time steps
component.strength = 1 # amplitude of the signal
component.freqs = 1/4 # frequency of signal (quarter^-1)
f.0 = 2 * pi # fundamental frequency (month^-1)
component.delay = 1 # delay of signal components (radians)

# Create the seasonal signal
signal_vec = - component.strength *
sin(component.freqs * f.0 * ts + component.delay)

## Filter EOFs dimensions and apply CCA
ndim <- 2 # number of dimension to filter in the EOFs
EOFset1 <- E$v[,1:ndim] # filter the loadings
cca_res = stats::cancor(EOFset1, signal_vec) # perform CCA

## Rotate the temporal loadings and spatial factors
## with CCA results
# Rotate the temporal loadings of the EOFs
ccavar1_t <- (as.matrix(EOFset1) %*% cca_res$xcoef)[,1]
# Rotate the spatial factor of the EOFs
ccavar1_x <- as.matrix(E$u[,1:ndim]) %*% cca_res$xcoef[,1]

