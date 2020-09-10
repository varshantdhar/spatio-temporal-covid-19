covid.states <- read.csv(url("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv"))
covid.counties <- read.csv(url("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv"))
states <- unique(covid.states$state)
counties <- unique(covid.counties$county)

dates <- unique(covid.states$date)
no_dates <- nlevels(dates)
data.states <- matrix(0, nrow=no_dates,ncol=nlevels(states))
colnames(data.states) <- states
row.names(data.states) <- unique(covid.states$date)

for (row in dates){
  data.states[row,c(as.character(covid.states$state[covid.states$date == row]))] <- covid.states$cases[covid.states$date == row]
}
# Cumulative Sums
dt <- matrix(0, nrow=no_dates-1,ncol=nlevels(states))
for (row in c(dates)[-no_dates]){
  dt[row,] <- data.states[row+1,] - data.states[row,]
}
colnames(dt) <- colnames(data.states)
dt <- rbind(data.states[1,],dt)
row.names(dt) <- row.names(data.states)
dt = dt[,-56]
# Data Cleaning complete
```

```{r}
N = dim(dt)[1] # number of observations
p = dim(dt)[2] # number of regions
# In this model:
# y_i = \Theta\xi(x_i)\eta_i + \epsilon_i
# \eta_i = \psi(x_i) + \xi_i
# \epsilon_i \sim N(0,\Sigma_0),    \xi_i \sim N(0,I).
x = c(1:N)/N

c = 100  # Gaussian process length scale parameter.  Smaller c will lead to slower changes in correlations.
d = 1
r = 1e-5
K = matrix(0, nrow=N, ncol=N)
for(ii in c(1:N)){
  for(jj in c(1:N)){
    dist_ii_jj = abs(x[ii]-x[jj])
    K[ii,jj] = d*exp(-c*(dist_ii_jj^2))
  }
}
K = K + diag(r*rep(1,times=N),nrow = N, ncol = N)
invK = solve(K)
logdetK = 2*sum(log(diag(chol(K))))

prior_params = c(K.c_prior = 1, K.logdetK = logdetK, sig.a_sig = 1, sig.b_sig = 0.1, hypers.a_phi = 1.5, hypers.b_phi = 1.5, hypers.a1 = 10, hypers.a2 = 10); 

settings = c(
  L = 10, # truncation level for dictionary of latent GPs
  k = 20, # latent factor dimension
  Niter = 10000, # number of Gibbs iterations to run
  sample_K_flag = 3, # choose to fix bandwidth parameter instead of sampling it
  latent_mean = 1, # indicate whether one wants to model a non-zero latent mean \mu(x)
  restart = 0)

y = t(dt)
# invK, K, inds_y are important function parameters
inds_y = matrix(1, nrow = p, ncol = N)
inds_y[y==0] = 0
inds_y = inds_y > 0 

BNP_covreg(y, prior_params, settings, invK, K, inds_y)

BNP_covreg <- function(y, prior_params, settings, invK, K, inds_y){
  p = dim(y)[1] # number of regions
  N = dim(y)[2] # number of observations
  sample_K_flag = settings["sample_K_flag"]
  latent_mean = settings["latent_mean"]
  k = settings["k"]
  L = settings["L"]
  Niter = settings["Niter"]
  
  # Sample hyperparams from prior:
  delta = c(rep(0, times=L))
  delta[1] = rgamma(1, shape=prior_params["hypers.a1"], scale = 1)
  delta[2:L] = rgamma(length(c(2:L)),shape=prior_params["hypers.a2"], scale=1)
  tau = exp(cumsum(log(delta)))
  phi = matrix(rgamma(n=p*L,shape=prior_params["hypers.a_phi"],scale=1),nrow = p, ncol = L)
  
  # Sample theta, eta, and Sigma initially as prior draws:
  theta = matrix(0, nrow=p, ncol=L)
  for(pp in c(1:p)){
    theta[pp,] = t(chol(diag(1/(phi[pp,]*tau))))%*%matrix(rnorm(L),nrow=L,ncol=1)
  }
  
  xi = matrix(rnorm(k*N),nrow=k,ncol=N)
  psi = matrix(0, nrow=k, ncol=N)
  eta = psi + xi
  
  # invSig_vec represents the diagonal elements of \Sigma_0^{-1}:
  invSig_vec = c(rgamma(n=p,shape = prior_params["sig.a_sig"],scale=1)/prior_params["sig.b_sig"])
  invSig = solve(y%*%t(y))
  
  # Sample initial GP cov K and GP latent functions zeta_i:
  K_ind = 1
  
  # Sample zeta_i using initialization scheme based on data and other
  # sampled params:
  zeta = array(rep(0,times = L*k*N), c(L,k,N))
  zeta = sample_zeta(y,theta,eta,invSig_vec,zeta,invK,inds_y);
  
  num_iters = 1
  nstart = 1
  
  for (nn in c(nstart:Niter)) {
    break
  }
}

sample_zeta <- function(y, theta, eta, invSig,zeta,invK,inds_y){
  p = dim(y)[1]
  N = dim(y)[2]
  L = dim(theta)[2]
  k = dim(eta)[1]
  
  AinvSig = matrix(0, nrow = N*k*L, ncol = N*p)
  AinvSigA = matrix(0, nrow=N*k*L, ncol=N*k*L)
  for (nn in c(1:N)){
    tmp1 = kronecker(theta, t(eta[,nn]))
    tmp2 = t(tmp1)%*%invSig
    print(dim(tmp2))
    print(dim(tmp1))
    AinvSig[((nn-1)*L*k+1):(nn*L*k),((nn-1)*p+1):(nn*p)] = tmp2
    AinvSigA[((nn-1)*L*k+1):(nn*L*k,(nn-1)*L*k+1:nn*L*k] = tmp2*tmp1
  }
  
  invbigK = kronecker(invK, diag(nrow=L*k))
  
  Sig = (invbigK + AinvSigA) / diag(nrow=N*k*L)
  m = Sig%*%(AinvSig%*%as.matrix(as.vector(y)))
  zeta_vec = m + t(chol(Sig))%*%as.matrix(as.vector(rnorm(n=N*k*L)))
  
  for(ii in c(1:N)){
    zeta_tmp = zeta_vec[((ii-1)*k*L + 1): (ii*k*L),]
    zeta[ , ,ii] = t(matrix(as.vector(zeta_tmp),nrow=k,ncol=L))
  }
  return(zeta)
}