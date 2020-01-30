
shape=6.2
mean=5.6
scale=mean/shape
sd=sqrt(shape*scale^2)
parms=c(mean,sd)
distribution="gamma"
curve=outbreak_cov
incubation_matrix = incubation_period_distribution_matrix(length(curve), distribution = distribution, parms = parms)
incubation_matrix_rl = incubation_period_distribution_matrix(length(curve), distribution = distribution, parms = parms,rl = TRUE)
fourier_filter_kernel = construct_fourier_filter_kernel(distribution = distribution, parms = parms)
frequency_matrix = find_frequency_matrix(length(curve), distribution = distribution, parms = parms)
matrix=frequency_matrix
incubation_length = determine_incubation_length(distribution, parms)


deconvolve_single_curve = function(curve, parms)
{
  #simple
  simple=deconvolve_infection_curve_simple(curve,parms[1])
  #random
  random=deconvolve_infection_curve_random(curve,generate_incubation_period, trials=35, distribution="gamma", parms=parms)
  #ridge
  ridge=deconvolve_infection_curve_ridge(curve,incubation_matrix)
  #RL
  rl=deconvolve_infection_curve_rl(curve,incubation_matrix_rl, random)
  #fourier filter - doesn't do well with shorter outbreaks. testing on simulated outbreaks still has strong correlation and low RMSE
  fourier=deconvolve_infection_curve_fourier(curve,fourier_filter_kernel)
  #frequency filter
  frequency=deconvolve_infection_curve_frequency(curve,frequency_matrix,incubation_length)
  
  #list of estimates to average
  deconvolutions=list(c(simple)
                      , random
                      , ridge
                      , rl
                      #, frequency
                      , c(frequency,0,0,0))
  
  average=deconvolve_infection_curve_average(deconvolutions)
  
  curves = list(curve=curve, simple=simple, random=random, ridge=ridge, rl=rl, fourier=fourier, frequency=frequency, average=average)
  return(curves)
}

deconvolve_single_curve(outbreak_cov, parms = parms)

plot(y)
y=outbreak_cov[1:29]
numlags = 8 # choose the model..
n = length(y)
ar1 = ar.ols(y, demean = F,intercept=T, order.max = numlags,aic = F) # standard OLS
plot(ar1$res)
scaled.res = scale(ar1$res)
#summ(scaled.res)
b1 = y[1:numlags]
R = 1000 # how many bootstrap samples do you want?
res.star = matrix(nrow = (n-numlags),ncol = R) 
obs.star = matrix(nrow = n,ncol = R) # will hold the bootstrapped series

for (i in 1:R){
  res.star[,i] = sample(na.omit(scaled.res),(n-numlags), replace = T)
  obs.star[1:numlags,i] = b1 # for the first obs we plug the original data
  for (j in (numlags+1):n){
    obs.star[j,i] = ar1$x.intercept + ar1$ar%*%obs.star[(j-1):(j-numlags),i] + res.star[(j-numlags),i]
  }}
# Sanity check:
for(i in 1:30){
  plot(obs.star[,i], ty = "l") 
}
obs.star=-obs.star

final_results=list()
for (i in 477:R) {
  print(i)
  obs.star[obs.star[,i]<0,i]=0
  results = deconvolve_single_curve(curve = obs.star[,i], parms = parms)
  final_results[[i]]=results$average
}
test=c(0,0,0,0,-1,2)
test[test<0]=1

for(i in 1:8){
  print(i)
  plot(final_results[[i]], main=i) 
}
final_results[[2]]
results$curve





