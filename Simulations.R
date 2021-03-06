

#model parameters
shape=6.2
mean=5.6
scale=mean/shape

sd=sqrt(shape*scale^2)

###### Simulations ######
easy_cov = simulate_trials(n=100, mean=mean, sd=sd, distribution="gamma", rl0="random", difficulty = c("easy"))
medium_cov = simulate_trials(n=100, mean=mean, sd=sd, distribution="gamma", rl0="random", difficulty = c("medium"))
hard_cov = simulate_trials(n=100, mean=mean, sd=sd, distribution="gamma", rl0="random", difficulty = c("hard"))


boxplot(easy_cov$easy$deconvolution_average_correlation
        , medium_cov$medium$deconvolution_average_correlation
        , hard_cov$hard$deconvolution_average_correlation
        , names=c("Easy","Medium","Hard")
        , main="100 simulated outbreaks, mean=5.6, sd=2.25, Correlation")

boxplot(easy_cov$easy$deconvolution_average_rmse
        , medium_cov$medium$deconvolution_average_rmse
        , hard_cov$hard$deconvolution_average_rmse
        , names=c("Easy","Medium","Hard")
        , main="100 simulated outbreaks, mean=5.6, sd=2.25, RMSE")

boxplot(  easy_cov$easy$deconvolution_simple_rmse
          , easy_cov$easy$deconvolution_random_rmse
          , easy_cov$easy$deconvolution_ridge_rmse
          , easy_cov$easy$deconvolution_rl_rmse
          , easy_cov$easy$deconvolution_fourier_rmse
          , easy_cov$easy$deconvolution_frequency_rmse
          , easy_cov$easy$deconvolution_average_rmse
          , names=c("Simple","Random","Ridge","RL","Fourier","Frequency","Average")
          , main="100 easy simulated outbreaks, mean=5.6, sd=2.25, RMSE")
boxplot(  medium_cov$medium$deconvolution_simple_rmse
          , medium_cov$medium$deconvolution_random_rmse
          , medium_cov$medium$deconvolution_ridge_rmse
          , medium_cov$medium$deconvolution_rl_rmse
          , medium_cov$medium$deconvolution_fourier_rmse
          , medium_cov$medium$deconvolution_frequency_rmse
          , medium_cov$medium$deconvolution_average_rmse
          , names=c("Simple","Random","Ridge","RL","Fourier","Frequency","Average")
          , main="100 medium simulated outbreaks, mean=5.6, sd=2.25, RMSE")
boxplot(  hard_cov$hard$deconvolution_simple_rmse
          , hard_cov$hard$deconvolution_random_rmse
          , hard_cov$hard$deconvolution_ridge_rmse
          , hard_cov$hard$deconvolution_rl_rmse
          , hard_cov$hard$deconvolution_fourier_rmse
          , hard_cov$hard$deconvolution_frequency_rmse
          , hard_cov$hard$deconvolution_average_rmse
          , names=c("Simple","Random","Ridge","RL","Fourier","Frequency","Average")
          , main="100 hard simulated outbreaks, mean=5.6, sd=2.25, RMSE")



###### Current Outbreak ######
#symptom curve - derived from confirmed cases - 1/25
outbreak_cov=c(rep(4,times=10), 0, rep(1,times=4),17,59,77,93,149,131,259,457,164)

plot(outbreak_cov)

#standard model inputs
parms=c(mean,sd)
distribution="gamma"
#simple
simple=deconvolve_infection_curve_simple(outbreak_cov,5.6)
#random
random=deconvolve_infection_curve_random(outbreak_cov,generate_incubation_period, trials=100, distribution="gamma", parms=parms)
#ridge
incubation_matrix = incubation_period_distribution_matrix(length(outbreak_cov), distribution = distribution, parms = parms)
ridge=deconvolve_infection_curve_ridge(outbreak_cov,incubation_matrix)
#RL
incubation_matrix = incubation_period_distribution_matrix(length(outbreak_cov), distribution = distribution, parms = parms,rl = TRUE)
rl=deconvolve_infection_curve_rl(outbreak_cov,incubation_matrix, random)
#fourier filter - doesn't do well with shorter outbreaks. testing on simulated outbreaks still has strong correlation and low RMSE
fourier_filter_kernel = construct_fourier_filter_kernel(distribution = distribution, parms = parms)
fourier=deconvolve_infection_curve_fourier(outbreak_cov,fourier_filter_kernel)
#frequency filter
frequency_matrix = find_frequency_matrix(length(outbreak_cov), distribution = distribution, parms = parms)
matrix=frequency_matrix
incubation_length = determine_incubation_length(distribution, parms)
frequency=deconvolve_infection_curve_frequency(outbreak_cov,frequency_matrix,incubation_length)

#list of estimates to average
deconvolutions=list(c(simple)
                    , random
                    , ridge
                    , rl
                    #, frequency
                    , c(frequency,0,0,0))

average=deconvolve_infection_curve_average(deconvolutions)

plot(outbreak_cov)
plot(average, main='Wuhan Outbreak Infection Curve')
for(method in deconvolutions){
  plot(method)
}
sum(simple)
sum(outbreak_cov)

symptom_onset_curve=outbreak_cov
median_incubation_period=5.6





# The deconvolution backs us up a full incubation period, so any deconvolution is 8 days behind current reporting. We need to forecast ahead those 8 days to get to "now." I haven't been able to get anything to forecast one day, if that. 




