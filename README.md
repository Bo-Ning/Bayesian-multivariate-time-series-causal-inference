# Bayesian-multivariate-time-series-causal-inference
R code for Ning, B., Ghosal, S. and Thomas, J. (2017)

We provide the code for the new Bayesian causal inference method proposed in our paper: Ning, B., Ghosal, S. and Thomas, J. (2017).
Here we will briefly introduce the model, algorithm from the paper, and provide descriptions for the R files.

### The multivariate basic structural time series model is:

<p align="center">
<img src="http://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20%5Cboldsymbol%20%7BY%7D_t%20%26%20%3D%20%5Cboldsymbol%20%5Cmu_t%20&plus;%20%5Cboldsymbol%20%5Cdelta_t%20&plus;%20%5Cboldsymbol%20X_t%20%5Cboldsymbol%20%5Cbeta%20&plus;%20%5Cboldsymbol%20%5Cepsilon_t%2C%20%5Cquad%20%5Cepsilon_t%20%5Csim%20%5Cmathcal%7BN%7D%28%5Cboldsymbol%200%2C%20%5Cboldsymbol%20%5CSigma%29%5C%5C%20%5Cboldsymbol%20%5Cmu_%7Bt&plus;1%7D%20%26%20%3D%20%5Cboldsymbol%20%5Cmu_t%20&plus;%20%5Cboldsymbol%20%5Ctau_t%20&plus;%20%5Cboldsymbol%20u_t%2C%20%5Cquad%20%5Cboldsymbol%20u_t%20%5Csim%20%5Cmathcal%7BN%7D%28%5Cboldsymbol%200%2C%20%5Cboldsymbol%20%5CSigma_u%29%5C%5C%20%5Cboldsymbol%20%5Ctau_%7Bt&plus;1%7D%20%26%20%3D%20%5Cboldsymbol%20D%20&plus;%20%5Cboldsymbol%20%5CPhi%20%28%5Cboldsymbol%20%5Ctau_t%20-%20%5Cboldsymbol%20D%29%20&plus;%20%5Cboldsymbol%20v_t%2C%20%5Cquad%20%5Cboldsymbol%20v_t%20%5Csim%20%5Cmathcal%7BN%7D%28%5Cboldsymbol%200%2C%20%5Cboldsymbol%20%5CSigma_v%29%5C%5C%20%5Cboldsymbol%20%5Cdelta_%7Bt&plus;1%7D%20%26%20%3D%20-%20%5Csum_%7Bj%3D0%7D%5E%7BS-2%7D%20%5Cboldsymbol%20%5Cdelta_%7Bt-j%7D%20&plus;%20%5Cboldsymbol%20w_t%2C%20%5Cquad%20%5Cboldsymbol%20w_t%20%5Csim%20%5Cmathcal%7BN%7D%28%5Cboldsymbol%200%2C%20%5Cboldsymbol%20%5CSigma_w%29%5C%5C%20%5Cend%7Balign*%7D">
</p>


### Two-stage algorithm
#### Stage 1: EMVS step
Use EMVS algorithm to estimate 
![beta](http://latex.codecogs.com/gif.latex?%5Cboldsymbol%20%5Cbeta) and obtain 
![EMVS](http://latex.codecogs.com/gif.latex?%5Ctilde%7B%5Cboldsymbol%20Y%7D_t%20%3D%20%5Cboldsymbol%20Y_t%20-%20%5Cboldsymbol%20X_t%20%5Chat%7B%5Cboldsymbol%20%5Cbeta%7D)

#### Stage 2: MCMC step
(a) Generate 
![alpha](http://latex.codecogs.com/gif.latex?%5Cboldsymbol%20%5Calpha_%7B1%3AT%7D)
using the Kalman-filter and simulation smoother method

(b) Generate 
![phi](http://latex.codecogs.com/gif.latex?%5Cboldsymbol%20%5CPhi)
using the Metropolis-Hastings algorithm

(c) Generate
![D](http://latex.codecogs.com/gif.latex?%5Cboldsymbol%20D)
from its posterior

(d) Generate
covariance matrices from their respective 
![G](http://latex.codecogs.com/gif.latex?%5Cmathcal%7BG%7D)-Wishart posterior 

(e) Go to Step (a) and repeat until the chain converges.

Skip Step (b) and (c) if no stationarity restriction is imposed on 
![tau](http://latex.codecogs.com/gif.latex?%5Cboldsymbol%20%5Ctau_t)

### Our R code includes:

example.R: uses to simulate a dataset. The simulation process the same as the one described in that paper.

BayesianCausalImpact folder: includes the code for runing the causal inference analysis

 -- two.stage.estimate.R: the main code use to run the two-stage algorithm
 
 -- estimate.counterfactual.R: the code for in Stage 1, the output is to obtain 
    ![Ytilde](http://latex.codecogs.com/gif.latex?%5Ctilde%7B%5Cboldsymbol%20Y%7D_t)
    
 -- EMVS.R: the code for conducting the EMVS algorithm including stationary/nonstationary/misspecified models
 
 -- MCMC.multivariate.ssm.R: the code for running MCMC in Stage 2, the outputs are the parameter posterior draws
 
 -- koopmanfilter.R: the code for sampling 
 ![alpha](http://latex.codecogs.com/gif.latex?%5Cboldsymbol%20%5Calpha_%7B1%3AT%7D)
 
 -- kalmflter.R: the code for Kalman-filter and Backward smoother
 
 -- stationaryRestrict.R & varp.R: the codes for making stationarity constraint on 
 ![tau](http://latex.codecogs.com/gif.latex?%5Cboldsymbol%20%5Ctau_t)
 
 -- MultiCausalImpact.R: the code for conducting our new causal inference method
 
 -- MungeMatrix.R: the code to convert a computationally singular matrix to non-singular
 
 -- ks.distance-plot.R: the code to generate plots of KS distances
 
 
