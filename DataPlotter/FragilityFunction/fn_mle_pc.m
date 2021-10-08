function [theta, beta] = fn_mle_pc(IM, num_gms, num_collapse)
% by Jack Baker, 10/9/2012

% Modified by Gemma Cremen, 1/25/2017, to avoid estimating negative median
% values for the fragility function

% Modified by Jack Baker, 1/25/2017, to update citation information

% Modified by Jaewon Saw and Adam Zsarnoczay, 2/10/2020, to improve initial 
% guess by method of moments; replace "theta" with "params" in mlefit function
% to avoid confusion with median "theta"; avoid negative beta values by
% applying penalty; remove check for negative theta in mlefit by using theta 
% in log space; and avoid zero likelihood values by replacing zeros with realmin


% This function fits a lognormal CDF to observed probability of collapse 
% data using optimization on the likelihood function for the data. 
% These calculations are based on equation 11 of the following paper:

% Baker, J. W. (2015). “Efficient analytical fragility function fitting 
% using dynamic structural analysis.” Earthquake Spectra, 31(1), 579-599.

% INPUTS:
% IM            1xn           IM levels of interest
% num_gms       1x1 or 1xn    number of ground motions used at each IM level
% num_collapse 	1xn           number of collapses observed at each IM level
% 
% OUTPUTS:
% theta         1x1           median of fragility function
% beta          1x1           lognormal standard deviation of fragility function

%% Initial guess for the fragility function parameters theta and beta
% ** Use method of moments **
x0 = [mean(log(IM)), std(log(IM))];

%% Run optimization
options = optimset('MaxFunEvals',1000, 'GradObj', 'off'); %maximum 1000 iterations, gradient of the function not provided
x = fminsearch(@mlefit, x0, options, num_gms, num_collapse, IM) ;
theta = exp(x(1)); % return theta in linear space
beta = x(2);

%% Objective function to be optimized
function [loglik] = mlefit(params, num_gms, num_collapse, IM)

% ** Penalize any negative beta with a very large loglik value **
if params(2) < 0
    loglik = 1e10;
    
else
    % estimated probabilities of collapse, given the current fragility function
    % parameter estimates
    p = normcdf(log(IM), (params(1)), params(2)); 

    % likelihood of observing num_collapse(i) collapses, given num_gms
    % observations, using the current fragility function parameter estimates
    likelihood = binopdf(num_collapse', num_gms', p'); %
    
    % ** Cannot have zero likelihood value, so replace every zero likelihood 
    % value with the smallest positive normalized fixed-point value **
    likelihood(likelihood == 0) = realmin;
    
    % sum negative log likelihood (we take the negative value because we want
    % the maximum log likelihood, and the function is searching for a minimum)
    loglik = -sum(log(likelihood));
end