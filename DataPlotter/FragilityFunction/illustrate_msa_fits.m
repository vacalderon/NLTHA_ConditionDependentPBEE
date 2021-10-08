% Example calculations to demonstrate the use of fragility fitting
% functions. These calculations are based on the following paper:
%
% Baker, J. W. (2015). “Efficient analytical fragility function fitting 
% using dynamic structural analysis.” Earthquake Spectra, 31(1), 579-599.
%
% Created by Jack Baker
% 2/4/2013
% Modified by Jack Baker, 1/25/2017, to update citation information


% example data: IM levels, number of analyses, and number of collapses
IM =            [0.2    0.3     0.4    0.6      0.7    0.8    0.9    1 ];
num_collapse =  [0      0       0      4        6      13     12     16];      
num_gms =       [40     40      40     40       40     40     40     40];


% estimate fragility function using MLE method (equation 11)
[theta_hat_mle, beta_hat_mle] = fn_mle_pc(IM, num_gms, num_collapse);

% estimate fragility function using MLE method (this version uses probit 
% regression, which produces equivalent answers to equation 11)
[theta_hat_probit, beta_hat_probit] = fn_mle_pc_probit(IM, num_gms, num_collapse);

% estimate fragility function using sum of squared errors between data and
% fragility (equation 12). This approach is not recommended, but is
% provided here for comparison purposes.
[theta_hat_sse, beta_hat_sse] = fn_sse_pc(IM, num_gms, num_collapse);


% compute fragility functions using estimated parameters
x_vals = 0.01:0.01:3; % IM levels to plot fragility function at

p_collapse_mle = normcdf((log(x_vals/theta_hat_mle))/beta_hat_mle); % compute fragility function using equation 1 and estimated parameters
p_collapse_probit = normcdf((log(x_vals/theta_hat_probit))/beta_hat_probit); % compute fragility function using equation 1 and estimated parameters
p_collapse_sse = normcdf((log(x_vals/theta_hat_sse))/beta_hat_sse); % compute fragility function using equation 1 and estimated parameters


%% plot resulting fragility functions

figure
plot(IM,num_collapse./num_gms, '^b', 'linewidth', 2)
hold on
plot(x_vals,p_collapse_mle, '-b', 'linewidth', 2)
plot(x_vals,p_collapse_probit, '--r', 'linewidth', 2)
plot(x_vals,p_collapse_sse, '-.k', 'linewidth', 1)
legh = legend('Observed fractions of collapse', 'Max likelihood fit', 'Probit regression fit', 'Sum of squared errors fit (do not use)', 'location', 'southeast');
set(legh, 'fontsize', 12)
hx = xlabel('IM', 'Fontsize', 14);
hy = ylabel('Probability of collapse', 'Fontsize', 14);
axis([0 3 0 1])







