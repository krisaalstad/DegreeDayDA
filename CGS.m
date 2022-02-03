function x=CGS(rho,mu,sigma,N)
%% x=CGS(d,c): Cholesky-based Gaussian Sampler
% Generates samples from a joint n-dimensional Gaussian allowing for 
% non-zero correlations between variables by using a Cholesky decomposition
% of the covariance matrix and pseudorandom sampling from the standard
% Gaussian.
% 
% This function can be used to generate spatially correlated Gaussian prior
% distributions for data assimialtion. It can also be used to generate
% Gaussian priors with correlations between different variables.
%
% Inputs: 
% rho = Correlation matrix (n x n) with entries ranging from -1 to 1
% mu = Mean vector (n x 1)
% sig = Standard deviation vector (n x 1)
% N = Number of (pseudo)random samples to generate.
% 
% Output:
% x = Matrix (n x N) of correlated random samples where n=number of 
% dimensions, N=number of samples.
%
% Note that the covariance matrix C is symmetric and given by 
%   C=[sigma_1^2 rho_12*sigma_1*sigma_2 ... rho_1n*sigma_1*sigma_n; 
%      rho_21*sigma_2*sigma_1 sigma_2^2 ... rho_2n*sigma_2*sigma_n;
%      ....
%     rho_n1*sigma_n*sigma_1 rho_n2*sigma_n*sigma_2 ... sigma_n^2]
% Thereby,
%   C=rho.*(sigma*(sigma')) where .* is the elementwise product, the outer
%   product of the standard deviation vector sigma*(sigma') yields a 
%   symmetric matrix, and rho is the correlation matrix.
% So C can be built up from the correlation matrix rho_ij and the standard
% deviation vector sigma_i. This is helpful because it's often easier or at
% least more familiar to deal with standard deviations and correlations
% rather than with covariances directly.
%
%{
% Simple 2-D example
rho=[1 0.8; 0.8 1]; % Diagnoals ("self"correlation) always have to be 1
mu=[0;1];
sig=[0;2];
N=1e4;
x=CGS(rho,mu,sig,N);
scatter(x(1,:),x(2,:),150,'.');
%}
% For a more complex example that also uses the Gaspari-Cohn (GC) function
% for a large 2-D grid, see dancing_blobs.m.
%% Main code:

n=size(sigma,1); % Number of dimensions.
% This can be the number of dimensions in space and/or time or more
% generally just the number of different variables.

% Control that sigma and mu are column vectors of the same size.
if (size(sigma,2)>1)||(size(mu,2)>1)||numel(mu)~=n||numel(sigma)~=n
    error('CGS only accepts sigma and mu as (n x 1) column vectors');
end

C=rho.*(sigma*(sigma')); % Build up the covariance matrix.
S=chol(C,'lower'); % Cholesky decomposition of the covariance: C=S*S'
% S is a lower triangular matrix loosely cooresponding to a standard
% deviation matrix.

z=randn(n,N); % Generate (n x N) random draws from a standard Gaussian N(0,1).
x=mu+S*z; % Generate random draws from the n-dimensional Gaussian N(mu,C).


end









