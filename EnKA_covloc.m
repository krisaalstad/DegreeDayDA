function [ T_a ] = EnKA_covloc( T, Yp, y, R, alpha, pert_stat, loc_Yp, loc_theta,c)
%% Efficient implementation of the ensemble Kalman analysis step with 
%% covariance localisation and perturbed observations.
%
% Dimensions: No = Number of observations to assimilate.
%             Np = Number of parameters to update.
%             Ne = Number of ensemble members. 
%  
% -----------------------------------------------------------------------
% Inputs: 
%
% T     => Np x Ne matrix containing an ensemble of Ne prior
%       parameter column vectors each with Np entries.
%        
% Yp    => No x Ne matrix containing an ensemble of Ne predicted
%         observation column vectors each with No entries.
% 
% y     => No x 1 matrix containing the batch of unperturbed observations.
%
% R     => Observation error covariance matrix. 
%
% alpha => Inflation factor for R in the case of multiple-DA (MDA). Set
%            alpha to 1 or empty for the standard Kalman analysis equation.
%
% pert_stat => True/false (0/1) status on scaling the observation
%                     perturbation in the case of MDA.
% loc_Yp =>  No x Ne matrix containing the (gridpoint) location of the predicted observations 
%        
% loc_theta => Np x Ne matrix containing the (gridpoint) location of the prior parameters     
% 
% c     => tapering distance
% 
% -----------------------------------------------------------------------
% Outputs:
% 
% T_a => Np x Ne matrix containing an ensemble of Ne posterior
%      (i.e. the analysis) column vectors each with Np entrieobs_times_enss.
%
% -----------------------------------------------------------------------
% N.B. The analysis can also be implemented in batch mode as an Ensemble Smoother
% if all the No perturbed observations and predicted observations within the batch
% are collected and stored in HX and Y. The standard sequential
% implementation is also possible, in such a case No would be the number of
% observations at the current time step. We do not follow the routine in
% Evensen (2003) for the perturbed observation EnKF because this requires
% storing a Ne x Ne identity matrix which is prohibitive for large (e.g. Ne
% = 10^6) ensembles. Instead scaled covariance matrices are computed directly
% in ensemble space.
% Code by Kristoffer Aalstad (December 2015, last revision June 2019).



%% Scheme:

%% Change obs perturbation to alpha scaled or not 
Ne=size(Yp,2); No=size(Yp,1); Np=size(T,1); 
% define observation covariance matrix from observation uncertainties in R
if numel(R)==1
    R=R.*eye(No);
elseif numel(R)==No
    R=diag(R);
end
% make sure alpha is defined
if isempty(alpha)
    alpha=1;
end
% calculate perturbations added to the predicted observations
% dependend on whether or not MDA is used
alpha_pert=(~pert_stat)+pert_stat*alpha;
perts=sqrt(alpha_pert).*sqrt(R)*randn(No,Ne);
% expand the observations into a No x Ne matrix
Y=repmat(y,1,Ne);

% Define useful shorthands to avoid repeated calculations:
A=T-mean(T,2);            % Np x Ne parameter anomaly matrix.
HE=Yp-mean(Yp,2);         % No x Ne predicted observation anomaly matrix.
HEt=HE';                  % Ne x No transposed predicted observation anomaly matrix. 
Inn=Y-(Yp+perts);         % No x Ne innovation matrix. Perturb pred obs, not obs in line with van Leeuwen 2020.

% Covariance matrices (scaled by the number of ensemble members)
C_AHE=A*HEt;               % Np x No parameter-predicted observation error covariance matrix. 
C_HEHE=HE*HEt;             % No x No predicted observation error covariance matrix. 
aC_DD=(Ne*alpha).*R;       % No x No observation error covariance matrix (scaled by alpha as well).

% add localization on the covariance matrices C_AHE and C_HEHE following
% Sakov and Bertino (2009) eq 7
    % calculate distance matrices
dist_C_AHE = sqrt( ( repmat(loc_theta,1,No) - repmat(loc_Yp',Np,1) ).^2);
dist_C_HEHE = sqrt( ( repmat(loc_Yp,1,No) - repmat(loc_Yp',No,1) ).^2);
    % calculate distance-based tapering matrix using Gaspari Cohn function
    % with one localization constant c
rho_C_AHE = GC(dist_C_AHE,c);
rho_C_HEHE = GC(dist_C_HEHE,c);
    % calculate localized covariance matrices
loc_C_AHE  =  C_AHE .* rho_C_AHE;
loc_C_HEHE = C_HEHE .* rho_C_HEHE;
    
% Kalman analysis step.
K = loc_C_AHE/(loc_C_HEHE + aC_DD);      % Kalman gain for the parameters.
T_a = T + K*Inn;                 % Analysis.

end
