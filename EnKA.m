function [ T_a ] = EnKA( T,Yp,y,R,alpha,pert_stat)
%% Efficient implementation of the ensemble Kalman analysis step. 
% with perturbed observations.
%
% Dimensions: No = Number of observations to assimilate.
%                      Np = Number of parameters to update.
%                      Ne = Number of ensemble members. 
%  
% -----------------------------------------------------------------------
% Inputs: 
%
% T     => Np x Ne matrix containing an ensemble of Ne prior
%       parameter column vectors each with Np entries.
%        
% Yp   => No x Ne matrix containing an ensemble of Ne predicted
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
%
% -----------------------------------------------------------------------
% Outputs:
% 
% T_a => Np x Ne matrix containing an ensemble of Ne posterior
%      (i.e. the analysis) column vectors each with Np entries.
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
Ne=size(Yp,2); No=size(Yp,1);
if numel(R)==1
    R=R.*eye(No);
elseif numel(R)==No
    R=diag(R);
end
if isempty(alpha)
    alpha=1;
end
alpha_pert=(~pert_stat)+pert_stat*alpha;
%Y=y*ones(1,Ne);+sqrt(alpha_pert).*sqrt(R)*randn(No,Ne); % Good to reperturb on each pass to reduce sampling error
%Y=y*ones(1,Ne);
Y=repmat(y,1,Ne);
perts=sqrt(alpha_pert).*sqrt(R)*randn(No,Ne);

% Define useful shorthands to avoid repeated calculations:
%A=anomaly(T);             % Np x Ne parameter anomaly matrix.
A=T-mean(T,2);
%HEold=anomaly(HX);       	  % No x Ne predicted observation anomaly matrix.
HE=Yp-mean(Yp,2);
%disp(HE-HEold);
HEt=HE';                  % Ne x No transposed predicted observation anomaly matrix. 
Inn=Y-(Yp+perts);                 % No x Ne innovation matrix. Perturb pred obs, not obs in line with van Leeuwen 2020.

% Covariance matrices (scaled by the number of ensemble members)
C_AHE=A*HEt;               % Np x No parameter-predicted observation error covariance matrix. 
C_HEHE=HE*HEt;             % No x No predicted observation error covariance matrix. 
%Apinv=pinv(A);
%ApinvA=Apinv*A;
%C_HEHE=HE*ApinvA*((HE*ApinvA)'); % From Evensen (2019)
aC_DD=(Ne*alpha).*R;       % No x No observation error covariance matrix (scaled by alpha as well).

% Kalman analysis step.
K=C_AHE/(C_HEHE+aC_DD);      % Kalman gain for the parameters.
T_a=T+K*Inn;                 % Analysis.
%toc;


end
