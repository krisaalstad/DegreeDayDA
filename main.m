%% DegreeDayDA: Degree-Day snow model with Data Assimilation
% 
% This script can be used to demonstrate or explore data assimilation (DA).
% It sets up an Ensemble Smoother (ES) or Ensemble Smoother with Multi Data
% Assimilation (ES-MDA) in a synthetic experiment, i.e. the observations
% used in the assimulation are generated by the model. 
% 
% We use a very simple model, a degree day model, to simulate the evolution
% of the snowpack through one hydrological year. A degree day model takes
% daily temperature (in degrees Celsius) and precipitation (in mm/day) as
% inputs to simulate the evolution of the snow water equivalent (SWE). Snow
% accumulation is diagnosed using a fixed (assumed known) threshold
% temperature for snowfall and by applying a snowfall multiplier (uncertain
% parameter) to the precipitation to account for potential biases and
% unresolved proceses like snowdrift. Snow ablation is diagnosed for
% positive degree days (days with a temperature > 0 degrees Celsius) by
% scaling positive temperatures with a degree-day factor (uncertain
% parameter) which parametrizes all the unresolved processes related to
% snowmelt and sublimation. We model the evolution of the snow depth (in mm
% w.e.) in a number of grid
% 
% The idea in a synthetic experiment is to:
% 
%   1) Perform a so-called "truth run" by setting some true value for the
%   parameters in the respective grid cells. For the rest of the experiment
%   (other than the validation), we pretend that we don't know the true
%   values of the two parameters (degree day factor, snowfall multiplier)
%   and states (SWE) other than in the final validation. 
%   2) Generate synthetic observations by perturbing the true SWE for a
%   subset of days in the truth runs to mimic reality with sparse and
%   imperfect observations.
%   3) Set a prior distribution on the uncertain parameters (remember, we
%   wouldn't know the truth in practice) and run an ensemble of degree-day
%   models for the entire water year with different values for the two
%   parameters.
%   4) Assimilate the synthetic observations at the end of the water year
%   to update the parameters. Then rerun the water year with the posterior
%   (updated) parameters to get the posterior SWE.
%   5) Perform some validation, i.e. compare the prior and posterior to the
%   truth, to see if the assimilation is working as intended. Potentially
%   compare the effect of using different: (i) Assimilation schemes, (ii)
%   Observation errors, (iii) Observation density, (iv) Observation types,
%   (v) Snow models, (vi) Localization routines.    



% points close to Finse, Norway. We use real climatic forcing from
% down-scaled weather data TopoSCALE


% THINGS TO DO TO IMPROVE CODE
% make it possible to define global parameter, a model parameter that is
% the same for the entire model domaine, and local parameter, a model
% parameter that can vary from grid point to grid point

% make it possible to define observations per grid point / location

% rename variables:
% x a => Dday_factor / Dday_factor_true
% x b => P_factor / P_factor_true
% x D => snowdepth
% x Dobs => snowdepth_obs
% x Dtrue => snowdepth_true
% here => obs_times
% x ysd => y_std
% x tobs => t_obs


%% start code

clearvars;

% set parameters

% Data assimilation parameters
Ne          = 1e2; % ensemble size 
Na          = 4;   % number of MDA iterations, Na=1 corresponds to ES

% model run 
t_start     = '01-Sep-2018';
t_end       = '01-Sep-2019';
N_gridp     = 2;

% parameter definition and truth run settings:
% a 'local' parameter can have a different value for each grid cell;
% a 'global' parameter has one singular value for all grid cells.
% (option to set P_factor to 'global' not implented yet)
Dday_factor_type    = 'global';

% set truth values 
% literature value for the ddf for snow 2.5 to 11.6 mm/d/K
% P_factor value must be positive
ddf_lit_values      = [6; 10]; % mm/d/K, if 
P_lit_values        = [0.5; 2];

if strcmp(Dday_factor_type, 'global') == 1
    Dday_factor_true = ddf_lit_values(1);
elseif strcmp(Dday_factor_type, 'local') == 1
    Dday_factor_true = ddf_lit_values(:);    
else
    disp('error: Dday_factor_type not set to ''local'' or ''global''');
end
P_factor_true = P_lit_values(:);


%  synthetic observations
y_std         = 20; % standard deviation of the error term added to the synthetic SWE observations
t_obs        = [datenum('01-Jan-2019'); datenum('01-Feb-2019');...
                datenum('01-Mar-2019'); datenum('01-Apr-2019');...
                datenum('01-May-2019'); datenum('01-Jun-2019');...
                datenum('01-Jul-2019'); datenum('01-Aug-2019')];



%% load climatic forcing

% load real-data forcing
load('TopoSCALE.mat');   % contains weather data for two grid points 
f.P=double(f.P).*f.P_sf; % unpack precipitation data
f.T=double(f.T).*f.T_sf; % unpack temperature data


% select model period and define model time t
t=datenum(t_start):1:datenum(t_end); t=t';

% translate hourly values into daily values
Nt=numel(t);
P=zeros(Nt,2);
T=P;
for j=1:Nt
   here=f.t>t(j)&f.t<=(t(j)+1);
   T(j,:)=mean(f.T(:,here),2);
   P(j,:)=sum(f.P(:,here),2); 
end
clear here

%% 1) Synthetic model run creating the 'TRUTH' 
%D=zeros(size(P));
% run the degree day model for defined period, climate input, and model
% parameters
if strcmp(Dday_factor_type, 'global') == 1
    Dday_factor_DDMinput = Dday_factor_true * ones(N_gridp, 1);
else
    Dday_factor_DDMinput = Dday_factor_true;
end
snowdepth_true = DDM(t,P,T,Dday_factor_DDMinput,P_factor_true);
% set result to truth 
% snowdepth_true=snowdepth;

%% 2) Generate synthetic observations 
% add observation error to truth
snowdepth_obs=snowdepth_true+y_std.*randn(size(snowdepth_true));
snowdepth_obs=max(snowdepth_obs,0); % Observed SWE can not be negative


% and select observations by defining the time (=day) of observation
No=numel(t_obs);
tmp=zeros(No,2);
obs_times=zeros(Nt,1,'logical');
for j=1:No
    obs_timesj = t==t_obs(j);
    obs_times = obs_times + obs_timesj;
    tmp(j,:) = snowdepth_obs(obs_timesj,:);
end
obs_times=logical(obs_times);
snowdepth_obs=tmp;
clear tmp

%{
figure(1); clf;
plot(t,snowdepth);
datetick('x','keepticks','keeplimits');
hold on;
scatter(t_obs,snowdepth_obs,150,'.');
%}


%% 3) Set a prior distribution on the uncertain parameters
% Dday_factor=exp(log(5)+1.*randn(2,Ne));
P_factor=exp(log(1)+1.*randn(2,Ne));

if strcmp(Dday_factor_type, 'global') == 1
    Dday_factor=exp(log(5)+1.*randn(1,Ne));
    Dday_factor_DDMinput = repmat(Dday_factor, [N_gridp, 1]);
elseif strcmp(Dday_factor_type, 'local') == 1
    Dday_factor=exp(log(5)+1.*randn(N_gridp,Ne));
    Dday_factor_DDMinput = Dday_factor;
else
    disp('ERROR: Dday_factor not defined because ''Dday_factor_type'' is either ''local'' nor ''global''')
end



%% 4) Run the model and assimilate the synthetic observations 

% run the data assimilation Na times (user-defined number of MDA
% iterations), and run the model one more time to calculate the model
% result with the final parameter ensembles (Na+1)
for ell=1:(Na+1)    
    % model run
    snowdepth=DDM(t,P,T,Dday_factor_DDMinput,P_factor);
    
    % Data Assimilation
    theta=[Dday_factor;P_factor];    % theta combines all parameters that can be updated
    
    % save ensemble prior and post (ensembles before and after the
    % assimilation(s), the intermediate results in MDA are not saved)
    if ell==1
        thetapri = theta;
        if strcmp(Dday_factor_type, 'global') == 1
            Dday_factor_pri = Dday_factor(1,:);
        elseif strcmp(Dday_factor_type, 'local') == 1
            Dday_factor_pri=Dday_factor(1:N_gridp,:);
        end
        P_factor_pri=P_factor(1:N_gridp,:);
        snowdepth_pri=snowdepth;
    elseif ell==(Na+1)
        thetapost=theta;
        if strcmp(Dday_factor_type, 'global') == 1
            Dday_factor_post = Dday_factor(1,:);
        elseif strcmp(Dday_factor_type, 'local') == 1
            Dday_factor_post = Dday_factor(1:N_gridp,:);
        end
        P_factor_post=P_factor(1:N_gridp,:);
        snowdepth_post=snowdepth;
    end
    
    % peform data assimilation; except for the Na+1 iteration
    if ell<=Na
        % perform anamorphosis on the parameters to assure Gaussian
        % distribution
        phi=log(theta);
        % define the model predicted observations matrix from the model results
        Yp=snowdepth(obs_times,:,:); % time, space, ensemble
        Yp=reshape(Yp,No*2,Ne);
        % define the observation matrix
        y=snowdepth_obs;
        y=y(:);
        % define observation covariance (measure for observation
        % uncertainty
        R=y_std.^2;
        alpha=Na; % number of MDA iterations
        pert_stat=1; % 1 for MDA, 0 for ES
        % update parameters with the EnKA function
        phi=EnKA(phi,Yp,y,R,alpha,pert_stat);
        % transform parameters back
        theta=exp(phi);
        % update model parameters
        if strcmp(Dday_factor_type, 'global') == 1
            Dday_factor =theta(1,:);
            Dday_factor_DDMinput = repmat(Dday_factor, [N_gridp, 1]);     
            P_factor    =theta(2:N_gridp+1,:);
        elseif strcmp(Dday_factor_type, 'local') == 1
            Dday_factor =theta(        1:N_gridp,  :);
            Dday_factor_DDMinput = Dday_factor;
            P_factor    =theta(N_gridp+1:2*N_gridp,:);
        end
    end
end

%return

%% 5) Perform some validation

figure(1); clf;
for loc=1:2
    subplot(2,2,loc);
    plot(t,squeeze(snowdepth_pri(:,loc,:)),'LineWidth',0.5,'Color',[0.8 0 0 0.1]); hold on;
    plot(t,squeeze(snowdepth_post(:,loc,:)),'LineWidth',0.5,'Color',[0 0 0.8 0.1]); hold on;
    pt(1)=plot(t,squeeze(mean(snowdepth_pri(:,loc,:),3)),'LineWidth',2,'Color',[0.8 0 0 1]); hold on;
    pt(2)=plot(t,squeeze(mean(snowdepth_post(:,loc,:),3)),'LineWidth',2,'Color',[0 0 0.8 1]); hold on;
    snowdepthtl=snowdepth_true(:,loc);
    pt(3)=plot(t,snowdepthtl,'-k','LineWidth',2);
    pt(4)=scatter(t_obs,snowdepth_obs(:,loc),100,'o','MarkerFaceColor',[0.7 0.7 0],'MarkerEdgeColor',[0 0 0]);
    %scatter(t_obs,snowdepth_obs(:,loc),250,'.','MarkerEdgeColor',[0 0 0]);
    leg=legend(pt(:),{'Prior','Post','Truth','Obs'},'Interpreter','Latex',...
        'FontSize',16,'Location','NorthWest');
    clear pt;
    yl=5*max(snowdepthtl); 
    ylim([0 2e3]);
    xlim([min(t) max(t)]);
    datetick('x','keepticks','keeplimits');
    axis xy;
    ylabel('SWE [mm]','Interpreter','Latex','FontSize',20);
    box on; %grid on;
    set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.005, 0.005]);
    set(groot, 'defaultAxesTickLabelInterpreter','Latex');
    set(gca,'TickLabelInterpreter','latex'); 
    title(sprintf('States in location %d',loc),'Interpreter','Latex','FontSize',16);
    
    subplot(2,2,2+loc);
    if strcmp(Dday_factor_type, 'global') == 1
        xpri=Dday_factor_pri(1,:); ypri=P_factor_pri(loc,:);
        xpost=Dday_factor_post(1,:); ypost=P_factor_post(loc,:);
        xtrue = Dday_factor_true * size(N_gridp); ytrue = P_factor_true;
    else
        xpri=Dday_factor_pri(loc,:); ypri=P_factor_pri(loc,:);
        xpost=Dday_factor_post(loc,:); ypost=P_factor_post(loc,:);
        xtrue = Dday_factor_true; ytrue = P_factor_true;
    end
    alphis=0.4;
    pt(1)=scatter(xpri,ypri,4e1,'o','MarkerFaceColor',[0.8 0 0],'MarkerFaceAlpha',alphis,...
        'MarkerEdgeColor',[0 0 0]); hold on;
    pt(2)=scatter(xpost,ypost,4e1,'s','MarkerFaceColor',[0 0 0.8],'MarkerFaceAlpha',alphis,...
        'MarkerEdgeColor',[0 0 0]); 
    pt(3)=scatter(xtrue(loc),ytrue(loc),250,'p','MarkerEdgeColor',[0 0 0],'LineWidth',2,...
        'MarkerFaceColor',[1 1 1]);
    leg=legend(pt(:),{'Prior','Post','Truth'},'Interpreter','Latex',...
        'FontSize',16,'Location','NorthEast');
    clear pt;
    xlim([0 30]);
    ylim([0 5]);
    axis xy;
    ylabel('Snowfall multiplier, $b$ [-]','Interpreter','Latex','FontSize',20);
    xlabel('Degree-day factor, $a$ [mm/K/d]','Interpreter','Latex','FontSize',20);
    box on; %grid on;
    set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.005, 0.005]);
    set(groot, 'defaultAxesTickLabelInterpreter','Latex');
    set(gca,'TickLabelInterpreter','latex'); 
    axis square;
    title(sprintf('Parameters in location %d',loc),'Interpreter','Latex','FontSize',16);
end

% print('-djpeg','results','-r300','-opengl');