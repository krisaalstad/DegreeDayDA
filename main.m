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
% parameter that can vary from grid point to grid point -- DONE

% make it possible to define observations per grid point / location  --
% DONE

% rename variables:
% x a => Dday_factor / Dday_factor_true
% x b => P_factor / P_factor_true
% x D => snowdepth
% x Dobs => snowdepth_obs
% x Dtrue => snowdepth_true
% x here => obs_times
% x ysd => y_std
% x tobs => t_obs
% --DONE

% make more than 2 grid points possible
% --DONE

% include spatial correlation in the priors
% and inlcude covariance localization
% -- DONE


%% set model parameters

clearvars;

% model run 
t_start     = '01-Sep-2018';
t_end       = '01-Sep-2019';
N_gridp     = 4; % min 2, max 6 grid points, add literature values to ddf and P_factor in case more than 6 points are wanted
save_figure = 0; % set save results figure: 0 or 1
fig_path    = pwd; % or define path for saved figs, e.g. '~/Documents/presentaties/workmeeting/2022-05-03_snowdepth-kickoff';

% Data assimilation (ES MDA) parameters
Ne          = 100; % ensemble size 
Na          = 4;   % number of MDA iterations, Na=1 corresponds to ES
% rng(1234)        % prescribe SEED of the random number generator. needed if you want to
                   % compare experimental setups with the same pick from Gaussian distribution.

% set localization options
corr_prior_distr    = 1;    % spatial correlation in the prior parameter distributions on/off (=1/0)
cov_localization    = 1;    % spatial covariance localization in the ES updates on/off (=1/0)
% define tapering distance in grid point units    
c_priorcorr         = 2;    % tapering distance in prior correlation and localization    
c_covloc            = 2;    % tapering distance in covariance localization 

% parameter definition and truth run settings:
% a 'local' parameter can have a different value for each grid cell;
% a 'global' parameter has one singular value for all grid cells.
% (option to set P_factor to 'global' not implented yet)
Dday_factor_type    = 'local'; % 'local' or 'global'
 
% set truth values 
% literature value for the ddf for snow 2.5 to 11.6 mm/d/K
% P_factor value must be positive
ddf_lit_values      = [6;    9;   10;  9.5;   7; 6.4]; % mm/d/K, if snowdepth_obs
P_lit_values        = [1; 1.25; 1.55; 1.75; 1.5; 1.3];
    % degree day factor
if strcmp(Dday_factor_type, 'global') == 1
    Dday_factor_true = ddf_lit_values(1);
elseif strcmp(Dday_factor_type, 'local') == 1
    Dday_factor_true = ddf_lit_values(1:N_gridp);    
else
    disp('error: Dday_factor_type not set to ''local'' or ''global''');
end
    % precipitation factor
P_factor_true = P_lit_values(1:N_gridp);


%  synthetic observations
y_std         = 20; % standard deviation of the error term added to the synthetic SWE observations, (default = 20)

% define the time of the observations for each grid point.
% fill up entries in t_obs with small numbers <= Nmaxobs
Nmaxobs     = 11;
t_obs       = zeros(Nmaxobs, N_gridp); 
t_obs(:,1)  = [ datenum('01-Jan-2019'),... % obs time first grid point
                datenum('01-Feb-2019'),...
                datenum('01-Mar-2019'),...
                datenum('15-Mar-2019'),...
                datenum('01-Apr-2019'),...
                datenum('22-Apr-2019'),...
                datenum('01-May-2019'),...
                datenum('15-May-2019'),...
                datenum('01-Jun-2019'),...
                datenum('01-Jul-2019'),...
                datenum('01-Aug-2019'),...
                ]; 
t_obs(:,2)  = [ datenum('15-Mar-2019'),... % obs time second grid point
                datenum('01-Apr-2019'),...                
                datenum('15-May-2019'),...
                4:Nmaxobs...
                ];                
t_obs(:,3)  = [ datenum('01-Apr-2019'),... % obs time third grid point
                2:Nmaxobs...
                ]; 
% t_obs(:,4)  = [ datenum('10-May-2019'),... % obs time fourth grid point
%                 2:Nmaxobs...
%                 ];           
            
% example dates, can be copied in above:
%                 datenum('01-Jan-2019'),... % obs time second grid point
%                 datenum('01-Feb-2019'),...
%                 datenum('01-Mar-2019'),...
%                 datenum('15-Mar-2019'),...
%                 datenum('01-Apr-2019'),...
%                 datenum('15-Apr-2019'),...
%                 datenum('01-May-2019'),...
%                 datenum('15-May-2019'),...
%                 datenum('01-Jun-2019'),...
%                 datenum('15-Jun-2019'),...
%                 datenum('01-Jul-2019'),...
%                 datenum('01-Aug-2019'),...



            
%% load climatic forcing

% load real-data forcing
%load('TopoSCALE.mat');   % contains weather data for two grid points 
% Fetch the input data if you haven't yet.
if ~exist('input','dir')
    disp('Fetching input data');
    urlis='https://www.dropbox.com/s/hi7ky0y0jeh2ob7/input.zip?dl=1';
    %'https://www.dropbox.com/s/9qckrp6tp7v6a81/input.zip?dl=1'; % Changed dl=0 to dl=1
    tarf='input.zip';
    websave(tarf,urlis);
    unzip(tarf,''); % Unzip
    paths; % Update path
end
load('input/TopoSCALE.mat');


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
clear here % clear temporary variable

% expand, or decrease, the climatic forcing to the number of defined grid
% points. At present the forcing is defined for 2 grid points
if N_gridp~=2
    if N_gridp == 1
        T = T(:,1);
        P = P(:,1);
    elseif N_gridp > 2
        T(:,3:N_gridp) = repmat(T(:,1),1,N_gridp-2);
        P(:,3:N_gridp) = repmat(P(:,1),1,N_gridp-2);
    end
end
        

%% 1) Synthetic model run creating the 'TRUTH' 

% define the input model parameters to have [N_gridp,1] even if the degree
% day parameter is the same for all grid points.
if strcmp(Dday_factor_type, 'global') == 1
    Dday_factor_DDMinput = Dday_factor_true * ones(N_gridp, 1);
else
    Dday_factor_DDMinput = Dday_factor_true;
end

% run the degree day model for defined period, climate input, and true model
% parameters and set the result as the truth
snowdepth_true = DDM(t,P,T,Dday_factor_DDMinput,P_factor_true);


%% 2) Generate synthetic observations 

% add observation error to truth
snowdepth_addederror=snowdepth_true+y_std.*randn(size(snowdepth_true));
snowdepth_addederror=max(snowdepth_addederror,0); % Observed SWE can not be negative

% and select observations by finding the model time (=day) of the
% observations
No              = zeros(N_gridp,1);                   % number of observations per grid point
obs_times       = zeros(Nt, N_gridp, 'logical');    % [Nt,Ngridp] array to select the model snow depth value from the true run
obs_times_ens   = zeros(Nt, N_gridp, Ne, 'logical');% [Nt,Ngridp, Ne] array to select the model snow depth values from the ensemble run
snowdepth_obs   = nan(Nmaxobs, N_gridp);            % [Nmaxobs, Ngridp] array containing the observations for every grid point as derived from the true run
for ii = 1:N_gridp
    No(ii) = sum(t_obs(:,ii)> Nmaxobs);
    obs_times(:,ii) = ismember(t, t_obs(1:No(ii), ii) ); % select the model time steps with an observation for grid point ii
    obs_times_ens(obs_times(:,ii),ii,:) = 1; % set observation time for all ensemble members of this grid point ii 
    snowdepth_obs(1:No(ii),ii) = snowdepth_addederror(obs_times(:,ii), ii); % select the snow depth observation at the observation time from the true run with added observation error
end
No_tot = sum(No,1);       % total number of available observations summed over all grid points

% keep track of observation location, needed for the localization
loc_obs = nan(size(snowdepth_obs));
for ii = 1:N_gridp
    loc_obs(1:No(ii),ii) = ii*ones(No(ii), 1);
end


% % control figure to see synthetic observations:
% figure(1); clf;
% plot(t,snowdepth_addederror);
% datetick('x','keepticks','keeplimits');
% hold on;
% plot(t,snowdepth_true, '--');
% scatter(t_obs(t_obs>Nmaxobs),snowdepth_obs(~isnan(snowdepth_obs)),550,'.');



%% 3) Set a prior distribution on the uncertain parameters

% set parameters spatially correlated prior distribution function CGS,
% based on distance between grid points and cut of length c_priorcorr
if corr_prior_distr == 1
    x   = (1:N_gridp)'; % spatial coordinate (dummy units). Column vectors
    d   = sqrt((repmat(x,1,N_gridp)-repmat(x',N_gridp,1)).^2); % distance matrix fot the used grid points   
    rho = GC(d,c_priorcorr); % Gaspari Cohn function
    sig = ones(N_gridp,1);  % standard dev
end

% degree day factor, global 
if strcmp(Dday_factor_type, 'global') == 1
    Dday_factor=exp(log(5)+1.*randn(1,Ne));
    Dday_factor_DDMinput = repmat(Dday_factor, N_gridp, 1);

% degree day factor, local    
elseif strcmp(Dday_factor_type, 'local') == 1
    if corr_prior_distr == 0 % option A: No spatial correlation
        Dday_factor=exp(log(5)+1.*randn(N_gridp,Ne));
        Dday_factor_DDMinput = Dday_factor;
    elseif corr_prior_distr == 1 % option B: Possible spatial correlation using the CGS function (change c = cut off distance)
        Dday_factor=exp(CGS(rho,log(5).*ones(N_gridp,1),sig,Ne));
        Dday_factor_DDMinput = Dday_factor;
    else
        disp('ERROR: Dday_factor not defined because ''corr_prior_distr'' is not either 1 or 0');
    end
    
else
    disp('ERROR: Dday_factor not defined because ''Dday_factor_type'' is either ''local'' nor ''global''')
end

% precipitation factor
if corr_prior_distr == 0    % option A: prior distribution uncorrelated between grid points     
    P_factor = exp(log(1)+ 1.* randn(N_gridp,Ne));   
elseif corr_prior_distr == 1 % optionB: spatial correlation between grid points using the CGS function
    P_factor = exp(CGS(rho,log(1).*ones(N_gridp,1),sig,Ne));
else
    disp('ERROR: P_factor not defined because ''corr_prior_distr'' is not either 1 or 0');
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
        Yp = zeros(No_tot, Ne);
        tmp2 = snowdepth(obs_times_ens); % time, space, ensemble; tmp lists all model predicted observations for all ensemble members in a [No_tot * Ne,1] array
        % sort [NoxNe,1] listed model predicted observations in a [No,Ne] matrix 
        % results in a matrix that lists all pred obs for gridp 1, followed by those for gridp2 etc, in rows with Ne colomns for all ensemble members   
        jj1 = 1;
        for ii=1:Ne
            jj2 = jj1+No_tot-1;
            Yp(:,ii) = tmp2(jj1:jj2);  % [Nobs x Ne] observation matrix with on each row the ensemble prediction for each observation
            jj1 = jj2+1;
        end
        % define the observation matrix
        y=snowdepth_obs(~isnan(snowdepth_obs)); % creates a [Nobs, 1] list of all observation with first all obs for gridp 1, followed by all obs for gridp 2, etc til gridp N_gridp
        % define observation covariance (measure for observation uncertainty)
        R=y_std.^2;
        alpha=Na; % number of MDA iterations
        pert_stat=1; % 1 for MDA, 0 for ES
        

        
        if  cov_localization == 0
            % update parameters with the EnKA function
            phi=EnKA(phi,Yp,y,R,alpha,pert_stat);
        elseif cov_localization == 1
            % define location matrix for the predicted observations
            loc_Yp  = loc_obs(~isnan(snowdepth_obs));
            % define location matrix for the parameters
            loc_theta = nan(size(theta,1),1);
            N_ddf = size(Dday_factor,1);
            N_Pfac = size(P_factor,1);
            for ii = 1:N_ddf
                loc_theta(ii) = ii;
            end
            for ii = 1:N_Pfac
                jj = N_ddf+ii;
                loc_theta(jj) = ii;
            end
            % update parameters with the EnKA_covloc function using covariance localization    
            phi=EnKA_covloc(phi,Yp,y,R,alpha,pert_stat, loc_Yp, loc_theta, c_covloc);
        else
            disp('ERROR: ''cov_localization'' is not either 0 or 1, no EnKA parameter updates');
        end
            
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

fig1 = figure(1); 
clf;
set(fig1, 'position', [75 25 1850 1000])
for loc=1:N_gridp
    subplot(2,N_gridp,loc);
    plot(t,squeeze(snowdepth_pri(:,loc,:)),'LineWidth',0.5,'Color',[0.8 0 0 0.1]); hold on;
    plot(t,squeeze(snowdepth_post(:,loc,:)),'LineWidth',0.5,'Color',[0 0 0.8 0.1]); hold on;
    pt(1)=plot(t,squeeze(mean(snowdepth_pri(:,loc,:),3)),'LineWidth',2,'Color',[0.8 0 0 1]); hold on;
    pt(2)=plot(t,squeeze(mean(snowdepth_post(:,loc,:),3)),'LineWidth',2,'Color',[0 0 0.8 1]); hold on;
    snowdepthtl=snowdepth_true(:,loc);
    pt(3)=plot(t,snowdepthtl,'-k','LineWidth',2);
    pt(4)=scatter(t_obs(:,loc),snowdepth_obs(:,loc),100,'o','MarkerFaceColor',[0.7 0.7 0],'MarkerEdgeColor',[0 0 0]);
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
    
    subplot(2,N_gridp,N_gridp+loc);
    if strcmp(Dday_factor_type, 'global') == 1
        xpri=Dday_factor_pri(1,:); ypri=P_factor_pri(loc,:);
        xpost=Dday_factor_post(1,:); ypost=P_factor_post(loc,:);
        xtrue = Dday_factor_true * ones(N_gridp,1); ytrue = P_factor_true;
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

if save_figure == 1
    % save figure
    print('-djpeg',[fig_path, '/results'],'-r300','-opengl');
end

%% plot histogram P_factor prior and posterior last grid point
% bar_matrix = [ypri; ypost]';
% f2 = figure(2);
% clf
% H = hist(bar_matrix,[0.25:0.25:5]);
% hist(bar_matrix,[0.25:0.25:5]);
% hold on;
% plot([ytrue(end), ytrue(end)], [0,max(max(H))], 'linewidth', 3, 'color', [0,0,0])
% xlim(gca,[0,5]);
% l1 = legend('prior', 'posterior', 'truth');
% title('P factor distribution last grid point')









