%% Example doing spatially distributed DA for a larger area around Finse
clearvars;


t_start     = '01-Sep-2018';
t_end       = '01-Sep-2019';
t=datenum(t_start):datenum(t_end); t=t';

%% Read in the forcing

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
load('input/TopoSCALE_spatial.mat');


%% Set up the domain
x=f.tps.xg;
y=f.tps.yg;
cn=f.tps.cn; % Cluster numbers of each grid cell from TopoSUB
z=f.tps.z;
z=z(cn);
Z=nan(size(f.tps.mask));
Z(f.tps.mask)=z;
figure(1); clf;
imagesc(x,y,Z);
axis xy;

x0=413782+2e3; y0=6724560-2e3; 
xl=[x0 x0+8e3];
yl=[y0-8e3 y0]; % 100 x 100 grid cells (of 100 m resolution) around Finse.

xlim(xl);
ylim(yl);

thesex=x>=min(xl)&x<=max(xl);
thesey=y>=min(yl)&y<=max(yl);
xg=x(thesex);
yg=y(thesey);
maskg=f.tps.mask(thesey,thesex);

% Set the spatial boundary

%% Generate forcing 

f.P=double(f.P).*f.P_sf; % unpack precipitation data
f.T=double(f.T).*f.T_sf; % unpack temperature data

% Translate hourly values into daily values
Nt=numel(t);
Ns=size(f.P,1);
P=zeros(Nt,Ns);
T=P;
for j=1:Nt
   here=f.t>t(j)&f.t<=(t(j)+1);
   T(j,:)=mean(f.T(:,here),2);
   P(j,:)=sum(f.P(:,here),2); 
end



