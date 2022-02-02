function rho=GC(d,c)
%% rho=GC(d,c): Gaspari-Cohn correlation function
% The fifth-order piecewise-rational function fo Gaspari and Cohn (1999) as
% expressed in eq. 16 of Rho et al. (2016). 
% Note that this routine has been vectorized, so it accepts inputs in the
% form of separation distances d as arrays (of any dimensions).
%
% Inputs:
% d : Absolute separation distance between points (grid cells or time steps).
% c : Localization constant (decay length-scale, at 2*c rho is 0).
%
% Example:
% d=0:0.01:5; 
% c=2; % Correlation goes to zero at 2*c=4
% rho=GC(d,c);
% plot(d,rho); 
%
% References:
% Gaspari and Cohn (199): QJRMS, https://doi.org/10.1002/qj.49712555417
% Roh et al. (2015): NPG, https://doi.org/10.5194/npg-22-723-2015


rho=zeros(size(d));
near=d<c;
mid=(d>=c)&(d<=2*c);
far=d>2*c;

% Nearby points
dn=d(near); 
dnc=dn./c;
rn=(-1/4).*(dnc).^5+(1/2).*(dnc).^4+(5/8).*(dnc).^3-(5/3).*(dnc).^2+1;
rho(near)=rn;

% Mid-range points
dm=d(mid);
dmc=dm./c;
rm=(1/12).*(dmc).^5-(1/2).*(dmc).^4+(5/8).*(dmc).^3+(5/3).*(dmc).^2-5.*(dmc)+4-(2/3).*(c./dm);
rho(mid)=rm;

% Far points
rf=0;
rho(far)=rf;


end