function snowdepth=DDM(t,P,T,Dday_factor,P_factor)
%% DDM: A simple degree day snowmelt model
% 

% Hyperparameters:
Tsnow=1; % Temperature threshold for snowfall (degrees C)
Tmelt=0; % Temperature threshold for snowmelt (degrees C)

Nt=numel(t); % Number of time steps.
Ns=size(P,2); % Number of points in space (can be just 1).
Ne=size(Dday_factor,2); % Number of ensemble members (can be just 1).
snowdepth_old=0; % Initial condition (no snow).

% Ensure arrays are compatible
Pe=repmat(P,1,1,Ne);
Te=repmat(T,1,1,Ne);
if Ne==1&&(size(Dday_factor,1)==Ns) 
    % In this case (deterministic run=single ensemble member) 
    % Pj and Tj will be 1 x Ns, so ensure Dday_factor, P_factor are 1 x Ns for
    % compatibility in the elementwise products for the Ablation and 
    % Accumulation terms.
    Dday_factor=Dday_factor'; 
    P_factor=P_factor';
elseif Ne>1&&Ns>1 
    % In this case Pj and Tj will be 1 x Ns x Ne
    % Make Dday_factor , P_factor size 1 x Ns x Ne
    tmp=zeros(1,Ns,Ne);
    tmp(1,:,:)=Dday_factor;
    Dday_factor=tmp;
    tmp(1,:,:)=P_factor;
    P_factor=tmp;
end
%{
if Ns>1&&Ne>1 % If many points in space and ensemble.
    Pe=repmat(P,1,1,Ne);
    Te=repmat(T,1,1,Ne);
elseif Ns==1&&Ne>1 % If just ensemble.
    Pe=repmat(P,1,Ne);
    Te=repmat(T,1,Ne);
elseif Ne==1 % If deterministic (single realizations).
    Pe=P;
    Te=T;
end
%}
snowdepth=zeros(size(Pe));

for j=1:Nt
    %Pj=squeeze(Pe(j,:,:)); % Precip for this time step
    %Tj=squeeze(Te(j,:,:)); % Air temperature for this time step
    Pj=Pe(j,:,:);
    Tj=Te(j,:,:);
    cansnow=Tj<=Tsnow; % Snow is possible.
    P_factor_j=P_factor;
    P_factor_j(~cansnow)=0;
    try
        Accumulation=P_factor_j.*Pj;
    catch
        disp('hi! error in accumulation calculation');
    end
    dd=Tj-Tmelt; % Degree day for this day.
    melting=dd>0;
    Dday_factor_j=Dday_factor;
    Dday_factor_j(~melting)=0;
    Ablation=Dday_factor_j.*dd;
    snowdepth_new=snowdepth_old+Accumulation-Ablation;
    snowdepth_new=max(snowdepth_new,0);
    snowdepth(j,:,:)=snowdepth_new; 
    snowdepth_old=snowdepth_new;
end



end