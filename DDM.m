function D=DDM(t,P,T,a,b)
%% DDM: A simple degree day snowmelt model
% 

% Hyperparameters:
Tsnow=1; % Temperature threshold for snowfall (degrees C)
Tmelt=0; % Temperature threshold for snowmelt (degrees C)

Nt=numel(t); % Number of time steps.
Ns=size(P,2); % Number of points in space (can be just 1).
Ne=size(a,2); % Number of ensemble members (can be just 1).
Dold=0; % Initial condition (no snow).

% Ensure arrays are compatible
Pe=repmat(P,1,1,Ne);
Te=repmat(T,1,1,Ne);
if Ne==1&&(size(a,1)==Ns) 
    % In this case (deterministic run=single ensemble member) 
    % Pj and Tj will be 1 x Ns, so ensure a, b are 1 x Ns for
    % compatibility in the elementwise products for the Ablation and 
    % Accumulation terms.
    a=a'; 
    b=b';
elseif Ne>1&&Ns>1 
    % In this case Pj and Tj will be 1 x Ns x Ne
    % Make a , b size 1 x Ns x Ne
    tmp=zeros(1,Ns,Ne);
    tmp(1,:,:)=a;
    a=tmp;
    tmp(1,:,:)=b;
    b=tmp;
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
D=zeros(size(Pe));

for j=1:Nt
    %Pj=squeeze(Pe(j,:,:)); % Precip for this time step
    %Tj=squeeze(Te(j,:,:)); % Air temperature for this time step
    Pj=Pe(j,:,:);
    Tj=Te(j,:,:);
    cansnow=Tj<=Tsnow; % Snow is possible.
    bj=b;
    bj(~cansnow)=0;
    try
        Accumulation=bj.*Pj;
    catch
        disp('hi!');
    end
    dd=Tj-Tmelt; % Degree day for this day.
    melting=dd>0;
    aj=a;
    aj(~melting)=0;
    Ablation=aj.*dd;
    Dnew=Dold+Accumulation-Ablation;
    Dnew=max(Dnew,0);
    D(j,:,:)=Dnew; 
    Dold=Dnew;
end



end