% Find a better directory name and comment this a bit better.
clearvars;
load('TopoSCALE.mat');
f.P=double(f.P).*f.P_sf;
f.T=double(f.T).*f.T_sf;

ysd=20;

% ddf for snow 2.5 to 11.6 mm/d/K
t=datenum('01-Sep-2018'):1:datenum('01-Sep-2019'); t=t';
Nt=numel(t);
P=zeros(Nt,2);
T=P;
for j=1:Nt
   here=f.t>t(j)&f.t<=(t(j)+1);
   T(j,:)=mean(f.T(:,here),2);
   P(j,:)=sum(f.P(:,here),2); 
end

atrue=[6; 6];%6; % mm/d/K
btrue=[0.5; 2];

%D=zeros(size(P));
D=DDM(t,P,T,atrue,btrue);
Dtrue=D;
Dobs=Dtrue+ysd.*randn(size(Dtrue));
Dobs=max(Dobs,0);
tobs=[datenum('01-Jan-2019'); datenum('01-Feb-2019');...%];%;...
    datenum('01-Mar-2019'); datenum('01-Apr-2019');...
    datenum('01-May-2019'); datenum('01-Jun-2019');...
    datenum('01-Jul-2019'); datenum('01-Aug-2019')];
No=numel(tobs);
tmp=zeros(No,2);
here=zeros(Nt,1,'logical');
for j=1:No
    herej=t==tobs(j);
    here=here+herej;
    tmp(j,:)=Dobs(herej,:);
end
here=logical(here);
Dobs=tmp;

%{
figure(1); clf;
plot(t,D);
datetick('x','keepticks','keeplimits');
hold on;
scatter(tobs,Dobs,150,'.');
%}

%% Ensemble DA
Ne=1e2; Na=4;
a=exp(log(5)+1.*randn(2,Ne));
b=exp(log(1)+1.*randn(2,Ne));

k=0;
for ell=1:(Na+1)
    D=DDM(t,P,T,a,b);
    theta=[a;b];
    if ell==1
        thetapri=theta;
        apri=a(1:2,:);
        bpri=b(1:2,:);
        Dpri=D;
    elseif ell==(Na+1)
        thetapost=theta;
        apost=a(1:2,:);
        bpost=b(1:2,:);
        Dpost=D;
    end
    if ell<=Na
        phi=log(theta);
        Yp=D(here,:,:); % time, space, ensemble
        Yp=reshape(Yp,No*2,Ne);
        y=Dobs;
        y=y(:);
        R=ysd.^2;
        alpha=Na;
        pert_stat=1;
        phi=fastpEnKF(phi,Yp,y,R,alpha,pert_stat);
        theta=exp(phi);
        a=theta(1:2,:);
        b=theta(3:4,:);
    end
end

return

%%
figure(1); clf;

for loc=1:2
    subplot(2,2,loc);
    plot(t,squeeze(Dpri(:,loc,:)),'LineWidth',0.5,'Color',[0.8 0 0 0.1]); hold on;
    plot(t,squeeze(Dpost(:,loc,:)),'LineWidth',0.5,'Color',[0 0 0.8 0.1]); hold on;
    pt(1)=plot(t,squeeze(mean(Dpri(:,loc,:),3)),'LineWidth',2,'Color',[0.8 0 0 1]); hold on;
    pt(2)=plot(t,squeeze(mean(Dpost(:,loc,:),3)),'LineWidth',2,'Color',[0 0 0.8 1]); hold on;
    Dtl=Dtrue(:,loc);
    pt(3)=plot(t,Dtl,'-k','LineWidth',2);
    pt(4)=scatter(tobs,Dobs(:,loc),100,'o','MarkerFaceColor',[0.7 0.7 0],'MarkerEdgeColor',[0 0 0]);
    %scatter(tobs,Dobs(:,loc),250,'.','MarkerEdgeColor',[0 0 0]);
    leg=legend(pt(:),{'Prior','Post','Truth','Obs'},'Interpreter','Latex',...
        'FontSize',16,'Location','NorthWest');
    clear pt;
    yl=5*max(Dtl); 
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
    xpri=apri(loc,:); ypri=bpri(loc,:);
    alphis=0.4;
    pt(1)=scatter(xpri,ypri,4e1,'o','MarkerFaceColor',[0.8 0 0],'MarkerFaceAlpha',alphis,...
        'MarkerEdgeColor',[0 0 0]); hold on;
    xpost=apost(loc,:); ypost=bpost(loc,:);
    pt(2)=scatter(xpost,ypost,4e1,'s','MarkerFaceColor',[0 0 0.8],'MarkerFaceAlpha',alphis,...
        'MarkerEdgeColor',[0 0 0]); 
    pt(3)=scatter(atrue(loc),btrue(loc),250,'p','MarkerEdgeColor',[0 0 0],'LineWidth',2,...
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

print('-djpeg','results','-r300','-opengl');