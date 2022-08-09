%% Danging blobs
% Demonstrates how to generate a samples from spatial field of correlated
% Gaussian random variables.
clearvars;

%% Set up the grid.
del=0.02; % Grid spacing (don't make it too small, then things blow up!)
x=-1:del:1; % x coordinates 
y=x; % y coordinates (make the domain square for simplicity)
[X,Y]=meshgrid(x,y); % Spatial grid.
xg=X(:); yg=Y(:); % Make the grid into a column vector
n=numel(xg); % Total number of points =numel(y)*numel(x)

%% Calculate correlations based on the GC function
c=0.5; % Cut-off distance for the GC function.
xgt=xg'; ygt=yg';
d=sqrt((repmat(xg,1,n)-repmat(xgt,n,1)).^2+...
    (repmat(yg,1,n)-repmat(ygt,n,1)).^2);
rho=GC(d,c);
% Note, the above does the same as the loop:
%{
for j=1:n
    d=sqrt((xg-xg(j)).^2+(yg-yg(j)).^2);
    rhoj=GC(d,c);
    rho(:,j)=rhoj;
end
%}

%% Generate some (pseudo)random samples
mu=zeros(n,1); % Mean vector
sig=ones(n,1); % Standard deviation vector
N=1e2;
zs=CGS(rho,mu,sig,N);


%% Visualize some dancing blobs
figure('units','normalized','outerposition',[0 0 1 1],...
    'visible','on')
nsqrt=sqrt(n); % We assume the grid is square of size sqrt(n)*sqrt(n)
fs=24;
for j=1:N
    zj=zs(:,j);
    zj=reshape(zj,nsqrt,nsqrt);
    imagesc(x,y,zj);
    xlabel('$x$','Interpreter','Latex','FontSize',fs);
    ylabel('$y$','Interpreter','Latex','FontSize',fs);
    title(sprintf('Ensemble member $%d$',j),...
        'Interpreter','Latex','FontSize',fs);
    axis xy;
    caxis([-2.*max(sig(:)) 3.*max(sig(:))]);
    c=colorbar;
    c.FontSize=fs;
    c.TickLabelInterpreter='latex';
    c.TickDirection='out';
    c.LineWidth=2;
    xlabel(c,'$z$','Interpreter','Latex','FontSize',fs);
    ax=gca; ax.LineWidth = 2;
    axis square;
    box on;
    set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.005, 0.005]);
    set(groot, 'defaultAxesTickLabelInterpreter','Latex');
    set(gca,'XTick',[],'YTick',[],'TickLabelInterpreter','latex');
    k=waitforbuttonpress();
end