%% POD of the observables on the Lorenz attractor
% "Ergodic Theory, Dynamic Mode Decomposition & computation of Koopman
% spectral properties" by Hassan Arbabi & Igor Mezic 2016
% sec 4

clc,clear
% Lorenz chaotic model 1963
sigma=10;
rho = 28;
beta = 8/3;

Lorenz = @(x) [sigma*(x(2)-x(1)); ...
               x(1)*(rho-x(3))-x(2);...
               x(1)*x(2)-beta*x(3)];
           
tspan = 0:.01:1000;
x0 = [0.1;0;0.1];

[tspan,X]= ode45(@(t,y)Lorenz(y),tspan,x0);


%% SVD on Hankel --> POD basis on the attractor

m = 15000;       % # of points on which functions are smapled - should be large enough to sample the whole limit cycle
n = 101;         % # of operator iterations - should be larger so that transients die

nd = m+n; % number of total observations fo Hankel
ntr = 1000;     % discard 10 sec transient
x = X(ntr+(1:nd),1)';  % taking x as observable

% make the Hankel matrix
c = x(1:n).';
r = x(n:n+m-1);
H = hankel(c,r).';

[U,S,V]=svd(H/sqrt(m),0);
disp('first six singular values')
diag(S(1:6,1:6))    % check out how they vary with changing m
%% POD basis as colorfield on the attractor
set(0,'defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex');

tdata = tspan(ntr+(1:m));
xc = x(1:m);
yc = X(ntr + (1:m),2)';
zc = X(ntr + (1:m),3)';

figure(122),clf
for im = 1:6
        col = U(:,im)';
        subplot(2,3,im)
        surface([xc;xc],[yc;yc],[zc;zc],[col;col],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',2);
        colormap('jet')

    box on
    set(gca,'xtick',[],'ytick',[],'ztick',[])
    colorbar
    max(col);

    view(-6,13)
    title(['$\tilde \psi_',num2str(im),'$'],'FontSize',14)
end
    xlabel('$z_1$')%,ylabel('y'),
    zlabel('$z_3$')
    set(gcf,'Position',[100 100 1044 450])


%% principal coordinates

figure(21),clf
plot(1:n,V(:,1:6))
xlabel('$i$')
xlim([1 100])
legend({'$v_1$','$v_2$','$v_3$','$v_4$','$v_5$','$v_6$'},'interpreter','latex','FontSize',14)




