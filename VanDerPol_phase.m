%% computing the asymptotic phase for trajectories of Van der Pol 
% "Ergodic Theory, Dynamic Mode Decomposition & computation of Koopman
% spectral properties" by Hassan Arbabi & Igor Mezic 2016
% sec III.B
clc,clear

% the Van der Pol oscillator model
mu = 0.3;
VDP = @(x) [x(end/2+1:end); ...
                       mu*(1-x(1:end/2).^2).*x(end/2+1:end)-x(1:end/2)];
                   
IC = 4*[1;1];   % initial condition
dt =0.1;        % time steps
tspan = 0:dt:200;
[T,Z1]= ode45(@(t,x)VDP(x),tspan,IC);
IC = [3,4];   % initial condition
[~,Z2]= ode45(@(t,x)VDP(x),tspan,IC);

%% we choose the observable f=z1+z2 and put it in a row vector
Data1 = (Z1(:,1)+Z1(:,2)).'; 
Data2 = (Z2(:,1)+Z2(:,2)).'; 
%  Hankel-DMD method
m = 250;       % number of points on which function is smapled 
n = 100;    % number of Koopman operator iterations 
index1 = 1:n;
index2 = n:n+m-1;
    % Hankel blocks ()
    c = Data1(index1).'; r = Data1(index2);
    H = hankel(c,r).';
    c = Data1(index1+1).'; r = Data1(index2+1);
    UH= hankel(c,r).';
    % Hankel blocks () on the second trajectory
    c = Data2(index1).'; r = Data2(index2);
    H2 = hankel(c,r).';
    c = Data2(index1+1).'; r = Data2(index2+1);
    UH2= hankel(c,r).';
    
[ HModes, Evalues, ExactModes,Norms ] = DMD.Exact_DMD( [H;H2],[UH;UH2],1e-10 );

iM = 1;     % the basic frequency should be picked  manually: depending on the observable f it might not be the first 
w0_VanderPol = ( log(Evalues(iM))./(1i*dt));     % the basic frequency - w_0
ws = ( log(Evalues(1:end))./(1i*dt));
BETA1 = angle(HModes(1:m,iM));
BETA2 = angle(HModes(m+(1:m),iM));


%%  plotting the phase on the trajectory
x = Z1(1:m,1).'; y = Z1(1:m,2).';
z = zeros(size(x));
col = BETA1';
set(0,'defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex');
figure(21),clf,set(gcf,'Position',[100 100  850 300])

subplot(1,2,2)
plot([Z1(1,1),Z2(1,1)],[Z1(1,2),Z2(1,2)],'ok','MarkerSize',5,'MarkerFaceColor','k')
hold on
surface([x;x],[y;y],[z;z],[col;col],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',2);
colormap('jet')
cb=colorbar; ylabel(cb,'~~$\theta$','interpreter','latex','FontSize',16,'rot',0)
title(['$\omega_0$=',num2str(abs(real(w0_VanderPol)))],'FontSize',14)
box on
ylim([-3 5]),xlim([-4,5])
set(gca,'xtick',[-3:2:5],'ytick',[-3:2:5])
xlabel('$z_1$','FontSize',14),ylabel('$z_2$','FontSize',14)



x = Z2(1:m,1).'; y = Z2(1:m,2).';
z = zeros(size(x));
col = BETA2';
surface([x;x],[y;y],[z;z],[col;col],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',2);

%% comparison with the harmonic average
% first compute the eigenfunction using harmonic averaging
w0 = 0.995;     % the basic frequency in harmonic averaging should be known a priori 
tdata = tspan(1:m);
fstar = DMD.HarmonicAverage(Data1(1:m),w0,tdata);     % this gives the eigenfunction at initial condition 
fstar = exp(1i*w0*(tdata))*fstar;          % -- extending to the whole trajectory using the fact phi(t)=exp(iw_0t)phi(0) -- see ref [23] in the paper

% adjust the initial phase 
Ratio = fstar(1)/HModes(1,iM);   
Phi0 = HModes(1:end/2,iM)*Ratio;
figure(21)
subplot(1,2,1)
cla
plot(tdata,real(fstar),'o')
hold on
plot(tdata,real(Phi0))
box on
ylim(3.5*[-1 1])
set(gca,'xtick',[0,max(tdata)],'XTickLabel',{'0',num2str(max(tdata))},'FontSize',10)
title('$real(\phi_0$) along the trajectory starting at $\mathbf{z}^1$','FontSize',15)
xlabel('$t$','FontSize',14)
legend({'Fourier average','Hankel-DMD'},'interpreter','latex','FontSize',13,'Location','Southwest')




