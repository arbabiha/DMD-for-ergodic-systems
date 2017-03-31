%% computing the asymptotic phase for trajectories of Van der Pol 
% "Ergodic Theory, Dynamic Mode Decomposition & computation of Koopman
% spectral properties" by Hassan Arbabi & Igor Mezic 2016
% sec III.B
clc,clear

% the Van der Pol oscillator model
mu = 0.3;
VDP = @(x) [x(end/2+1:end); ...
                       mu*(1-x(1:end/2).^2).*x(end/2+1:end)-x(1:end/2)];
                   
IC = -4*[1;1];   % initial condition
dt =0.1;        % time steps
tspan = 0:dt:100;
[T,Z]= ode45(@(t,x)VDP(x),tspan,IC);

%% we choose the observable f=z1+z2 and put it in a row vector
Data = (Z(:,1)+Z(:,2)).'; 

%  Hankel-DMD method
m = 350;       % number of points on which function is smapled 
n = 100;    % number of Koopman operator iterations 


[ HModes, Evalues, Norms ] = DMD.Hankel_DMD( Data,n,m );


w0_Hankel = ( log(Evalues(1))./(1i*dt));     % the dominant frequency - w_0

BETA = angle(HModes(:,1));

%%  plotting the phase on the trajectory
x = Z(1:m,1).'; y = Z(1:m,2).';
z = zeros(size(x));
col = BETA';

figure(21),clf
subplot(1,2,2)
surface([x;x],[y;y],[z;z],[col;col],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',2);
colormap('jet')
cb=colorbar; ylabel(cb,'$~~\theta$','interpreter','latex','FontSize',16,'rot',0)
title(['$\omega_0$=',num2str(abs(real(w0_Hankel)))],'interpreter','latex','FontSize',12)
% axis([-3,5,-3,5])
box on
set(gca,'xtick',[-3:2:5],'ytick',[-3:2:5])
xlabel('$z_1$','interpreter','latex','FontSize',12),ylabel('$z_2$','interpreter','latex','FontSize',14)

%% comparison with the harmonic average
% first compute the eigenfunction using harmonic averaging
w0 = 0.995;     % the basic frequency in harmonic averaging should be known a priori
tdata = tspan(1:m);
fstar = HarmonicAverage(Data(1:m),w0,tdata);     % this gives the eigenfunction at initial condition
fstar = exp(1i*w0*(tdata))*fstar;          % -- extending to the whole trajectory using the fact phi(t)=exp(iw_0t)phi(0)

% adjust the initial phase 
Ratio = fstar(1)/HModes(1,1);   
Phi0 = HModes(:,1)*Ratio;
    
figure(21),
subplot(1,2,1)
cla
plot(tdata,real(fstar),'o')
hold on
plot(tdata,real(Phi0))
box on
set(gca,'xtick',[0,max(tdata)],'XTickLabel',{'0',num2str(max(tdata))},'FontSize',10)
title('real($\phi_0$)','interpreter','latex','FontSize',12)
xlabel('$t~(sec)$','interpreter','latex','FontSize',12)
legend({'Fourier average','Hankel-DMD'},'interpreter','latex','FontSize',11)
set(gcf,'Position',[200 100  1100 420])
 