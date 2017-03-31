%% Koopman Spectral Properties for quasi-periodic lid-driven cavity flow
% at Re=16k
% "Ergodic Theory, Dynamic Mode Decomposition & computation of Koopman
% spectral properties" by Hassan Arbabi & Igor Mezic 2016
% sec V
clc,clear

% load the data
load cavity16k_qp.mat

% t : time
% KE: kinetic energy of the flow
% G : values of stream function measured at some random points in the
% flow domain 
% KnownFreqs: Koopman Frequencies computed from FFT+Harmonic averaging

%% Hankel-DMD 

m = 6000;       % # of points on which functions are smapled - should be large enough to sample the whole attractor
n = 500;        % # of operator iterations - hould be larger than the dimension of invaraint subspace

Data = [KE;G(:,1:length(KE))];  % the data matrix that goes into Hankel-DMD: each row is measurements on an observable


[ ProjectedModes, DEv, Norms ] = DMD.Hankel_DMD( Data([1,5],:) ,n,m );

Freqs = real(10*1i* log(DEv(1:20)));    % turn eigenvalues into frequency

 %% first we parameterize the torus
f1= abs(real(Freqs(9)))     % the basic frequencies - should be picked manually
f2= abs(real(Freqs(3)))

Time = t(1:m);

xtorus = mod(f1*Time,2*pi);
ytorus = mod(f2*Time,2*pi);
 
 % visualize trajectory unfolded on a torus
figure(99),clf
plot(xtorus(1:end),ytorus(1:end),'.')
hold on
plot(xtorus(1),ytorus(1),'rx')
plot(xtorus(end),ytorus(end),'kx')
axis square
axis(2*pi*[0 1 0 1])
set(gca,'xtick',[0 2*pi],'ytick',[0 2*pi])
set(gca,'TickLabelInterpreter', 'latex','XTickLabel',{'0','$2\pi$'},'YTickLabel',{'0','$2\pi$'},'FontSize',14);
xlabel('$\theta=\omega_1 t$','interpreter','latex')
ylabel('$\phi=\omega_2 t$','interpreter','latex')
title({'trajectory on torus.',' this plot can be used to see if sampling is dense enough'},'FontSize',10)
%% plotting the eigenfunctions on the unfolding of the torus
% we have the value on the trajectory 
% we interpolate it onto the whole torus

[X0,Y0]=meshgrid(0:.05:2*pi);
figure(102),clf
EF=zeros([size(X0),6]);

IM = [1,2,4,6,8,10];      % efunctions come in complex conjugate pairs except the almost constant efunction

for im = 1:6
    Mode = ProjectedModes(:,IM(im))./max(abs(real(ProjectedModes(:,IM(im)))));  % normalizing eigenfunctions
    subplot(2,3,im)
    FInterp = scatteredInterpolant(xtorus,ytorus,Mode,'linear'); 
    EF2 = FInterp(X0,Y0);
    contourf(X0,Y0,real(EF2),50,'LineStyle','None')
    colorbar
    caxis([-1 1])
    colormap('jet')
    set(gca,'xtick',[0 2*pi],'ytick',[0 2*pi])
    set(gca,'TickLabelInterpreter', 'latex','XTickLabel',{'',''},'YTickLabel',{'',''},'FontSize',14);
    xlabel('$\theta$','interpreter','latex')
    ylabel('$\phi$','interpreter','latex')    

end

