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
% domain
% KnownFreqs: Koopman Frequencies computed from FFT+Harmonic averaging

%% analysis of kinetic energy

m = 6000;       % # of points on which functions are smapled - should be large enough to sample the whole attractor
n = 500;        % # of operator iterations - hould be larger than the dimension of invaraint subspace
nd = m+n;     % total number of data points

index1 = 1:n;
index2 = n:n+m;


% we use Kinetic Energy and 4th row of G - feel free to use other rows
% we form 4 Hankel matrices

% kinetic energy 
c = KE(index1).'; r = KE(index2);
H1 = hankel(c,r).';

% G(4,:) 
c = G(4,index1).'; r = G(4,index2);
H2 = hankel(c,r).';

% kinetic energy - time-shifted forward
c = KE(index1+1).'; r = KE(index2+1);
UH1 = hankel(c,r).';

% G(4,:) - time-shifted forward
c = G(4,index1+1).'; r = G(4,index2+1);
UH2 = hankel(c,r).';

% normalizing the observables
a12 = norm(H1(:,1))/norm(H2(:,1));

X = [H1 , -a12 * H2];      % data matrices for Exact DMD
Y = [UH1, -a12 *UH2];


[ProjectedModes,DEv,~,~ ] = DMD.ExactDMD( X,Y,200 );

Freqs = real(10*1i* log(DEv(1:20)));    % turn eigenvalues into frequency


 %% first we parameterize the torus
f1= abs(real(Freqs(9)))     % the basic frequencies - should be picked manually
f2= abs(real(Freqs(3)))

Time = t(1:m+1);

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

IM = [1,3,9,14,15,20];      % choose other eigenfunctions

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
    EF(:,:,im)=EF2; 

end

