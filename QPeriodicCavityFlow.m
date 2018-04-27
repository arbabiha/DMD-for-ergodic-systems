%% Koopman Spectral Properties for quasi-periodic lid-driven cavity flow
% at Re=16k
% "Computation of Koopman spectrum for systems with ergodic attractors" 
% by Hassan Arbabi & Igor Mezic 2016
% section 6.3
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


[ ProjectedModes, DEv, Norms ] = DMD.Hankel_DMD( Data([1,5],:) ,n,m );      % we choos two observables = two rows of Data matrix

Freqs = real(10*1i* log(DEv(1:20)));    % turn eigenvalues into frequency

 %% first we parameterize the torus
 disp('basic frequencies:')
f1= abs(real(Freqs(9)))     % the basic frequencies - should be picked manually
f2= abs(real(Freqs(3)))

Time = t(1:m);

xtorus = mod(f1*Time,2*pi);
ytorus = mod(f2*Time,2*pi);
 
 % visualize trajectory unfolded on a torus
 set(0,'defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex');
figure(6),clf
plot(xtorus(1:end),ytorus(1:end),'.')
axis square, axis(2*pi*[0 1 0 1])
xlabel('$\theta$'), ylabel('$\eta$')

title({'trajectory on torus:','this plot can be used to see if sampling is dense enough'},'FontSize',14)
%% plotting the eigenfunctions on the parameterized torus
% we have the value on the trajectory 
% we interpolate it onto the whole torus
[X0,Y0]=meshgrid(0:.05:2*pi);
figure(36),title('Koopman eigenfunctions')
EF=zeros([size(X0),6]);

IM = [1,2,8,12,14,16];      % this choice of modes gives figure 7 in the paper - feel free to check other modes

for im = 1:6
    Mode = ProjectedModes(:,IM(im))./max(abs(real(ProjectedModes(:,IM(im)))));  % normalizing eigenfunctions
    subplot(2,3,im)
    FInterp = scatteredInterpolant(xtorus,ytorus,Mode,'linear'); 
    EF2 = FInterp(X0,Y0);
    contourf(X0,Y0,real(EF2),50,'LineStyle','None')
    colorbar, caxis([-1 1]), colormap('jet')
    set(gca,'xtick',[0 2*pi],'ytick',[0 2*pi])
    set(gca,'XTickLabel',{'',''},'YTickLabel',{'',''},'FontSize',14);  

end
% adding title and label 
subplot(2,3,2), 
xlabel('$\theta$'),ylabel('$\eta$'), title('Koopman eigenfunctions')
set(gcf,'Position',[100 100 890 395])

