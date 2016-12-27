%% computation of Koopman Spectral Properties for periodic lid-driven 
% cavity flow art Re=13k
% "Ergodic Theory, Dynamic Mode Decomposition & computation of Koopman
% spectral properties" by Hassan Arbabi & Igor Mezic 2016
% sec III
clc,clear

% load the data
load cavity13k_p.mat

% t : time
% KE: kinetic energy of the flow
% G : values of stream function measured at some random points in the
% domain
% KnownFreqs: Koopman Frequencies computed from FFT+Harmonic averaging

%% Hankel-DMD method

n = 100;       % # of operator iterations - should be larger than the dimension of invaraint subspace
m = 100;       % # of points on which functions are smapled - should be large enough to sample the whole limit cycle

nd = m+n;
Data1 = KE(1:nd);     % we use kinetic energy - feel free to use any row of G


% the Hankel matrix
c = Data1(1:n).';
r = Data1(n:n+m);
H = hankel(c,r).';

% applying SVD-enhanced DMD
[ HModes,Hank,~ ] = DMD_Schmid(H);

Freqs = ( log(Hank)./(1i*.1));

K = [0,1,2,3,4,6];      % the harmonic numbers - better be figured out manually

LineForm = {'-','--','-.',':','o','s'};



% plot eigenfunctions  along the trajectory

figure(12),clf, hold on
for im = 1:6
    TheMode = real(HModes(:,2*im-1));
    plot((1:m+1),TheMode./max(TheMode),LineForm{im},'LineWidth',1.25)
    Title{im}=['$\omega=',num2str(abs(real(Freqs(2*im-1)))),'$'];
    Modes(im,:)=HModes(:,2*im-1);
end
xlim([1,100])
ylim([-1.2 1.2])
legend(Title,'interpreter','latex','FontSize',12,'Location','SouthWest')
box on


