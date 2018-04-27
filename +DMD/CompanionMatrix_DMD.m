function [ KModes,KEv,Norms ] = CompanionMatrix_DMD( Data )
% Dynamic Mode Decomposition as presented by 
% "Spectral analysis of nonlinear flows" by Rowley et al., Journal of FLuid
% Mechanics, 2009


% inputs : 
% Data - the data matrix: each column is a a set of measurements done at
% each instant - the sampling frequency is assumed to be constant

% outputs: 
% 1 - KModes- Koopman or Dynamic Modes: structures that scale exponentially
% with time - the exponent is given by Koopman or Dynamic Eigenvalues and
% could be imaginary as well

% 2- KEv - Koopman or Dynamic Eigenvalues: the exponents and frequencies
% that make up the time-varaiation of data in time

% 3- Norms - Euclidean (vector) norm of each mode, used to sort the data




X=Data(:,1:end-1);

c = pinv(X)*Data(:,end);

m = size(Data,2)-1;
C = spdiags(ones(m,1),-1,m,m);
C(:,end)=c;                  %% companion matrix

[P,D]=eig(full(C));                           %% Ritz evalue and evectors

KEunsrtd = diag(D);                     %% Koopman Evalues
KMunsrtd = X*P;                         %% Koopman Modes

Cnorm = sqrt(sum(abs(KMunsrtd).^2,1));  %% Euclidean norm of Koopman Modes
[Norms, ind]=sort(Cnorm,'descend');      %% sorting based on the norm

KEv= KEunsrtd(ind);           %% Koopman Eigenvalues

KModes = KMunsrtd(:,ind);                %% Koopman Modes



end


%=========================================================================%
% Hassan Arbabi - 08-17-2015
% Mezic research group
% UC Santa Barbara
% arbabiha@gmail.com
%=========================================================================%