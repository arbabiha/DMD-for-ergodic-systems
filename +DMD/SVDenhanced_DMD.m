function [ DModes,DEv,Norm ] = SVDenhanced_DMD( Data,varargin )
% Dynamic Mode Decomposition as presented by 
% "Dynamic Mode Decomposition of numerical and experimental data" by Peter
% J. Schmid, Journal of Fluid Mechanics, 2010

% inputs : 
% Data - the data matrix: each column is a a set of measurements done at
% each instant - the sampling frequency is assumed to be constant

% or 

% ( Data,Tol ) - with Tol (optional) being the threshold for filtering thru SVD - the
% default value is 1e-10

% outputs: 
% 1 - DModes- Koopman or Dynamic Modes: structures that scale exponentially
% with time - the exponent is given by Koopman or Dynamic Eigenvalues and
% could be imaginary as well

% 2- DEv - Koopman or Dynamic Eigenvalues: the exponents and frequencies
% that make up the time-varaiation of data in time

% 3- Norm - The norm here is set to the energy contribution of each mode to
% the last snapshot of data and then used to sort the modes




tic
if isempty(varargin)
    Tol=1e-10;
else
    Tol=varargin{1};
end

disp(['Tolerance used for filtering in DMD:',num2str(Tol)])

X=Data(:,1:end-1);
Y=Data(:,2:end);


[U,S,V]=svd(X,0);
r = find(diag(S)>Tol,1,'last');
disp(['DMD subspace dimension:',num2str(r)])

U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);

Atilde = (U'*(Y*V)) / S; 

[w,lambda]=eig(Atilde);

DModes = U*w;
DEv = diag(lambda);

if nargout>2
    b = pinv(DModes)*Data(:,end);     % terminal coordinates in the Koopman subspace
    [Norm,Index]=sort(abs(b),'descend');
    DEv = DEv(Index);
    DModes = DModes(:,Index);
    disp('modes sorted based on energy contribution to the last snapshot')
end

disp(['time to compute DMD:',num2str(toc)])

end

%=========================================================================%
% Hassan Arbabi - 08-17-2015
% Mezic research group
% UC Santa Barbara
% arbabiha@gmail.com
%=========================================================================%
