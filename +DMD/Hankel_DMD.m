function [ HModes, HEvalues, Norms ] = Hankel_DMD( Data,n,m,varargin )
%HANKEL-Dynamic Mode Decompsition as presented in
% "Ergodic theory, dynamic mode decomposition and computation of Koopman
% spectral properties" by H. Arbabi and I. Mezic, arXiv:1611.06664

% inputs : 
% Data - the data matrix: each row is a time series data on an observable
% m,n: number of columns and rows, respectively, in the data Hankel
% matrices: size(Data,2)>m+n and preferably m>>n 

% or 

% ( Data,n,m,tol )- with Tol (optional) being the threshold for filtering thru SVD - the
% default value is 1e-10




% outputs: 
% 1 - HModes- Dynamic modes which approximate the Koopman eigenfunctions

% 2- HEvalues - Koopman or Dynamic Eigenvalues: the exponents and frequencies
% that make up the time-varaiation of data in time

% 3- Norms - The L2-norm of Koopman eigenfunction contribution to the
% observable


% setting the SVD hard threshold for DMD
if isempty(varargin)
    Tol=1e-10;
else
    Tol=varargin{1};
end

index1 = 1:n;
index2 = n:n+m-1;

X = []; Y=[];

for ir = 1:size(Data,1)
   
    % Hankel blocks ()
    c = Data(ir,index1).'; r = Data(ir,index2);
    H = hankel(c,r).';
    c = Data(ir,index1+1).'; r = Data(ir,index2+1);
    UH= hankel(c,r).';
    
    % scaling of Hankel blocks
    if ir>1
        alpha = norm(H(:,1))/norm(X(:,1));
        H = alpha*H;
        UH= alpha* UH;
    end
    
    % the data matrices fed to exact DMD
    X=[X,H]; Y=[Y,UH];
    
end
    % the DMD
   [HModes,HEvalues,~,Norms ] = DMD.Exact_DMD( X,Y,Tol );

    
end

