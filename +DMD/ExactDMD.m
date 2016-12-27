function [ ProjectedModes,DEv,ExactModes,Norm ] = ExactDMD( X,Y,k )
% Dynamic Mode Decomposition as presented by 
% "On Dynamic Mode Decomposition: theory and applications" by Tu et al.,
% arXiv, 2013



% inputs : 
% Data Sets X and Y- should have the same size 
% each column of X is a set of measurements done at an instant 
% each column of Y is the image of the corresponding columns from X
% k - number of basis elemnts used to realize DMD operator

% outputs: 
% 1 - Projected Dynamic Modes
% 2 - Exact Dynamic Modes
% 3 - Dynamic Eigenvalues
% 4 - Norms - Euclidean (vector) norm of each mode, used to sort the data



[U,S,V]=svd(X,'econ');

U = U(:,1:k); V=V(:,1:k); S=S(1:k,1:k);
disp(['DMD subspace dimension:',num2str(k)])

Atilde = ((U'*Y) *V )* diag((1./diag(S)));


[W,DEv]=eig(Atilde);
DEv = diag(DEv);

ProjectedModes = U*W;

ExactModes = bsxfun(@times,1./DEv.',((Y*V)*S^(-1))*W);



if nargout>3
   b = pinv(ProjectedModes)*X(:,1);     % terminal coordinates in the Koopman subspace
    [Norm,Index]=sort(abs(b),'descend');
    DEv = DEv(Index);
    ProjectedModes = ProjectedModes(:,Index);
    disp('modes sorted based on energy contribution to the last snapshot')
end

end


%=========================================================================%
% Hassan Arbabi - 08-17-2015
% Mezic research group
% UC Santa Barbara
% arbabiha@gmail.com
%=========================================================================%