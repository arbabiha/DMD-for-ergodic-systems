function [ Average ] = HarmonicAverage( Data,w,t,Filter)
% computes the harmonic average using the Filter

HarmonicWeight = 2*exp(-1i*w'*t);            

% using the specified filter
if exist('Filter','var')
   switch Filter
       case 'Hamming'
           HarmonicWeight=bsxfun(@times,HarmonicWeight,hamming(length(t))');
       case 'Hann'
           HarmonicWeight=bsxfun(@times,HarmonicWeight,hann(length(t))');
       case 'exponential'       % the exponential weighting defined by Das & Yorke 2015
           HarmonicWeight=bsxfun(@times,HarmonicWeight,ExpWeight(length(t))*length(t));
       otherwise
           disp('filter not defined ... no filter used')
   end 
end

Average = HarmonicWeight*Data'./length(t);     

end