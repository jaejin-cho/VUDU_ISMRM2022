function [ res ] = mfft( x, dim )
%%%%%%%%%%%%%%%%%%%%%%%% mfft  %%%%%%%%%%%%%%%%%%%%
% made by JaeJin Cho            2016.12.01  
% 
% fft operater
% [ res ] = mfft( DATA, dim )
% DATA    : data
% dim     : dimension number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res = 1/sqrt(size(x,dim))*fftshift(fft(ifftshift(x,dim),[],dim),dim);

end

