function Img_real = RealDiffusion_lowRes(complex_data,PF,ExtraCropPElines,POCS_flag)

% Kawin Setsompop, Oct 5th 2015

% INPUT: 
%complex_data: this is a 5D vectors (x,y,z,diffusion,pigeon (slider))




%% POCS reconstruction of the partial Fourier data
disp('START: lowRes Real diffussion')

Np = size(complex_data,2);

kdata = mrir_fDFT_freqencode(mrir_fDFT_phasencode(complex_data));
%kdata = mrir_fDFT_freqencode(mrir_fDFT_phasencode(complex_data(:,:,end/2,2,2)));
kdata = kdata(:,[round(Np*(1-PF)+ExtraCropPElines):end],:,:,:);

[Nf, Np_acq, Nslc, Ndiff, Ns] = size(kdata);
recon_hf = zeros(Nf, Np, Nslc, Ndiff, Ns);

tic
% can actually do the whole matrix at once but seem to be slower than
% do this for loop
for kk = 1:Ndiff
    
    for jj = 1:Ns
        recon_hf(:, :, :, kk, jj) = lowRes_PhaseRemoval(kdata(:, :, :, kk, jj), [Nf, Np, Np_acq],POCS_flag);
    end
end
toc



disp('done with Real diffussion')

Img_real = real(recon_hf);

 