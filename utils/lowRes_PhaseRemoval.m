function [recon_hf recon_hf_complex] = lowRes_PhaseRemoval(kdata, Nsize, POCS_flag)

% POCS_flag: if set to 1, then do POCS for p.f., o.w. just zero fill

 Nf     = Nsize(1);
 Np     = Nsize(2);
 Np_acq = Nsize(3);
 
 [Nf Np_acq Nslc Ndiff Ns] = size(kdata);
 
 FilterCropRatio = 1; %1; % can be above one 
 
 %% Phase estimation from the acquired symmetric k-space data
 sym_acq_region  = (Np - Np_acq + 1):Np_acq;
 acq_region      = (Np - Np_acq + 1):Np;
 
 kdata_zf        = zeros(Nf, Np, Nslc, Ndiff, Ns);
 kdata_zf(:,acq_region,:,:,:) = kdata;
    
 hamming_filterPE  = zeros(Nf, Np);
 PE_filterLength = length(sym_acq_region)*FilterCropRatio;
%  PE_filterIndex = [-floor(PE_filterLength/2):floor(PE_filterLength/2)] + round(Np/2);
 PE_filterIndex = [-floor(PE_filterLength/2)+1:floor(PE_filterLength/2)] + round(Np/2);

 hamming_filterPE(:,  PE_filterIndex) = repmat(hamming(numel(PE_filterIndex)).', [Nf, 1]);
%  hamming_filterPE(:,  PE_filterIndex) = repmat(medfilt2(numel(PE_filterIndex)).', [Nf, 1]);
 
 % assume data is same resolution in PE and RO but FOV can be different.- do same amount of Hamming in RO -so get higher SNR phase estimate by
 % smooting this way as well
 hamming_filterRO = zeros(Nf, Np);
 RO_filterLength = length(sym_acq_region)*(Nf/Np)*FilterCropRatio;
%  RO_filterIndex = [-floor(RO_filterLength/2):floor(RO_filterLength/2)] + round(Nf/2);
 RO_filterIndex = [-floor(RO_filterLength/2)+1:floor(RO_filterLength/2)] + round(Nf/2);
 
 hamming_filterRO(RO_filterIndex,:) = repmat(hamming(numel(RO_filterIndex)), [1, Np]);
 
 hamming_filter =  hamming_filterPE.* hamming_filterRO;
 
 image_recon     = mrir_iDFT_freqencode(mrir_iDFT_phasencode(repmat(hamming_filter,[ 1 1 Nslc Ndiff Ns]).*kdata_zf)); 
 phase_recon     = angle(image_recon); % estimating the image phase 
    
 recon_hf     = mrir_iDFT_freqencode(mrir_iDFT_phasencode(kdata_zf)); 
 recon_hf    = recon_hf.*exp(-1i*phase_recon);


 
 if POCS_flag == 1
     recon_last      = recon_hf;
     maxiter = 10; 
     %% POCS reconstruction
     for ii = 1:maxiter

         % project onto the data constraint set
         %kdata_recon = 1/sqrt(Np*Nf)*fftshift(fft2(ifftshift(recon_hf)));
         kdata_recon = mrir_fDFT_freqencode(mrir_fDFT_phasencode(recon_hf));

         kdata_recon(:, acq_region,:,:,:) = kdata;
         %recon_hf    = sqrt(Np*Nf)*ifftshift(ifft2(fftshift(kdata_recon)));
         recon_hf    = mrir_iDFT_freqencode(mrir_iDFT_phasencode(kdata_recon));

         % project onto the phase constraint set
         recon_hf    = abs(recon_hf).*exp(1i*phase_recon);

         % error check
         if(1)
             relchange   = norm(recon_last(:) - recon_hf(:))/norm(recon_hf(:));
             recon_last  = recon_hf;
             fprintf('Iter = %d, relative change = %4.e \n', ii, relchange);
         end

     end

     recon_hf_complex = recon_hf;

     if (1) % do real diffussion
         %kdata_recon = 1/sqrt(Np*Nf)*fftshift(fft2(ifftshift(recon_hf)));
         kdata_recon =  mrir_fDFT_freqencode(mrir_fDFT_phasencode(recon_hf));
         kdata_recon(:, acq_region,:,:,:) = kdata;

         %recon_hf    = sqrt(Np*Nf)*ifftshift(ifft2(fftshift(kdata_recon)));
         recon_hf    = mrir_iDFT_freqencode(mrir_iDFT_phasencode(kdata_recon));
         recon_hf = recon_hf.*exp(-1i*phase_recon);
     end

 end

 
    
    
    
 