%--------------------------------------------------------------------------
%% load data
%--------------------------------------------------------------------------
load('data.mat');
[nx,~,nz,nc,nd] = size(kdata_ap);
ny = size(csm,2);
nd = 1;

sgnl_ap = mifft(kdata_ap,1);
sgnl_pa = mifft(kdata_pa,1);
%--------------------------------------------------------------------------
%% image reconstruction loop
%--------------------------------------------------------------------------
idata_ap    = zeros([nx,ny,nz,1,nd]);
idata_pa    = zeros([nx,ny,nz,1,nd]);
FTyC1       = zeros([param.etl*nc,ny*rz]);
FTyC2       = zeros([param.etl*nc,ny*rz]);
lambda1     = 1e-3;
lambda2     = 2.0;
winSize     = [7,7];
    
for id = 1:nd
    fprintf('diffusion direction : %2d \n',id);
    for iz = 1:nz        
        %--------------------------------------------------------------------------
        %% Loop
        %--------------------------------------------------------------------------
        zind = iz;
        for iter = 1:100
            % res in the prev loop
            im_prev = cat(4,idata_ap(:,:,iz,:,id),idata_pa(:,:,iz,:,id));
            for ix = 0.25*nx+1:0.75*nx
                %--------------------------------------------------------------------------
                %% skip background region
                %--------------------------------------------------------------------------
                if max(abs(vec(csm(ix,:,iz,:)))) == 0
                    continue;
                end
                %--------------------------------------------------------------------------
                %% define A and b
                %--------------------------------------------------------------------------
                % rhs - acquired k-space data
                rhs1 = vec(sgnl_ap(ix,:,iz,:,id));
                rhs2 = vec(sgnl_pa(ix,:,iz,:,id));                      
                
                % define A matrix
                for ir = 1:rz
                    ind_rz = (ir-1)*ny + 1 : ir*ny ;
                    % Define B0 effect
                    B0m = vec( b0(ix,:,zind(ir))).';
                    BM1 = exp( sqrt(-1) * 2 * pi * param.time_ap * B0m);
                    BM2 = exp( sqrt(-1) * 2 * pi * param.time_pa * B0m);
                    for ch = 1:nc
                        ind0 = (ch-1)*param.etl + 1 : ch*param.etl ;
                        ind1 = (ch-1)*param.etl + 1 : ch*param.etl ;
                        % Assign Coil Component
                        FTyC1(ind0,ind_rz) = (param.FT1.*BM1) * diag(csm(ix,:,zind(ir),ch));
                        FTyC2(ind1,ind_rz) = (param.FT2.*BM2) * diag(csm(ix,:,zind(ir),ch));
                    end
                end                   
                %--------------------------------------------------------------------------
                %% CG
                %--------------------------------------------------------------------------
                idata_ap(ix,:,:,:,id) = reshape(double(conjgrad(FTyC1'*FTyC1 + lambda1.*eye(ny*rz,ny*rz), FTyC1'*rhs1,vec(idata_ap(ix,:,:,:,id)),100)),[1,ny,rz]);
                idata_pa(ix,:,:,:,id) = reshape(double(conjgrad(FTyC2'*FTyC2 + lambda1.*eye(ny*rz,ny*rz), FTyC2'*rhs2,vec(idata_pa(ix,:,:,:,id)),100)),[1,ny,rz]);
            end
            
            %--------------------------------------------------------------------------
            %% Low Rank Contraint
            %--------------------------------------------------------------------------
            for ir = 1:rz
                img_cat = cat(3,idata_ap(:,:,ind_rz(ir),:,id),idata_pa(:,:,ind_rz(ir),:,id));
                img_cat = cat(3,img_cat,conj(img_cat));
                
                Am      = im2row(mfft2(sq(img_cat)),winSize);
                [U,S,V] = svd(Am,'econ');
                keep    = 1:floor(lambda2*prod(winSize));
                Am      = U(:,keep) * S(keep,keep) * V(:,keep)';
                k_pocs  = Row2im(Am,[nx,ny,4],winSize);
                
                idata_ap(:,:,ind_rz(ir),:,id)  = mifft2(k_pocs(:,:,1));
                idata_pa(:,:,ind_rz(ir),:,id)  = mifft2(k_pocs(:,:,2));
            end            
            
            %--------------------------------------------------------------------------
            %% Check break point
            %--------------------------------------------------------------------------
            rmse_curr = rmse(im_prev,cat(4,idata_ap(:,:,iz,:,id),idata_pa(:,:,iz,:,id)));
%             fprintf('iter: %3d,     update: %5f \n',    iter,   rmse_curr);
            if rmse_curr < 1
                break;
            end            
        end
    end
end
%--------------------------------------------------------------------------
%% Combine two shot
%--------------------------------------------------------------------------
for id = 1:nd
    idata_ap(:,:,:,id) = RealDiffusion_lowRes(idata_ap(:,:,:,id),1,1,0);
    idata_pa(:,:,:,id) = RealDiffusion_lowRes(idata_pa(:,:,:,id),1,1,0);
end
idata_combine = idata_ap + idata_pa;

%--------------------------------------------------------------------------
%% display averaged b0 image & DWI
%--------------------------------------------------------------------------
img_b0 = mean(idata_combine(0.25*end+1:0.75*end,:,floor(nz/2+1),ind_b0),4);
img_b1 = mean(idata_combine(0.25*end+1:0.75*end,:,floor(nz/2+1),ind_b1),4);

figure, imagesc(rot90(abs(img_b0)),[0,1e-3]), colormap gray, axis off image;
figure, imagesc(rot90(abs(img_b1)),[0,5e-4]), colormap gray, axis off image;

%--------------------------------------------------------------------------
%% save images
%--------------------------------------------------------------------------
save_img = @(f,x,m) imwrite(cat(3,x,x,x)./ m,f);
save_img('img_b0.png',  rot90(abs(img_b0)), 1e-3);
save_img('img_dwi.png', rot90(abs(img_b1)), 5e-4);
