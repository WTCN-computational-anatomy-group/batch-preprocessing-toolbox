function [tau,lambda] = estimate_model_parameters(V,modality)
N = numel(V);

% Estimate noise precision (tau)
if strcmp(modality,'CT')
    vx0 = vxsize(V{1}{1}.mat);
    img = get_img(V,1,1,modality);
    
    SmoSuf = ComputeSmoSuf_Gauss(img,[],vx0,modality);
    fwhm   = sqrt(4*log(2)*(SmoSuf(2)/SmoSuf(1))/(sum(SmoSuf([4 6 8]))/sum(SmoSuf([3 5 7]))));
    fwhm   = sqrt(max(fwhm^2 - 0.5, 4*log(2)/pi)); 
    ff     = compute_ff_Gauss(fwhm,vx0); 
    tau    = {1/ff};   
else
    tau = cell(1,N);
    for n=1:N
        I    = numel(V{n});
        tau1 = zeros([2 I],'single');
        for i=1:I    
            vx0 = vxsize(V{n}{i}.mat);
            img = get_img(V,n,i,modality);
            
            % From noise in the background
            sd        = spm_noise_estimate(V{n}{i}.fname);       
            tau1(1,i) = 1/(sd*prod(vx0))^2;

            % From estimating FWHM
            SmoSuf    = ComputeSmoSuf_Gauss(img,[],vx0);
            fwhm      = sqrt(4*log(2)*(SmoSuf(2)/SmoSuf(1))/(sum(SmoSuf([4 6 8]))/sum(SmoSuf([3 5 7]))));
            fwhm      = sqrt(max(fwhm^2 - 0.5, 4*log(2)/pi)); 
            ff        = compute_ff_Gauss(fwhm,vx0); 
            tau1(2,i) = 1/ff;  
        end

        % Mean of both sd estimates
        tau{n} = mean(tau1,1);
    end
end

% Set regularisation parameter (lambda)
if strcmp(modality,'CT')
    lambda = 1e-2;
else
    k      = 0.4878; % Approximated from IXI high-res images
    lambda = zeros([N 1],'single');    
    for n=1:N        
        I  = numel(V{n});
        mu = zeros(1,I);
        for i=1:I 
            img   = get_img(V,n,i,modality);
            mu(i) = mean(reshape(img(isfinite(img)),[],1));                   
        end
        sd        = k*mean(mu);
        b         = sqrt(sd.^2./2);
        lambda(n) = 1./b;
    end
end
%==========================================================================

%==========================================================================
function img = get_img(V,n,i,modality)
img       = single(V{n}{i}.private.dat(:,:,:));   
msk       = msk_modality(img,modality);
img(~msk) = NaN;
%==========================================================================