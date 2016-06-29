function pcasl_cbf = cbfmap_base_pCASL(mbfilename,reffilename)
%nii = load_untouch_nii('m_MBpCASLsat_cp_v409_offset90_MB.nii');

% filename = 'p2d_MBpCASLsat_cp_v409_BH3';
% nii = load_untouch_nii([filename '_MB.nii']);
nii = load_untouch_nii(mbfilename);
delta = mean(nii.img(:,:,:,4:2:end),4)-mean(nii.img(:,:,:,5:2:end),4);

%filename = 'ep2d_MBpCASLsat_cp_v409_offset90';
%nii = load_untouch_nii([filename '_ref.nii']);
nii = load_untouch_nii(reffilename);
pcasl_cbf_m0 = nii.img;

mask = zeros(size(delta,1), size(delta,2), size(delta,3));
th = mean(pcasl_cbf_m0(:));
mask(pcasl_cbf_m0 > th) = 1;
%figure, imagesc(tile3d(delta),[-10 50])

tmp = zeros(size(delta,1), size(delta,2), size(delta,3), size(delta,4));
 
w = 1.8;  % laveling time
pld = 1.7; % post-labeling delay
T1b = 1.650;

for ii = 1:5
    tmp(:,:,ii,:) = delta(:,:,ii,:)/(T1b*(exp(-w/T1b) - exp(-(pld+0.05*(ii-1)+w)/T1b)));
    tmp(:,:,ii+5,:) = delta(:,:,ii+5,:)/(T1b*(exp(-w/T1b) - exp(-(pld+0.05*(ii-1)+w)/T1b)));
    tmp(:,:,ii+10,:) = delta(:,:,ii+10,:)/(T1b*(exp(-w/T1b) - exp(-(pld+0.05*(ii-1)+w)/T1b)));
    tmp(:,:,ii+15,:) = delta(:,:,ii+15,:)/(T1b*(exp(-w/T1b) - exp(-(pld+0.05*(ii-1)+w)/T1b)));
    if size(delta,3) == 25 
        tmp(:,:,ii+20,:) = delta(:,:,ii+20,:)/(T1b*(exp(-w/T1b) - exp(-(pld+0.05*(ii-1)+w)/T1b)));
    end
end


ptmp1 = mean(tmp,4);
%figure, imagesc(tile3d(ptmp1),[0 150])

% for ii = 1:size(tmp,3)
%     alpha(:,:,ii) = (ptmp1(:,:,ii)*0.95)./tmp1(:,:,ii);
% end
% figure, imagesc(tile3d(flipdim(alpha(:,:,end:-1:1).*mask(:,:,end:-1:1),2)),[0 1])

alpha = 0.85;
%pcasl_cbf = 60*100*ptmp1./(2*alpha.*double(pcasl_cbf_m0));
pcasl_cbf = 60*100*0.9*tmp./(2*alpha.*double(pcasl_cbf_m0));
pcasl_cbf = pcasl_cbf.*mask;


%figure, imagesc(tile3d(flipdim(pcasl_cbf(:,:,end:-1:1),2)),[0 100]), colormap(gray), axis image, axis off

nii.hdr.dime.dim = [3 size(pcasl_cbf,1) size(pcasl_cbf,2) size(pcasl_cbf,3) 1 1 1 1];
nii.img = pcasl_cbf;
save_untouch_nii(nii,[dirname(mbfilename) '/CBF_pCASL'])
