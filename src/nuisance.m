%
% Nuisance inputs: unsmoothed, but processed 4D func, brain_mask, white_test3.nii 
% ROI, rp motion params, fsl motion outliers text files
% output: confounds
%
%
unsmoothed_nii_path = '/home/tseytlin/Data/EMBARC/embarc_CU_TX0001_1R1_mri_fmriraw_20110912/dicom_bold_resting1/wuTX0001_1R1_bold_resting1.nii'
white_mask_roi_path = '/home/tseytlin/Production/software/validation/ROI_EMBARC_scaled2mm/white_test3.nii';
brain_mask_file_path = '/home/tseytlin/Data/EMBARC/embarc_CU_TX0001_1R1_mri_fmriraw_20110912/dicom_bold_resting1/wuTX0001_1R1_bold_resting1_fix_brain_mask.nii'
rp_movement_file_path = '/home/tseytlin/Data/EMBARC/embarc_CU_TX0001_1R1_mri_fmriraw_20110912/dicom_bold_resting1/rp_TX0001_1R1_bold_resting1.txt';
%fd_alone_file_path = 'fd_alone.txt'
%fd_dm_file_path = 'fd_dm.txt'

% white matter ROI
mask_hdr = spm_vol(white_mask_roi_path);
mask_img = spm_read_vols(mask_hdr);

% constants being generated for hi-pass/lo-pass filter
f = (0:180-1)*(0.5/180);
BPF = ((0.01<abs(f))& (abs(f) < 0.08));
BPF2 = ((0.009<abs(f))& (abs(f) < 0.6));

% read in movement parameters
mov_params = importdata(rp_movement_file_path);

% read the output text file
%fd2 = load(fd_dm_file_path);

% if exists load the other text file
%if exist(fd_alone_file_path,'file')   
%	fd = load(fd_alone_file_path);
%end


% we need to use data before smoothing in preprocessing, so I'll need to add another output
% to preprocess() function for unsmoothed data and feed it to this script
uhdr = spm_vol(unsmoothed_nii_path);
uimg = spm_read_vols(uhdr);
NumScans = 180;

% Read in brain mask (BET output from preprocess())
mask3_hdr = spm_vol(brain_mask_file_path);
mask3_img = spm_read_vols(mask3_hdr);


n = 0;
store = 0;
nn = 0;


list = 0;
list2 = 0;
list3 = 0;

for x = 1:79
    for y = 1:95
        for z = 1:68
            if (mask_img(x, y, z) >= 1)
                % create a list of coordinates of white matter
                n= n+1;
                store(1, n) = x;
                store(2, n) = y;
                store(3, n) = z;
                
            elseif mask3_img(x,y,z) == 1
                % create a list of standard deviations per each
                % voxel and voxel locations
                nn = nn+1;
                list(nn, 1) = std(uimg(x, y, z, :));
                list(nn, 2) = x;
                list(nn, 3) = y;
                list(nn, 4) = z;
            end
        end
    end
end

% create an array of timeseries in white matter
array = zeros(n, NumScans);
for q = 1:n
	array(q, :) = uimg(store(1, q), store(2, q), store(3, q), :);
end

% output of above
% a spreadsheet with each column being a voxel in ROI
% and each row is a timeseries


% resort gray matter matrix by sd
list2 = sortrows(list, 1);

% top 2% of voxels with hightest SD
% in list 3
len = length(list2);
lenny = int32(len - (len/50));
for q = lenny:len
    list3(q-lenny+1,1) = list2(q,2);
    list3(q-lenny+1,2) = list2(q,3);
    list3(q-lenny+1,3) = list2(q,4);
end

l = len-lenny;


% timeseries of the top 2% of voxels (new_array)
new_array = zeros(l, NumScans);
for q = 1:l
  new_array(q, :) = uimg(list3(q, 1), list3(q, 2), list3(q, 3), :);
end

% do principal component analysis (PCA)
third_array = 0;
third_array = [array' new_array'];
third_array = third_array';
[M, N] = size(third_array);
mn = mean(third_array, 2);

third_array = third_array - repmat(mn, 1, N);
Y = third_array/(sqrt(N-1));
[u, S, PC] = svd(Y);
S = diag(S);
V = S.*S;
p_confounds = PC(:, 1:5);
m_confounds = mov_params;

% analyze motion params and do FFT (furier) 
for k = 1:6
   spektrum =(fft(m_confounds(:,k)'-mean(m_confounds(:,k)')))/180;
   spektrum1 = BPF2.*spektrum;
   m_confounds(1:180, k) = real(ifft((spektrum1)));
end


% create derivatives (always 6 motion parameters)
for c = 1:6
    m_confounds(1, c+6) = 0;
    for t = 2:NumScans 
    	m_confounds(t, c+6) = m_confounds(t, c) - m_confounds(t-1, c);
    end
end

confounds = [m_confounds, p_confounds, ones(180, 1)];
o_confounds = [m_confounds, p_confounds];



